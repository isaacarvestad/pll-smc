#include "phylo_forest.h"


PhyloForest::PhyloForest(const std::vector<std::pair<std::string, std::string>> sequences,
                         const unsigned int sequence_lengths) {
  setup_pll(sequences.size(), sequence_lengths);
  setup_sequences_pll(sequences, sequence_lengths);
}

void PhyloForest::setup_pll(const unsigned int leaf_node_count, const unsigned int sequence_lengths) {
  const unsigned int node_count = 2 * leaf_node_count - 1;
  const unsigned int inner_node_count = leaf_node_count - 1;
  const unsigned int branch_count = node_count - 1;

  // Number of states which a nucleotide can have {A,C,G,T}
  const unsigned int nucleotide_states = 4;
  //
  const unsigned int substitution_model_count = 4;

  // One scale buffer per inner node
  const unsigned int scale_buffer_count = inner_node_count;

  // "Discretized category rates from a gamma distribution with alpha shape 1" -
  // libpll
  const unsigned int rate_category_count = 4;
  double rate_categories[4] =
    { 0.13695378267140107,
      0.47675185617665189,
      0.99999999997958422,
      2.38629436117236260
    };

  const double nucleotide_frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };
  double substitution_parameters[6] = { 1, 1, 1, 1, 1, 1 };
  unsigned int parameter_indices[4] = { 0, 0, 0, 0 };

  partition = pll_partition_create(leaf_node_count,
                                   inner_node_count,
                                   nucleotide_states,
                                   sequence_lengths,
                                   substitution_model_count,
                                   // One probability matrix per branch
                                   branch_count,
                                   rate_category_count,
                                   scale_buffer_count,
                                   PLL_ATTRIB_ARCH_AVX);

  pll_set_frequencies(partition, 0, nucleotide_frequencies);
  pll_set_subst_params(partition, 0, substitution_parameters);
  pll_set_category_rates(partition, rate_categories);
}

void PhyloForest::setup_sequences_pll(std::vector<std::pair<std::string, std::string>> sequences,
                                      const unsigned int sequence_lengths) {

  for (unsigned int i = 0; i < sequences.size(); i++) {
    std::string label = sequences[i].first;
    std::string sequence = sequences[i].second;

    pll_set_tip_states(partition, i, pll_map_nt, sequence.data());

    pll_rnode_s * node = new pll_rnode_s
      { .label = strdup(label.c_str()),
        .length = 0.0f,
        .node_index = forest_node_count,
        .clv_index = forest_node_count,
        .scaler_index = PLL_SCALE_BUFFER_NONE,
        .pmatrix_index = 0,
        .left = nullptr,
        .right = nullptr,
        .parent = nullptr
      };

    forest_node_count++;
    roots.push_back(node);
  }

}


pll_rnode_s* PhyloForest::connect(int i, int j, double b1, double b2) {
  assert(roots.size() > 1 && "Expected more than one root");
  assert(i != j && "Cannot connect, this would make a loop");
  assert(b1 >= 0 && b2 >= 0 && "Branch length can't be negative");
  assert(i >= 0 && i < roots.size() &&
         j >= 0 && j < roots.size() && "Index out of bounds");

  pll_rnode_s * left = roots[i];
  pll_rnode_s * right = roots[j];

  left->length = b1;
  left->pmatrix_index = forest_branch_count;

  right->length = b2;
  right->pmatrix_index = forest_branch_count + 1;

  forest_branch_count += 2;

  pll_rnode_s * combined = new pll_rnode_s
    { .label = (char*) "",
      .length = 0.0f,
      .node_index = forest_node_count,
      .clv_index = forest_node_count,
      .scaler_index = (int) forest_internal_node_count,
      .pmatrix_index = 0,
      .left = left,
      .right = right,
      .parent = nullptr
    };

  forest_node_count++;
  forest_internal_node_count++;

  // Remove children
  remove_roots(i, j);
  // Add new internal node
  roots.push_back(combined);

  pll_operation_t operation =
    { .parent_clv_index = combined->clv_index,
      .parent_scaler_index = combined->scaler_index,
      .child1_clv_index = combined->left->clv_index,
      .child2_clv_index = combined->right->clv_index,
      .child1_matrix_index = combined->left->pmatrix_index,
      .child2_matrix_index = combined->right->pmatrix_index,
      .child1_scaler_index = combined->left->scaler_index,
      .child2_scaler_index = combined->right->scaler_index,
    };
  pll_operation_t operations[1] = { operation };
  unsigned int matrix_indices[2] = { operation.child1_matrix_index,
                                     operation.child2_matrix_index };
  double branch_lengths[2] = { combined->left->length,
                               combined->right->length };

  unsigned int parameter_indices[4] = { 0, 0, 0, 0 };

  pll_update_prob_matrices(partition,
                           parameter_indices,
                           matrix_indices,
                           branch_lengths,
                           2);

  pll_update_partials(partition, operations, 1);

  return combined;
}

/**
   TODO: compute q and delta
 */
double PhyloForest::likelihood_factor(pll_rnode_s* root) {
  assert(root->left && root->right && "Root cannot be a leaf");

  unsigned int parameter_indices[4] = { 0, 0, 0, 0 };

  double l_merged = pll_compute_root_loglikelihood(partition, root->clv_index, 0, parameter_indices, NULL);
  double l_left = pll_compute_root_loglikelihood(partition, root->left->clv_index, 0, parameter_indices, NULL);
  double l_right = pll_compute_root_loglikelihood(partition, root->right->clv_index, 0, parameter_indices, NULL);

  return l_merged / (l_left * l_right);
}

void PhyloForest::remove_roots(int i, int j) {
  assert(i != j && "Expected different indices");
  assert(i >= 0 && i < roots.size() &&
         j >= 0 && j < roots.size() && "Index out of bounds");

  if (i < j) {
    roots.erase(roots.begin() + j, roots.begin() + j + 1);
    roots.erase(roots.begin() + i, roots.begin() + i + 1);
  } else {
    roots.erase(roots.begin() + i, roots.begin() + i + 1);
    roots.erase(roots.begin() + j, roots.begin() + j + 1);
  }
}
