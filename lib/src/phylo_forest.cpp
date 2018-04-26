#include "phylo_forest.h"

PhyloForest::PhyloForest(const std::vector<std::pair<std::string, std::string>> sequences,
                         const unsigned int sequence_lengths) {
  setup_pll(sequences.size(), sequence_lengths);
  setup_sequences_pll(sequences, sequence_lengths);
}

PhyloForest::PhyloForest(const PhyloForest &original) {
  forest_branch_count = original.forest_branch_count;
  forest_node_count = original.forest_node_count;
  forest_internal_node_count = original.forest_internal_node_count;

  partition_manager = new PartitionManager(*original.get_partition_manager());

  forest_height = original.forest_height;

  for (auto &r : original.get_roots()) {
    phylo_tree_node* root = new phylo_tree_node;
    *root = *r;
    roots.push_back(root);
  }
}

PhyloForest::~PhyloForest() {
  delete(partition_manager);
}

void PhyloForest::shallow_copy(const PhyloForest &original) {
  forest_branch_count = original.get_forest_branch_count();
  forest_node_count = original.get_forest_node_count();
  forest_internal_node_count = original.get_forest_internal_node_count();

  partition_manager->shallow_copy(*original.get_partition_manager());

  forest_height = original.get_forest_height();

  roots.clear();
  for (auto &r : original.get_roots()) {
    phylo_tree_node* root = new phylo_tree_node(*r);
    roots.push_back(root);
  }
}

void PhyloForest::setup_pll(const unsigned int leaf_node_count, const unsigned int sequence_lengths) {
  const unsigned int node_count = 2 * leaf_node_count - 1;
  const unsigned int inner_node_count = leaf_node_count - 1;
  const unsigned int branch_count = node_count - 1;

  // Number of states which a nucleotide can have {A,C,G,T}
  const unsigned int nucleotide_states = 4;
  //
  const unsigned int substitution_model_count = 1;

  // One scale buffer per inner node
  const unsigned int scale_buffer_count = inner_node_count;

  // From libpll/wiki "Discretized category rates from a gamma distribution with alpha shape 1"
  const unsigned int rate_category_count = 4;
  double rate_categories[4] =
    { 0.13695378267140107,
      0.47675185617665189,
      0.99999999997958422,
      2.38629436117236260
    };

  const double nucleotide_frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
  double substitution_parameters[6] = { 1, 1, 1, 1, 1, 1 };

  partition_manager = new PartitionManager(leaf_node_count,
                                           inner_node_count,
                                           nucleotide_states,
                                           sequence_lengths,
                                           substitution_model_count,
                                           branch_count,
                                           rate_category_count,
                                           scale_buffer_count,
                                           PLL_ATTRIB_ARCH_SSE);
  assert(partition_manager);

  pll_set_frequencies(partition_manager->get_partition(), 0, nucleotide_frequencies);
  pll_set_category_rates(partition_manager->get_partition(), rate_categories);
  pll_set_subst_params(partition_manager->get_partition(), 0, substitution_parameters);
}

void PhyloForest::setup_sequences_pll(std::vector<std::pair<std::string, std::string>> sequences,
                                      const unsigned int sequence_lengths) {

  for (unsigned int i = 0; i < sequences.size(); i++) {
    std::string label = sequences[i].first;
    std::string sequence = sequences[i].second;

    pll_set_tip_states(partition_manager->get_partition(), i, pll_map_nt, sequence.data());

    phylo_tree_node * node = new phylo_tree_node
      { .label = label,
        .length = 0.0,
        .height = forest_height,
        .node_index = forest_node_count,
        .clv_index = forest_node_count,
        .scaler_index = PLL_SCALE_BUFFER_NONE,
        .pmatrix_index = 0,
        .left = nullptr,
        .right = nullptr,
        .parent = nullptr,
        .ln_likelihood = 0.0
      };

    unsigned int parameter_indices[4] = { 0, 0, 0, 0 };
    node->ln_likelihood = pll_compute_root_loglikelihood(partition_manager->get_partition(),
                                                         node->clv_index, node->scaler_index,
                                                         parameter_indices, NULL);

    forest_node_count++;
    roots.push_back(node);
  }
}


phylo_tree_node* PhyloForest::connect(int i, int j, double height_delta) {
  assert(roots.size() > 1 && "Expected more than one root");
  assert(i != j && "Cannot connect, this would make a loop");
  assert(height_delta >= 0 && "Height change can't be negative");
  assert(i >= 0 && i < roots.size() &&
         j >= 0 && j < roots.size() && "Index out of bounds");

  phylo_tree_node * left = roots[i];
  phylo_tree_node * right = roots[j];

  forest_height += height_delta;

  left->length = forest_height - left->height;
  assert(left->length >= 0 && "Branch length can't be negative");
  left->pmatrix_index = forest_branch_count;

  right->length = forest_height - right->height;
  assert(right->length >= 0 && "Branch length can't be negative");
  right->pmatrix_index = forest_branch_count + 1;

  forest_branch_count += 2;

  phylo_tree_node * combined = new phylo_tree_node
    { .label = "",
      .length = 0.0,
      .height = forest_height,
      .node_index = forest_node_count,
      .clv_index = forest_node_count,
      .scaler_index = (int) forest_internal_node_count,
      .pmatrix_index = 0,
      .left = left,
      .right = right,
      .parent = nullptr,
      .ln_likelihood = 0.0
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

  pll_update_prob_matrices(partition_manager->get_partition(),
                           parameter_indices,
                           matrix_indices,
                           branch_lengths,
                           2);

  pll_update_partials(partition_manager->get_partition(), operations, 1);

  combined->ln_likelihood = pll_compute_root_loglikelihood(partition_manager->get_partition(),
                                                           combined->clv_index, combined->scaler_index,
                                                           parameter_indices, NULL);

  return combined;
}

double PhyloForest::likelihood_factor(phylo_tree_node* root) {
  assert(root->left && root->right && "Root cannot be a leaf");

  double ln_m = root->ln_likelihood;
  double ln_l = root->left->ln_likelihood;
  double ln_r = root->right->ln_likelihood;

  assert(ln_m <= 0 && ln_l <= 0 && ln_r <= 0 && "Likelihood can't be more than 100%");

  return ln_m - (ln_l + ln_r);
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

void PhyloForest::destroy_tree(phylo_tree_node* root) {
  phylo_tree_node* left = root->left;
  phylo_tree_node* right = root->right;

  delete(root);
  if (left) destroy_tree(left);
  if (right) destroy_tree(right);
}
