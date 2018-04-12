#include "phylo_forest.h"

PhyloForest::PhyloForest(const std::vector<std::pair<std::string, std::string>> sequences,
                         const unsigned int sequence_lengths,
                         const pll_partition_t* reference_partition)
  : reference_partition(reference_partition)
  , forest_height(0.0) {
  setup_sequences_pll(sequences, sequence_lengths);
}

PhyloForest::PhyloForest(const PhyloForest& original) {
  reference_partition = original.reference_partition;
  forest_height = original.forest_height;
  roots = original.roots;
}

PhyloForest& PhyloForest::operator=(const PhyloForest& original) {
  if (this == &original) return *this;

  reference_partition = original.reference_partition;
  forest_height = original.forest_height;
  roots = original.roots;

  return *this;
}

PhyloForest::~PhyloForest() {
}

void PhyloForest::setup_sequences_pll(std::vector<std::pair<std::string, std::string>> sequences,
                                      const unsigned int sequence_lengths) {

  for (unsigned int i = 0; i < sequences.size(); i++) {
    std::string label = sequences[i].first;

    phylo_tree_node * node = new phylo_tree_node
      { .label = label,
        .height = forest_height,

        .clv = reference_partition->clv[i],
        .scale_buffer = nullptr,

        .child_edge_l = nullptr,
        .child_edge_r = nullptr
      };
    roots.push_back(node);
  }
}


phylo_tree_node* PhyloForest::connect(int i, int j, double height_delta) {
  assert(roots.size() > 1 && "Expected more than one root");
  assert(i != j && "Cannot connect, this would make a loop");
  assert(height_delta >= 0 && "Height change can't be negative");
  assert(i >= 0 && i < roots.size() &&
         j >= 0 && j < roots.size() && "Index out of bounds");

  const pll_partition_t* p = reference_partition;

  forest_height += height_delta;
  assert(forest_height < 100);

  phylo_tree_node* child_left = roots[i];
  phylo_tree_node* child_right = roots[j];

  phylo_tree_edge* edge_left = new phylo_tree_edge
    { .child = child_left,
      .length = forest_height - child_left->height,
      .pmatrix = (double *) std::malloc(p->states *
                                        p->states_padded *
                                        p->rate_cats *
                                        sizeof(double)
                                        ) // TODO: double check 'displacement' is not needed
    };
  phylo_tree_edge* edge_right = new phylo_tree_edge
    { .child = child_right,
      .length = forest_height - child_right->height,
      .pmatrix = (double *) std::malloc(p->states *
                                        p->states_padded *
                                        p->rate_cats *
                                        sizeof(double)
                                        ) // TODO: double check 'displacement' is not needed
    };

  unsigned int sites_alloc = p->asc_bias_alloc ?
    p->sites + p->states :
    p->sites;

  unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS) ?
    sites_alloc * p->rate_cats :
    sites_alloc;

  phylo_tree_node* parent = new phylo_tree_node
    { .label = "",
      .height = forest_height,
      .clv = (double *) std::malloc(sites_alloc *
                                    p->states_padded *
                                    p->rate_cats *
                                    sizeof(double)
                                    ),
      .scale_buffer = (unsigned int *) std::malloc(scaler_size * sizeof(unsigned int)),
      .child_edge_l = edge_left,
      .child_edge_r = edge_right
    };

  const unsigned int matrix_indices[1] = { 0 };
  const unsigned int param_indices[4] = { 0, 0, 0, 0 };

  int left_edge_pmatrix_result =
    pll_core_update_pmatrix(&parent->child_edge_l->pmatrix,
                            p->states,
                            p->rate_cats,
                            p->rates,
                            &parent->child_edge_l->length,
                            matrix_indices,
                            param_indices,
                            p->prop_invar,
                            p->eigenvals,
                            p->eigenvecs,
                            p->inv_eigenvecs,
                            1,
                            p->attributes
                            );
  assert(left_edge_pmatrix_result == PLL_SUCCESS);

  int right_edge_pmatrix_result =
    pll_core_update_pmatrix(&parent->child_edge_r->pmatrix,
                            p->states,
                            p->rate_cats,
                            p->rates,
                            &parent->child_edge_r->length,
                            matrix_indices,
                            param_indices,
                            p->prop_invar,
                            p->eigenvals,
                            p->eigenvecs,
                            p->inv_eigenvecs,
                            1,
                            p->attributes
                            );
  assert(right_edge_pmatrix_result == PLL_SUCCESS);

  pll_core_update_partial_ii(p->states,
                             (p->asc_bias_alloc ? p->sites + p->states : p->sites),
                             p->rate_cats,
                             parent->clv,
                             parent->scale_buffer,
                             child_left->clv,
                             child_right->clv,
                             edge_left->pmatrix,
                             edge_right->pmatrix,
                             child_left->scale_buffer,
                             child_right->scale_buffer,
                             p->attributes
                             );

  // Remove children
  remove_roots(i, j);
  // Add new internal node
  roots.push_back(parent);

  return parent;
}

double compute_ln_likelihood(double* clv, unsigned int* scale_buffer, const pll_partition_t* p) {
  const unsigned int parameter_indices[4] = { 0, 0, 0, 0 };

  return pll_core_root_loglikelihood(p->states,
                                     p->sites,
                                     p->rate_cats,

                                     clv,
                                     scale_buffer,

                                     p->frequencies,
                                     p->rate_weights,
                                     p->pattern_weights,
                                     p->prop_invar,
                                     p->invariant,
                                     parameter_indices,
                                     NULL,
                                     p->attributes);
}

double PhyloForest::likelihood_factor(phylo_tree_node* root) {
  assert(root->child_edge_l && root->child_edge_r && "Root cannot be a leaf");

  phylo_tree_node* left = root->child_edge_l->child;
  phylo_tree_node* right = root->child_edge_r->child;

  double ln_m = compute_ln_likelihood(root->clv, root->scale_buffer, reference_partition);
  double ln_l = compute_ln_likelihood(left->clv, left->scale_buffer, reference_partition);
  double ln_r = compute_ln_likelihood(right->clv, right->scale_buffer, reference_partition);

  assert(ln_m <= 0 && ln_l && ln_r && "Likelihood can't be more than 100%");

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
  phylo_tree_node* left = root->child_edge_l->child;
  phylo_tree_node* right = root->child_edge_r->child;

  delete(root->child_edge_l);
  delete(root->child_edge_r);
  delete(root->clv);
  delete(root->scale_buffer);
  delete(root);
  if (left) destroy_tree(left);
  if (right) destroy_tree(right);
}
