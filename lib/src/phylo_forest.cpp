#include "phylo_forest.h"

PhyloForest::PhyloForest(
    const std::vector<std::pair<std::string, std::string>> sequences,
    const unsigned int sequence_lengths,
    const pll_partition_t *reference_partition,
    PLLBufferManager *const pll_buffer_manager)
    : reference_partition(reference_partition),
      pll_buffer_manager(pll_buffer_manager), forest_height(0.0) {
  setup_sequences_pll(sequences, sequence_lengths);
}

PhyloForest::PhyloForest(const PhyloForest &original)
    : pll_buffer_manager(original.pll_buffer_manager) {
  reference_partition = original.reference_partition;
  forest_height = original.forest_height;
  roots = original.roots;
}

PhyloForest &PhyloForest::operator=(const PhyloForest &original) {
  if (this == &original)
    return *this;

  reference_partition = original.reference_partition;
  forest_height = original.forest_height;
  roots = original.roots;

  return *this;
}

PhyloForest::~PhyloForest() {}

double compute_ln_likelihood(double *clv, unsigned int *scale_buffer,
                             const pll_partition_t *p) {
  const unsigned int parameter_indices[4] = {0, 0, 0, 0};

  return pll_core_root_loglikelihood(
      p->states, p->sites, p->rate_cats,

      clv, scale_buffer,

      p->frequencies, p->rate_weights, p->pattern_weights, p->prop_invar,
      p->invariant, parameter_indices, NULL, p->attributes);
}

void PhyloForest::setup_sequences_pll(
    std::vector<std::pair<std::string, std::string>> sequences,
    const unsigned int sequence_lengths) {

  const pll_partition_t *p = reference_partition;

  for (unsigned int i = 0; i < sequences.size(); i++) {
    std::string label = sequences[i].first;

    unsigned int sites_alloc =
        p->asc_bias_alloc ? p->sites + p->states : p->sites;

    unsigned int clv_size =
        sites_alloc * p->states_padded * p->rate_cats * sizeof(double);
    unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
                                   ? sites_alloc * p->rate_cats
                                   : sites_alloc;

    std::shared_ptr<PhyloTreeNode> node = std::make_shared<PhyloTreeNode>(
        pll_buffer_manager, nullptr, nullptr, label, forest_height, clv_size,
        scaler_size);
    delete (node->clv);
    delete (node->scale_buffer);
    node->clv = p->clv[i];
    node->scale_buffer = nullptr;
    node->ln_likelihood = compute_ln_likelihood(node->clv, nullptr, p);

    roots.push_back(node);
  }
}

std::shared_ptr<PhyloTreeNode> PhyloForest::connect(int i, int j,
                                                    double height_delta) {
  assert(roots.size() > 1 && "Expected more than one root");
  assert(i != j && "Cannot connect, this would make a loop");
  assert(height_delta >= 0 && "Height change can't be negative");
  assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
         "Index out of bounds");

  const pll_partition_t *p = reference_partition;

  forest_height += height_delta;
  assert(forest_height < 100);

  std::shared_ptr<PhyloTreeNode> child_left = roots[i];
  std::shared_ptr<PhyloTreeNode> child_right = roots[j];

  double left_length = forest_height - child_left->height;
  double right_length = forest_height - child_right->height;

  unsigned int pmatrix_size =
      p->states * p->states_padded * p->rate_cats * sizeof(double);

  std::shared_ptr<PhyloTreeEdge> edge_left = std::make_shared<PhyloTreeEdge>(
      pll_buffer_manager, child_left, left_length, pmatrix_size);

  std::shared_ptr<PhyloTreeEdge> edge_right = std::make_shared<PhyloTreeEdge>(
      pll_buffer_manager, child_right, right_length, pmatrix_size);

  unsigned int sites_alloc =
      p->asc_bias_alloc ? p->sites + p->states : p->sites;

  unsigned int clv_size =
      sites_alloc * p->states_padded * p->rate_cats * sizeof(double);
  unsigned int scaler_size = (p->attributes & PLL_ATTRIB_RATE_SCALERS)
                                 ? sites_alloc * p->rate_cats
                                 : sites_alloc;
  unsigned int scale_buffer_size = scaler_size * sizeof(double);

  std::shared_ptr<PhyloTreeNode> parent = std::make_shared<PhyloTreeNode>(
      pll_buffer_manager, edge_left, edge_right, "", forest_height, clv_size,
      scale_buffer_size);

  const unsigned int matrix_indices[1] = {0};
  const unsigned int param_indices[4] = {0, 0, 0, 0};

  int left_edge_pmatrix_result = pll_core_update_pmatrix(
      &parent->edge_l->pmatrix, p->states, p->rate_cats, p->rates,
      &parent->edge_l->length, matrix_indices, param_indices, p->prop_invar,
      p->eigenvals, p->eigenvecs, p->inv_eigenvecs, 1, p->attributes);
  assert(left_edge_pmatrix_result == PLL_SUCCESS);

  int right_edge_pmatrix_result = pll_core_update_pmatrix(
      &parent->edge_r->pmatrix, p->states, p->rate_cats, p->rates,
      &parent->edge_r->length, matrix_indices, param_indices, p->prop_invar,
      p->eigenvals, p->eigenvecs, p->inv_eigenvecs, 1, p->attributes);
  assert(right_edge_pmatrix_result == PLL_SUCCESS);

  pll_core_update_partial_ii(
      p->states, (p->asc_bias_alloc ? p->sites + p->states : p->sites),
      p->rate_cats, parent->clv, parent->scale_buffer, child_left->clv,
      child_right->clv, edge_left->pmatrix, edge_right->pmatrix,
      child_left->scale_buffer, child_right->scale_buffer, p->attributes);

  parent->ln_likelihood =
      compute_ln_likelihood(parent->clv, parent->scale_buffer, p);

  assert(parent->ln_likelihood <= 0 && "Likelihood can't be more than 100%");

  // Remove children
  remove_roots(i, j);
  // Add new internal node
  roots.push_back(parent);

  return parent;
}

double PhyloForest::likelihood_factor(std::shared_ptr<PhyloTreeNode> root) {
  assert(root->edge_l && root->edge_r && "Root cannot be a leaf");

  std::shared_ptr<PhyloTreeNode> left = root->edge_l->child;
  std::shared_ptr<PhyloTreeNode> right = root->edge_r->child;

  double ln_m = root->ln_likelihood;
  double ln_l = left->ln_likelihood;
  double ln_r = right->ln_likelihood;

  assert(ln_m <= 0 && ln_l <= 0 && ln_r <= 0 &&
         "Likelihood can't be more than 100%");

  return ln_m - (ln_l + ln_r);
}

void PhyloForest::remove_roots(int i, int j) {
  assert(i != j && "Expected different indices");
  assert(i >= 0 && i < roots.size() && j >= 0 && j < roots.size() &&
         "Index out of bounds");

  if (i < j) {
    roots.erase(roots.begin() + j, roots.begin() + j + 1);
    roots.erase(roots.begin() + i, roots.begin() + i + 1);
  } else {
    roots.erase(roots.begin() + i, roots.begin() + i + 1);
    roots.erase(roots.begin() + j, roots.begin() + j + 1);
  }
}
