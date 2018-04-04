#include "partition_manager.h"

PartitionManager::PartitionManager(unsigned int tips,
                                   unsigned int clv_buffers,
                                   unsigned int states,
                                   unsigned int sites,
                                   unsigned int rate_matrices,
                                   unsigned int prob_matrices,
                                   unsigned int rate_cats,
                                   unsigned int scale_buffers,
                                   unsigned int attributes)
: tips(tips),
  clv_buffers(clv_buffers),
  states(states),
  sites(sites),
  rate_matrices(rate_matrices),
  prob_matrices(prob_matrices),
  rate_cats(rate_cats),
  scale_buffers(scale_buffers),
  attributes(attributes),
  asc_bias_alloc((attributes & (PLL_ATTRIB_AB_MASK | PLL_ATTRIB_AB_FLAG)) > 0),
  sites_alloc(asc_bias_alloc ? sites + states : sites),
  states_padded(get_states_padded(states, attributes)),
  scaler_size(sites_alloc * rate_cats),

  inner_clv_buffers_count(sites_alloc * states_padded * rate_cats),
  clv_buffers_count(tips + clv_buffers),

  inner_prob_matrices_count(states * states_padded * rate_cats),
  prob_matrices_count(prob_matrices),

  rate_cats_count(rate_cats),
  rate_weights_count(rate_cats),

  inner_subst_params_count((states * states - states) / 2),
  subst_params_count(rate_matrices),

  inner_scale_buffers_count(scaler_size),
  scale_buffers_count(scale_buffers),

  inner_frequencies_count(states_padded),
  frequencies_count(rate_matrices),

  prop_invars_count(rate_matrices),
  invars_count(rate_matrices),

  pattern_weights_count(sites_alloc),
  eigen_decomp_valids_count(rate_matrices),

  inner_eigen_vectors_count(states * states_padded),
  eigen_vectors_count(rate_matrices),

  inner_invar_eigen_vectors_count(states * states_padded),
  invar_eigen_vectors_count(rate_matrices),

  inner_eigen_values_count(states_padded),
  eigen_values_count(rate_matrices),

  inner_tip_chars_count(sites_alloc),
  tip_chars_count(tips),

  char_map_count(PLL_ASCII_SIZE),
  ttlookup_count(1024 * rate_cats),
  tip_map_count(PLL_ASCII_SIZE)
{
  partition_size = compute_partition_size();
  allocate_partition();
}


PartitionManager::PartitionManager(const PartitionManager &original)
  : tips(original.tips),
    clv_buffers(original.clv_buffers),
    states(original.states),
    sites(original.sites),
    rate_matrices(original.rate_matrices),
    prob_matrices(original.prob_matrices),
    rate_cats(original.rate_cats),
    scale_buffers(original.scale_buffers),
    attributes(original.attributes),
    asc_bias_alloc((original.attributes & (PLL_ATTRIB_AB_MASK | PLL_ATTRIB_AB_FLAG)) > 0),
    sites_alloc(original.asc_bias_alloc ? original.sites + original.states : original.sites),
    states_padded(get_states_padded(original.states, original.attributes)),
    scaler_size(original.sites_alloc * original.rate_cats),

    inner_clv_buffers_count(original.sites_alloc * original.states_padded * original.rate_cats),
    clv_buffers_count(original.tips + original.clv_buffers),

    inner_prob_matrices_count(original.states * original.states_padded * original.rate_cats),
    prob_matrices_count(original.prob_matrices),

    rate_cats_count(original.rate_cats),
    rate_weights_count(original.rate_cats),

    inner_subst_params_count((original.states * original.states - original.states) / 2),
    subst_params_count(original.rate_matrices),

    inner_scale_buffers_count(original.scaler_size),
    scale_buffers_count(original.scale_buffers),

    inner_frequencies_count(original.states_padded),
    frequencies_count(original.rate_matrices),

    prop_invars_count(original.rate_matrices),
    invars_count(original.rate_matrices),

    pattern_weights_count(original.sites_alloc),
    eigen_decomp_valids_count(original.rate_matrices),

    inner_eigen_vectors_count(original.states * original.states_padded),
    eigen_vectors_count(original.rate_matrices),

    inner_invar_eigen_vectors_count(original.states * original.states_padded),
    invar_eigen_vectors_count(original.rate_matrices),

    inner_eigen_values_count(original.states_padded),
    eigen_values_count(original.rate_matrices),

    inner_tip_chars_count(original.sites_alloc),
    tip_chars_count(original.tips),

    char_map_count(PLL_ASCII_SIZE),
    ttlookup_count(1024 * original.rate_cats),
    tip_map_count(PLL_ASCII_SIZE),

    partition_size(original.partition_size)
{
  partition = (pll_partition_t *) malloc(sizeof(struct pll_partition) + original.partition_size);
  std::memcpy(partition, original.partition, sizeof(struct pll_partition) + original.partition_size);

  setup_partition_pointers();
}


PartitionManager::~PartitionManager()
{
  free(partition);
}

void PartitionManager::shallow_copy(const PartitionManager &original) {
  std::memcpy(partition, original.partition, sizeof(struct pll_partition) + original.partition_size);

  setup_partition_pointers();
}

unsigned int PartitionManager::get_states_padded(unsigned int states, unsigned int attributes) {
  /* Code taken and modified from pll.c  */
  unsigned int states_padded = states;
#ifdef HAVE_SSE3
  if (attributes & PLL_ATTRIB_ARCH_SSE && PLL_STAT(sse3_present)) {
    states_padded = (states+1) & 0xFFFFFFFE;
  }
#endif
#ifdef HAVE_AVX
  if (attributes & PLL_ATTRIB_ARCH_AVX && PLL_STAT(avx_present)) {
    states_padded = (states+3) & 0xFFFFFFFC;
  }
#endif
#ifdef HAVE_AVX2
  if (attributes & PLL_ATTRIB_ARCH_AVX2 && PLL_STAT(avx2_present)) {
    states_padded = (states+3) & 0xFFFFFFFC;
  }
#endif

  return states_padded;
}

void PartitionManager::allocate_partition() {
  partition = (pll_partition_t*) malloc(sizeof(struct pll_partition) + partition_size);
  assert(partition && "Partition could not be allocated");

  partition->tips = tips;
  partition->clv_buffers = clv_buffers;
  partition->states = states;
  partition->sites = sites;
  partition->rate_matrices = rate_matrices;
  partition->prob_matrices = prob_matrices;
  partition->rate_cats = rate_cats;
  partition->scale_buffers = scale_buffers;
  partition->attributes = attributes;
  partition->alignment = PLL_ALIGNMENT_CPU;
  partition->states_padded = states_padded;
  partition->maxstates = 0;
  partition->asc_bias_alloc = asc_bias_alloc;

  setup_partition_pointers();

  for (int i = 0; i < sites; ++i)
    partition->pattern_weights[i] = 1;
  for (int i = sites; i < sites_alloc; ++i)
    partition->pattern_weights[i] = 0;
  for (int i = 0; i < rate_cats; ++i)
    partition->rate_weights[i] = 1.0/rate_cats;
}

void PartitionManager::setup_partition_pointers() {
  unsigned int delta = sizeof(struct pll_partition);

  partition->clv = (double **) ((char *) partition + delta);
  delta += clv_buffers_count * sizeof(double *);
  for (int i = 0; i < clv_buffers_count; i++) {
    partition->clv[i] = (double *) ((char *) partition + delta + (i * inner_clv_buffers_count));
  }
  delta += clv_buffers_count * inner_clv_buffers_count * sizeof(double);

  partition->pmatrix = (double **) ((char *) partition + delta);
  delta += prob_matrices_count * sizeof(double *);
  for (int i = 0; i < prob_matrices_count; i++) {
    partition->pmatrix[i] = (double *) ((char *) partition + delta + (i * inner_prob_matrices_count));
  }
  delta += prob_matrices_count * inner_prob_matrices_count * sizeof(double);

  partition->rates = (double *) ((char *) partition + delta);
  delta += rate_cats_count * sizeof(double);

  partition->rate_weights = (double *) ((char *) partition + delta);
  delta += rate_weights_count * sizeof(double);

  partition->subst_params = (double **) ((char *) partition + delta);
  delta += subst_params_count * sizeof(double *);
  for (int i = 0; i < subst_params_count; i++) {
    partition->subst_params[i] = (double *) ((char *) partition + delta + (i * inner_subst_params_count));
  }
  delta += subst_params_count * inner_subst_params_count * sizeof(double);

  partition->scale_buffer = (unsigned int **) ((char *) partition + delta);
  delta += scale_buffers * sizeof(unsigned int *);
  for (int i = 0; i < scale_buffers_count; i++) {
    partition->scale_buffer[i] = (unsigned int *) ((char *) partition + delta + (i * inner_scale_buffers_count));
  }
  delta += scale_buffers_count * inner_scale_buffers_count * sizeof(unsigned int);

  partition->frequencies = (double **) ((char *) partition + delta);
  delta += frequencies_count * sizeof(double *);
  for (int i = 0; i < frequencies_count; i++) {
    partition->frequencies[i] = (double *) ((char *) partition + delta + (i * inner_frequencies_count));
  }
  delta += frequencies_count * inner_frequencies_count * sizeof(double);

  partition->prop_invar = (double *) ((char *) partition + delta);
  delta += prop_invars_count * sizeof(double);

  partition->invariant = (int *) ((char *) partition + delta);
  delta += invars_count * sizeof(int);

  partition->pattern_weights = (unsigned int *) ((char *) partition + delta);
  delta += pattern_weights_count * sizeof(unsigned int);

  partition->eigen_decomp_valid = (int *) ((char *) partition + delta);
  delta += eigen_decomp_valids_count * sizeof(int);

  partition->eigenvecs = (double **) ((char *) partition + delta);
  delta += eigen_vectors_count * sizeof(double *);
  for (int i = 0; i < eigen_vectors_count; i++) {
    partition->eigenvecs[i] = (double *) ((char *) partition + delta + (i * inner_eigen_vectors_count));
  }
  delta += eigen_vectors_count * inner_eigen_vectors_count * sizeof(double);

  partition->inv_eigenvecs = (double **) ((char *) partition + delta);
  delta += invar_eigen_vectors_count * sizeof(double *);
  for (int i = 0; i < invar_eigen_vectors_count; i++) {
    partition->inv_eigenvecs[i] = (double *) ((char *) partition + delta + (i * inner_invar_eigen_vectors_count));
  }
  delta += invar_eigen_vectors_count * inner_invar_eigen_vectors_count * sizeof(double);

  partition->eigenvals = (double **) ((char *) partition + delta);
  delta += eigen_values_count * sizeof(double *);
  for (int i = 0; i < eigen_values_count; i++) {
    partition->eigenvals[i] = (double *) ((char *) partition + delta + (i * inner_eigen_values_count));
  }
  delta += eigen_values_count * inner_eigen_vectors_count * sizeof(double);

  partition->tipchars = (unsigned char **) ((char *) partition + delta);
  delta += tip_chars_count * sizeof(unsigned char *);
  for (int i = 0; i < tip_chars_count; i++) {
    partition->tipchars[i] = (unsigned char *) ((char *) partition + delta + (i * inner_tip_chars_count));
  }
  delta += tip_chars_count * inner_tip_chars_count * sizeof(unsigned char);

  partition->charmap = (unsigned char *) ((char *) partition + delta);
  delta += char_map_count * sizeof(unsigned char);

  partition->ttlookup = (double *) ((char *) partition + delta);
  delta += ttlookup_count * sizeof(double);

  partition->tipmap = (unsigned int *) ((char *) partition + delta);
  delta += tip_map_count * sizeof(unsigned int);
}

int PartitionManager::compute_partition_size() {
  return
    clv_buffers_count * sizeof(double *) + clv_buffers_count * inner_clv_buffers_count * sizeof(double) +
    prob_matrices_count * sizeof(double *) + prob_matrices_count * inner_prob_matrices_count * sizeof(double) +
    rate_cats_count * sizeof(double) +
    rate_weights_count * sizeof(double) +
    subst_params_count * sizeof(double *) + subst_params_count * inner_subst_params_count * sizeof(double) +
    scale_buffers_count * sizeof(unsigned int *) + scale_buffers_count * inner_scale_buffers_count * sizeof(unsigned int) +
    frequencies_count * sizeof(double *) + frequencies_count * inner_frequencies_count  * sizeof(double) +
    prop_invars_count * sizeof(double) +
    invars_count * sizeof(int) +
    pattern_weights_count * sizeof(unsigned int) +
    eigen_decomp_valids_count * sizeof(int) +
    eigen_vectors_count * sizeof(double *) + eigen_vectors_count * inner_eigen_vectors_count * sizeof(double) +
    invar_eigen_vectors_count * sizeof(double *) + invar_eigen_vectors_count * inner_invar_eigen_vectors_count * sizeof(double) +
    eigen_values_count * sizeof(double *) + eigen_values_count * inner_eigen_values_count * sizeof(double) +
    tip_chars_count * sizeof(double *) + tip_chars_count * inner_tip_chars_count * sizeof(unsigned char) +
    char_map_count * sizeof(unsigned char) +
    ttlookup_count * sizeof(double) +
    tip_map_count * sizeof(unsigned int);
}
