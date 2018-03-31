#ifndef LIB_PLL_SMC_PARTITION_MANAGER_H
#define LIB_PLL_SMC_PARTITION_MANAGER_H

#include <cstring> // Include for memcpy

#include <libpll/pll.h>

class PartitionManager {
  unsigned int partition_size;
  pll_partition_t* partition;

  // Partition variables
  const unsigned int tips;
  const unsigned int clv_buffers;
  const unsigned int states;
  const unsigned int sites;
  const unsigned int rate_matrices;
  const unsigned int prob_matrices;
  const unsigned int rate_cats;
  const unsigned int scale_buffers;
  const unsigned int attributes;

  const unsigned int asc_bias_alloc;
  const unsigned int sites_alloc;
  const unsigned int states_padded;
  const unsigned int scaler_size;

  // Partition pointer counts
  const unsigned int inner_clv_buffers_count;
  const unsigned int clv_buffers_count;

  const unsigned int inner_prob_matrices_count;
  const unsigned int prob_matrices_count;

  const unsigned int rate_cats_count;
  const unsigned int rate_weights_count;

  const unsigned int inner_subst_params_count;
  const unsigned int subst_params_count;

  const unsigned int inner_scale_buffers_count;
  const unsigned int scale_buffers_count;

  const unsigned int inner_frequencies_count;
  const unsigned int frequencies_count;

  const unsigned int prop_invars_count;
  const unsigned int invars_count;

  const unsigned int pattern_weights_count;
  const unsigned int eigen_decomp_valids_count;

  const unsigned int inner_eigen_vectors_count;
  const unsigned int eigen_vectors_count;

  const unsigned int inner_invar_eigen_vectors_count;
  const unsigned int invar_eigen_vectors_count;

  const unsigned int inner_eigen_values_count;
  const unsigned int eigen_values_count;

  const unsigned int inner_tip_chars_count;
  const unsigned int tip_chars_count;

  const unsigned int char_map_count;
  const unsigned int ttlookup_count;
  const unsigned int tip_map_count;


  /**
     Depending on the CPU architecture PLL might have to pad the CLV states.
  */
  unsigned int get_states_padded(unsigned int states, unsigned int attributes);

  void allocate_partition();

  /**
     Set the partition pointers to the correct memory locatsions.
   */
  void setup_partition_pointers();

  /**
     Compute the size of the PLL partition.
   */
  int compute_partition_size();
 public:

  PartitionManager(unsigned int tips,
                   unsigned int clv_buffers,
                   unsigned int states,
                   unsigned int sites,
                   unsigned int rate_matrices,
                   unsigned int prob_matrices,
                   unsigned int rate_cats,
                   unsigned int scale_buffers,
                   unsigned int attributes);

  PartitionManager(const PartitionManager &original);

  ~PartitionManager();

  pll_partition_t * get_partition() const { return partition; };
};

#endif
