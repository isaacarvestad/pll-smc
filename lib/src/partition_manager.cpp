#include "partition_manager.h"

PartitionManager::PartitionManager(unsigned int tips,
                                   unsigned int clv_buffers,
                                   unsigned int states,
                                   unsigned int sites,
                                   unsigned int rate_matrices,
                                   unsigned int prob_matrices,
                                   unsigned int rate_cats,
                                   unsigned int scale_buffers,
                                   unsigned int attributes) {
  partition = pll_partition_create(tips,
                                   clv_buffers,
                                   states, sites,
                                   rate_matrices,
                                   prob_matrices,
                                   rate_cats,
                                   scale_buffers,
                                   attributes);

}


PartitionManager::PartitionManager(const PartitionManager &original) {
  partition = pll_partition_create(original.partition->tips,
                                   original.partition->clv_buffers,
                                   original.partition->states,
                                   original.partition->sites,
                                   original.partition->rate_matrices,
                                   original.partition->prob_matrices,
                                   original.partition->rate_cats,
                                   original.partition->scale_buffers,
                                   original.partition->attributes);

  copy_partition(partition, original.partition);
}


PartitionManager::~PartitionManager()
{
  free(partition);
}

void PartitionManager::shallow_copy(const PartitionManager &original) {
  copy_partition(partition, original.partition);
}

void PartitionManager::copy_partition(pll_partition_t *to, const pll_partition_t *from) {
  assert(to->asc_bias_alloc == from->asc_bias_alloc);
  assert(to->sites == from->sites);
  assert(to->states == from->states);

  size_t sites_alloc = from->asc_bias_alloc ? from->sites + from->states : from->sites;

  assert(to->tips == from->tips);
  assert(to->clv_buffers == from->clv_buffers);
  assert(to->states_padded == from->states_padded);
  assert(to->rate_cats == from->rate_cats);

  for (int i = 0; i < from->tips + from->clv_buffers; i++) {
    memcpy(to->clv[i], from->clv[i], sites_alloc * from->states_padded * from->rate_cats * sizeof(double));
  }

  assert(to->prob_matrices == from->prob_matrices);

  size_t displacement = (from->states_padded - from->states) * (from->states_padded) * sizeof(double);
  memcpy(to->pmatrix[0], from->pmatrix[0],
         from->prob_matrices *
         from->states *
         from->states_padded *
         from->rate_cats *
         sizeof(double) +
         displacement);

  assert(to->rate_matrices == from->rate_matrices);

  for (int i = 0; i < from->rate_matrices; i++) {
    memcpy(to->eigenvecs[i], from->eigenvecs[i],
           from->states *
           from->states_padded *
           sizeof(double));

    memcpy(to->inv_eigenvecs[i], from->inv_eigenvecs[i],
           from->states *
           from->states_padded *
           sizeof(double));

    memcpy(to->eigenvals[i], from->eigenvals[i],
           from->states_padded *
           sizeof(double));

    memcpy(to->subst_params[i], from->subst_params[i],
           (from->states * from->states - from->states)/2 *
           sizeof(double));

    memcpy(to->frequencies[i], from->frequencies[i],
           from->states_padded *
           sizeof(double));
  }

  memcpy(to->rates, from->rates,
         from->rate_cats *
         sizeof(double));

  memcpy(to->rate_weights, from->rate_weights,
         from->rate_cats *
         sizeof(double));

  memcpy(to->prop_invar, from->prop_invar,
         from->rate_matrices *
         sizeof(double));

  memcpy(to->pattern_weights, from->pattern_weights,
         sites_alloc *
         sizeof(unsigned int));

  size_t scaler_size = (from->attributes & PLL_ATTRIB_RATE_SCALERS) ? sites_alloc * from->rate_cats : sites_alloc;

  for (int i = 0; i < from->scale_buffers; i++) {
    memcpy(to->scale_buffer[i], from->scale_buffer[i],
           scaler_size *
           sizeof(unsigned int));
  }
}
