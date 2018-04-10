#ifndef LIB_PLL_SMC_PARTITION_MANAGER_H
#define LIB_PLL_SMC_PARTITION_MANAGER_H

#include <cstring> // Include for memcpy

#include <libpll/pll.h>

class PartitionManager {
  unsigned int partition_size;
  pll_partition_t* partition;

  void copy_partition(pll_partition_t *to, const pll_partition_t *from);

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

  /**
     Copies a partition manager without allocating new memory for the underlying
     PLL partition.
   */
  void shallow_copy(const PartitionManager &original);

  pll_partition_t * get_partition() const { return partition; };
};

#endif
