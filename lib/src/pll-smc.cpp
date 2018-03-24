#include "pll-smc.h"

std::string library_test() {
  pll_partition_t* partition = pll_partition_create(4,2,4,6,1,5,4,2,PLL_ATTRIB_ARCH_AVX);

  return std::to_string(partition->tips);
}
