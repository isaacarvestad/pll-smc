#ifndef LIB_PLL_SMC_PLL_BUFFER_MANAGER_H
#define LIB_PLL_SMC_PLL_BUFFER_MANAGER_H

#include <stack>

/**
   A struct which keeps track of allocated but unused PLL data buffers.
 */
struct PLLBufferManager {
  std::stack<double *> clv_buffer;
  std::stack<double *> pmatrix_buffer;
  std::stack<unsigned int *> scale_buffer_buffer;
};

#endif
