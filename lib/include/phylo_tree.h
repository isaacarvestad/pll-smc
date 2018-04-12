#ifndef LIB_PLL_SMC_PHYLO_TREE_H
#define LIB_PLL_SMC_PHYLO_TREE_H

#include <string>
#include <vector>

struct phylo_tree_node;

struct phylo_tree_edge {
  phylo_tree_node* child;

  double length;

  double* pmatrix;
};

struct phylo_tree_node {
  std::string label;
  double height;

  double* clv;
  unsigned int* scale_buffer;

  phylo_tree_edge* child_edge_l;
  phylo_tree_edge* child_edge_r;
};

#endif
