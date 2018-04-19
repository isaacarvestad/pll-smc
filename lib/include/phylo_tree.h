#ifndef LIB_PLL_SMC_PHYLO_TREE_H
#define LIB_PLL_SMC_PHYLO_TREE_H

#include <memory>
#include <string>
#include <vector>

#include "pll_buffer_manager.h"

class PhyloTreeNode;

/**
   An edge in a PhyloTree points to a child node and keeps track of a pmatrix
   buffer.
 */
class PhyloTreeEdge {

  /**
     The PLLBufferManager to use when recyling the pmatrix.
   */
  PLLBufferManager *manager;

public:
  /**
     Constructs an edge. Attempts to re-use an existing pmatrix buffer if there
     is one available in 'manager'.
   */
  PhyloTreeEdge(PLLBufferManager *manager, std::shared_ptr<PhyloTreeNode> child,
                double length, unsigned int pmatrix_size);

  /**
     Destroys the edge and adds the pmatrix buffer to the 'manager'.
   */
  ~PhyloTreeEdge();

  std::shared_ptr<PhyloTreeNode> child;

  double length;
  double *pmatrix;
};

/**
   A node in a PhyloTree is either a leaf or points to two edges. A node also
   keeps track of a PLL clv buffer and a scale buffer.
 */
class PhyloTreeNode {

  /**
     The PLLBufferManager to use when recycling the pmatrix.
   */
  PLLBufferManager *manager;

public:
  /**
     Constructs a node. Attempts to re-use an exiting clv buffer and scale
     buffer from the 'manager'.
   */
  PhyloTreeNode(PLLBufferManager *manager,
                std::shared_ptr<PhyloTreeEdge> edge_l,
                std::shared_ptr<PhyloTreeEdge> edge_r, std::string label,
                double height, unsigned int clv_size,
                unsigned int scale_buffer_size);

  /**
     Destroys the node and adds the clv buffer and scale buffer to the
     'manager'.
   */
  ~PhyloTreeNode();

  std::shared_ptr<PhyloTreeEdge> edge_l;
  std::shared_ptr<PhyloTreeEdge> edge_r;

  std::string label;
  double height;
  double ln_likelihood;

  double *clv;
  unsigned int *scale_buffer;
};

#endif
