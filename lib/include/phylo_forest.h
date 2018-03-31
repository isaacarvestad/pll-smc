#ifndef PHYLO_FOREST_H
#define PHYLO_FOREST_H

#include <string>
#include <vector>

#include <libpll/pll.h>

#include "partition_manager.h"

/**
   A PhyloForest contains a pointer to a 'pll_partition_t', keeps track of pll
   state and performs likelihood calculations.
 */
class PhyloForest {
  unsigned int forest_branch_count = 0;
  unsigned int forest_node_count = 0;
  unsigned int forest_internal_node_count = 0;

  PartitionManager* partition_manager;

  std::vector<pll_rnode_s*> roots;

  /**
     Sets up PLL partition with the correct evolutionary model.
   */
  void setup_pll(const unsigned int leaf_node_count, const unsigned int sequence_lengths);

  /**
     Sets the sequences of the partition.
   */
  void setup_sequences_pll(std::vector<std::pair<std::string, std::string>> sequences,
                           const unsigned int sequence_lengths);


  /**
     Removes two root nodes with index i, j respectively from the root vector.
   */
  void remove_roots(int i, int j);

  /**
     Recursively destroys a tree starting from the root.
   */
  void destroy_tree(pll_rnode_s* root);

 public:

  /**
     Creates a PhyloForest instance using a vector of pairs of (label, sequence)
     and a constant specifying a length of every sequence.
   */
  PhyloForest(const std::vector<std::pair<std::string, std::string>> sequences,
              const unsigned int sequence_lengths);

  /**
     Copy constructor
   */
  PhyloForest(const PhyloForest &original);

  /**
     Remove PLL partition. Does not delete the root vector trees since they may
     be needed by other particles. This will leak memory.
   */
  ~PhyloForest();

  /**
     Connects two root nodes with index i and j respectively in the root vector
     and branch lengths {b1, b2}. Updates the internal PLL partition to reflect
     this change.

     The new internal node created by connecting root nodes i and j is added to
     end of the root vector and it's children are removed from the root vector.

     Returns the new root node.
   */
  pll_rnode_s* connect(int i, int j, double b1, double b2);

  /**
     Computes the likelihood factor (see equation 2.31)
   */
  double likelihood_factor(pll_rnode_s* root);

  /**
     Return the current root nodes of the forrest.
   */
  std::vector<pll_rnode_s*> get_roots() const { return roots; }

  /**
     Returns a pointer to the partition manager.
   */
  const PartitionManager* get_partition_manager() const { return partition_manager; };

  /**
     Number of root nodes in the forrest.
   */
  unsigned int root_count() const { return roots.size(); }
};

#endif
