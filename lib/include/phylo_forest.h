#ifndef PHYLO_FOREST_H
#define PHYLO_FOREST_H

#include <string>
#include <vector>

#include <libpll/pll.h>

#include "partition_manager.h"

/**
   Similar to a pll_rnode struct but with a height field which is the sum of the
   branch distances from the node to any leaf which it covers.
 */
struct phylo_tree_node {
  std::string label;
  double length;
  double height;
  unsigned int node_index;
  unsigned int clv_index;
  int scaler_index;
  unsigned int pmatrix_index;
  phylo_tree_node* left;
  phylo_tree_node* right;
  phylo_tree_node* parent;
};

/**
   A PhyloForest contains a pointer to a 'pll_partition_t', keeps track of pll
   state and performs likelihood calculations.
 */
class PhyloForest {
  unsigned int forest_branch_count = 0;
  unsigned int forest_node_count = 0;
  unsigned int forest_internal_node_count = 0;

  PartitionManager* partition_manager;

  double forest_height = 0.0;
  std::vector<phylo_tree_node*> roots;

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
  void destroy_tree(phylo_tree_node* root);

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
     Connects two root nodes with index 'i' and 'j' respectively in the root
     vector. Updates the internal PLL partition to reflect this change. Branch
     lengths are calculated from the height of the new internal node given by
     the parameter 'height'.

     The new internal node created by connecting root nodes i and j is added to
     end of the root vector and it's children are removed from the root vector.

     Returns the new root node.
   */
  phylo_tree_node* connect(int i, int j, double height);

  /**
     Computes the likelihood factor (see equation 2.31)
   */
  double likelihood_factor(phylo_tree_node* root);

  /**
     Return the current root nodes of the forrest.
   */
  std::vector<phylo_tree_node*> get_roots() const { return roots; }

  /**
     Returns a pointer to the partition manager.
   */
  const PartitionManager* get_partition_manager() const { return partition_manager; };

  /**
     Number of root nodes in the forrest.
   */
  unsigned int root_count() const { return roots.size(); }

  /**
     Maximum tree height in the forest.
   */
  double get_forest_height() const { return forest_height; };
};

#endif
