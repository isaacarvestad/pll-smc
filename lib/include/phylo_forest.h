#ifndef PHYLO_FOREST_H
#define PHYLO_FOREST_H

#include <string>
#include <vector>
#include <memory>

#include <libpll/pll.h>

#include "pll_buffer_manager.h"
#include "phylo_tree.h"

class PhyloForest {
  const pll_partition_t* reference_partition;
  PLLBufferManager* const pll_buffer_manager;

  double forest_height;
  std::vector<std::shared_ptr<PhyloTreeNode>> roots;

  /**
     Sets the sequences of the partition.
   */
  void setup_sequences_pll(std::vector<std::pair<std::string, std::string>> sequences,
                           const unsigned int sequence_lengths);


  /**
     Removes two root nodes with index i, j respectively from the root vector.
   */
  void remove_roots(int i, int j);

 public:

  /**
     Creates a PhyloForest instance using a vector of pairs of (label, sequence)
     and a constant specifying a length of every sequence.
   */
  PhyloForest(const std::vector<std::pair<std::string, std::string>> sequences,
              const unsigned int sequence_lengths,
              const pll_partition_t* reference_partition,
              PLLBufferManager* const pll_buffer_manager);

  /**
     Copy constructor
   */
  PhyloForest(const PhyloForest &original);

  /**
     Copy reference partition, forest height and root vector.
   */
  PhyloForest& operator=(const PhyloForest& original);

  /**
     Remove PLL partition. Does not delete the root vector trees since they may
     be needed by other particles. This will leak memory.
   */
  ~PhyloForest();

  /**
     Copies a forest with a shallow copy of the forest's partition manager.
   */
  void shallow_copy(const PhyloForest &original);

  /**
     Connects two root nodes with index 'i' and 'j' respectively in the root
     vector. Updates the internal PLL partition to reflect this change. Branch
     lengths are calculated from the height of the new internal node given by
     the parameter 'height'.

     The new internal node created by connecting root nodes i and j is added to
     end of the root vector and it's children are removed from the root vector.

     Returns the new root node.
   */
  std::shared_ptr<PhyloTreeNode> connect(int i, int j, double height);

  /**
     Computes the likelihood factor (see equation 2.31)
   */
  double likelihood_factor(std::shared_ptr<PhyloTreeNode> root);

  /**
     Return the current root nodes of the forrest.
   */
  std::vector<std::shared_ptr<PhyloTreeNode>> get_roots() const { return roots; }

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
