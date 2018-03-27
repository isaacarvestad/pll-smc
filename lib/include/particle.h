#ifndef PARTICLE_SMC
#define PARTICLE_SMC

#include <string>
#include <vector>
#include <random>

#include "phylo_forest.h"

/**
   An SMC particle contains a weight as well as an underlying forest state.
 */
class Particle {
  PhyloForest forest;

 public:
  double weight;

  /**
     Constructs a particle with a weight, a vector of sequences and the length
     of each sequence.
   */
  Particle(double weight,
           const std::vector<std::pair<std::string, std::string>> sequences,
           const unsigned int sequence_lengths);


  /**
     Proposes an update to the particle by following the proposal distribution.
   */
  void propose();

  /**
     Returns the current roots of the particles forest.
   */
  std::vector<pll_rnode_s*> get_roots() { return forest.get_roots(); };
};

#endif
