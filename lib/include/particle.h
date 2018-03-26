#ifndef PARTICLE_SMC
#define PARTICLE_SMC

#include <string>
#include <vector>

#include "phylo_forest.h"

/**
   An SMC particle contains a weight as well as an underlying forest state.
 */
class Particle {
  double weight;
  PhyloForest forest;

 public:
  /**
     Constructs a particle with a weight and vector of sequences.
   */
  Particle(double weight, const std::vector<std::string> sequences);
};

#endif
