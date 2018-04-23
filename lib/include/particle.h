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
  PhyloForest* forest;
  std::mt19937 mt_generator;

 public:
  double weight;
  double normalized_weight;

  /**
     Constructs a particle with a weight, a vector of sequences and the length
     of each sequence.
   */
  Particle(double weight,
           const std::vector<std::pair<std::string, std::string>> sequences,
           const unsigned int sequence_lengths);

  /**
     Copy constructor. Copies the particles forest but creates a new random
     generator.
   */
  Particle(const Particle &original);

  /**
     Frees the particle.
   */
  ~Particle();

  /**
     Copies a particle with a shallow copy of the particle's forest.
   */
  void shallow_copy(const Particle &original);

  /**
     Proposes an update to the particle by following the proposal distribution.
   */
  void propose();

  /**
     Returns the current roots of the particles forest.
   */
  std::vector<phylo_tree_node*> get_roots() const { return forest->get_roots(); };

  /**
     Returns the particles forest.
   */
  PhyloForest* get_forest() const { return forest; };
};

#endif
