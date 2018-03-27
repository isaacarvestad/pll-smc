#ifndef LIB_PLL_SMC_H
#define LIB_PLL_SMC_H

#include <string>
#include <vector>
#include <libpll/pll.h>

#include "particle.h"

/**
   Creates a vector with 'count' number of particles, each using the given
   vector of sequences.

   Each particle starts with a weight of 1/'count'.
 */
std::vector<Particle> create_particles(const unsigned int count,
                                       const std::vector<std::pair<std::string, std::string>> sequences);

/**
   Runs the Sequential Monte Carlo algorithm on a vector of input
   particles. Returns the resulting particles.
 */
std::vector<Particle> run_smc(std::vector<Particle> &particles,
                              const unsigned int iterations);

/**
   Resamples the particles based on their weights using multinomial resampling.
 */
void resample(Particle &particle);

/**
   Proposes an update to a partical using the particals propose method.
 */
void propose(Particle &particle);

/**
   Normalizes the weight of a particle where 'sum' is the sum of all partical
   weights.
 */
void normalize_weight(Particle &particle, double sum);

#endif
