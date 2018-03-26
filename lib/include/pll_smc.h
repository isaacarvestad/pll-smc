#ifndef LIB_PLL_SMC_H
#define LIB_PLL_SMC_H

#include <string>
#include <vector>
#include <libpll/pll.h>

#include "particle.h"

/**
   Creates a vector of 'count' number of particles, each using the given vector
   of sequences.

   Each particle starts with a weight of 1/'count'.
 */
std::vector<Particle> create_particles(const unsigned int count,
                                       std::vector<std::string> sequences);

/**
   Runs the Sequential Monte Carlo algorithm on a vector of input
   particles. Returns the resulting particles.
 */
std::vector<Particle> run_smc(std::vector<Particle> particles,
                              const unsigned int iterations);

void resample(Particle &particle);
void propose(Particle &particle);
void weight(Particle &particle);

#endif
