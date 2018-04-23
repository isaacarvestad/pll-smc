#ifndef LIB_PLL_SMC_H
#define LIB_PLL_SMC_H

#include <float.h>
#include <libpll/pll.h>

#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "particle.h"
#include "phylo_tree.h"

/**
   Runs the Sequential Monte Carlo algorithm with a number of particles. Returns
   the resulting particles.
 */
std::vector<Particle *>
run_smc(const unsigned int particle_count,
        const std::vector<std::pair<std::string, std::string>> sequences);

/**
   Resamples the particles based on their weights using multinomial resampling.
 */
void resample(std::vector<Particle *> &particles, const unsigned int iteration);

/**
   Proposes an update to a partical using the particals proposal method.
 */
void propose(std::vector<Particle *> &particles, const unsigned int iteration);

/**
   Normalizes the weight of the particle.
 */
void normalize_weights(std::vector<Particle *> &particles,
                       const unsigned int iteration);

#endif
