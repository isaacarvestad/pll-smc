#include "pll_smc.h"

std::vector<Particle> create_particles(const unsigned int count,
                                       const std::vector<std::pair<std::string, std::string>> sequences)
{
  const double initial_weight = 1.0f / (double) count;

  std::vector<Particle> particles(count, Particle(initial_weight, sequences));

  return particles;
}

std::vector<Particle> run_smc(std::vector<Particle> &particles,
                              const unsigned int iterations)
{
  for (int i = 0; i < iterations; i++) {
    for (auto &p : particles) {
      resample(p);
      propose(p);
      weight(p);
    }
  }

  return particles;
}

void resample(Particle &particle) {
}

void propose(Particle &particle) {
}

void weight(Particle &particle) {
}
