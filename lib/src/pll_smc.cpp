#include "pll_smc.h"

std::vector<Particle> create_particles(const unsigned int count,
                                       const std::vector<std::pair<std::string, std::string>> sequences)
{
  assert(sequences.size() > 0 && "Expected at least one sequence");
  const unsigned int sequence_lengths = sequences[0].second.length();
  for (auto &s : sequences) {
    assert(s.second.length() == sequence_lengths && "Sequence lengths do not match");
  }

  const double initial_weight = 1.0f / (double) count;

  std::vector<Particle> particles(count, Particle(initial_weight, sequences, sequence_lengths));

  return particles;
}

std::vector<Particle> run_smc(std::vector<Particle> &particles,
                              const unsigned int iterations)
{
  for (int i = 0; i < iterations; i++) {
    resample(particles);
    propose(particles);
    normalize_weights(particles);
  }

  return particles;
}

void resample(std::vector<Particle> &particles) {
}

void propose(std::vector<Particle> &particles) {
  for (auto &p : particles) {
    p.propose();
  }
}

void normalize_weights(std::vector<Particle> &particles) {
  double sum = 0.0f;
  for (auto &p : particles) {
    sum += p.weight;
  }

  for (auto &p : particles) {
    p.weight /= sum;
  }
}
