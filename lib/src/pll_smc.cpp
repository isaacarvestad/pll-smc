#include "pll_smc.h"
#include <iostream>

/**
   Approximates the exponential branch length. See
   https://en.wikipedia.org/wiki/Binomial_coefficient#Bounds_and_asymptotic_formulas
 */
double approx_rate(int x) {
  double n = (double) x;
  double k = 2;

  return (pow(n/k - 0.5, k) * exp(k)) / (sqrt(2 * M_PI * k));
}

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
                              const unsigned int sequence_count)
{
  const int iterations = sequence_count - 1;

  for (int i = 0; i < iterations; i++) {
    double rate = approx_rate(i - sequence_count + 1);

    std::cerr << "Iteration: " << i << " using rate: " << rate << std::endl;

    resample(particles);
    propose(particles, rate);
    normalize_weights(particles);
  }

  std::random_device random;
  std::mt19937 generator(random());
  for (auto &p : particles) {
    std::exponential_distribution<> exponential_dist(1);
    p.get_roots()[0]->length = exponential_dist(generator);
  }

  return particles;
}

void resample(std::vector<Particle> &particles) {
}

void propose(std::vector<Particle> &particles, const double rate) {
  for (auto &p : particles) {
    p.propose(rate);
  }
}

void normalize_weights(std::vector<Particle> &particles) {
  double sum = 0.0;
  for (auto &p : particles) {
    sum += p.weight;
  }

  for (auto &p : particles) {
    p.normalized_weight = p.weight / sum;
  }
}
