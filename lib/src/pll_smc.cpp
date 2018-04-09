#include "pll_smc.h"

/**
   Approximates the exponential branch length. See
   https://en.wikipedia.org/wiki/Binomial_coefficient#Bounds_and_asymptotic_formulas
 */
double approx_rate(int x) {
  double n = (double) x;
  double k = 2;

  return (pow(n/k - 0.5, k) * exp(k)) / (sqrt(2 * M_PI * k));
}

std::vector<Particle*> create_particles(const unsigned int count,
                                        const std::vector<std::pair<std::string, std::string>> sequences)
{
  assert(sequences.size() > 0 && "Expected at least one sequence");
  const unsigned int sequence_lengths = sequences[0].second.length();
  for (auto &s : sequences) {
    assert(s.second.length() == sequence_lengths && "Sequence lengths do not match");
  }

  const double initial_weight = 1.0 / (double) count;

  std::vector<Particle*> particles(count, nullptr);
  for (auto &p : particles) {
    p = new Particle(initial_weight, sequences, sequence_lengths);
  }

  return particles;
}

std::vector<Particle*> run_smc(std::vector<Particle*> &particles,
                               const unsigned int sequence_count)
{
  const int iterations = sequence_count - 1;

  for (int i = 0; i < iterations; i++) {
    std::cerr << "Iteration " << i << std::endl;
    double rate = approx_rate(i - sequence_count + 1);

    resample(particles, i);
    propose(particles, rate, i);
    normalize_weights(particles, i);
  }

  std::random_device random;
  std::mt19937 generator(random());
  for (auto &p : particles) {
    std::exponential_distribution<double> exponential_dist(1);
    p->get_roots()[0]->length = exponential_dist(generator);
  }

  return particles;
}

void resample(std::vector<Particle*> &particles, const unsigned int iteration) {
  int offset = iteration % 2 == 0 ? 0 : particles.size() / 2;

  std::vector<double> normalized_weights;
  for (int i = offset; i < particles.size() / 2 + offset; i++) {
     normalized_weights.push_back(particles[i]->normalized_weight);
  }

  std::random_device random;
  std::discrete_distribution<double> dist(normalized_weights.begin(), normalized_weights.end());

  for (int i = particles.size() / 2 - offset; i < particles.size() - offset; i++) {
    int index = dist(random);

    particles[i]->shallow_copy(*particles[index + offset]);
  }
}

void propose(std::vector<Particle*> &particles, const double rate, const unsigned int iteration) {
  int offset = iteration % 2 == 0 ? particles.size() / 2 : 0;

  for (int i = offset; i < particles.size() / 2 + offset; i++) {
    particles[i]->propose(rate);
  }
}

void normalize_weights(std::vector<Particle*> &particles, const unsigned int iteration) {
  int offset = iteration % 2 == 0 ? particles.size() / 2 : 0;

  double max = __DBL_MIN__;
  for (int i = offset; i < particles.size() / 2 + offset; i++) {
    if (particles[i]->weight > max) max = particles[i]->weight;
  }

  double sum = 0.0;
  for (int i = offset; i < particles.size() / 2 + offset; i++) {
    particles[i]->normalized_weight = particles[i]->weight - max;
    sum += exp(particles[i]->normalized_weight);
  }

  for (int i = offset; i < particles.size() / 2 + offset; i++) {
    particles[i]->normalized_weight = exp(particles[i]->normalized_weight) / sum;
  }
}
