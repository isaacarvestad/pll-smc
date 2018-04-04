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
  std::cerr << "Creating particles" << std::endl;

  assert(sequences.size() > 0 && "Expected at least one sequence");
  const unsigned int sequence_lengths = sequences[0].second.length();
  for (auto &s : sequences) {
    assert(s.second.length() == sequence_lengths && "Sequence lengths do not match");
  }

  const double initial_weight = 1.0f / (double) count;

  std::vector<Particle*> particles(count, nullptr);
  for (auto &p : particles) {
    p = new Particle(initial_weight, sequences, sequence_lengths);
  }

  std::cerr << "Particles successfully created" << std::endl;

  return particles;
}

std::vector<Particle*> run_smc(std::vector<Particle*> &particles,
                               const unsigned int sequence_count)
{
  const int iterations = sequence_count - 1;

  std::cerr << "Running SMC with " << iterations << " iterations" << std::endl;

  for (int i = 0; i < iterations; i++) {
    std::cerr << "Iteration " << i << std::endl;
    double rate = approx_rate(i - sequence_count + 1);

    resample(particles);
    propose(particles, rate);
    normalize_weights(particles);
  }

  std::random_device random;
  std::mt19937 generator(random());
  for (auto &p : particles) {
    std::exponential_distribution<> exponential_dist(1);
    p->get_roots()[0]->length = exponential_dist(generator);
  }

  std::cerr << "Finished running SMC" << std::endl;

  return particles;
}

void resample(std::vector<Particle*> &particles) {
  std::vector<double> normalized_weights;
  for (auto &p : particles) {
    normalized_weights.push_back(p->normalized_weight);
  }

  std::random_device random;
  std::discrete_distribution<double> dist(normalized_weights.begin(), normalized_weights.end());

  std::vector<Particle*> new_particles;
  for (int i = 0; i < particles.size(); i++) {
    int index = dist(random);

    Particle* p = new Particle(*particles[index]);
    new_particles.push_back(p);
  }

  for (int i = 0; i < particles.size(); i++) {
    delete(particles[i]);
  }

  particles = new_particles;
}

void propose(std::vector<Particle*> &particles, const double rate) {
  for (auto &p : particles) {
    p->propose(rate);
  }
}

void normalize_weights(std::vector<Particle*> &particles) {
  double max = __DBL_MIN__;
  for (auto &p : particles) {
    if (p->weight > max) max = p->weight;
  }

  double sum = 0.0;
  for (auto &p : particles) {
    p->normalized_weight = p->weight - max;
    sum += exp(p->normalized_weight);
  }

  for (auto &p : particles) {
    p->normalized_weight = exp(p->normalized_weight) / sum;
  }
}
