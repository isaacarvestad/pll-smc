#include "pll_smc.h"

double compute_rate(int leaves, int iteration) {
  int n = leaves - iteration + 1;
  int np = n - 1;

  return ((double) (n * np)) / 2.0;
}

/**
   Creates a vector with 'count' number of particles, each using the given
   vector of sequences.

   Each particle starts with a weight of 1/'count'.
 */
std::vector<Particle*> create_particles(const unsigned int count,
                                        const std::vector<std::pair<std::string, std::string>> sequences,
                                        const pll_partition_t* reference_partition,
                                        PLLBufferManager* const pll_buffer_manager)
{
  assert(sequences.size() > 0 && "Expected at least one sequence");
  const unsigned int sequence_lengths = sequences[0].second.length();
  for (auto &s : sequences) {
    assert(s.second.length() == sequence_lengths && "Sequence lengths do not match");
  }

  const double initial_weight = log(1.0 / (double) count);

  Particle particle = Particle(initial_weight, sequences, sequence_lengths, reference_partition, pll_buffer_manager);
  std::vector<Particle*> particles(count * 2, nullptr);
  for (auto &p : particles) {
    p = new Particle(particle);
  }

  return particles;
}

/**

 */
const pll_partition_t* create_reference_partition(const std::vector<std::pair<std::string, std::string>> sequences) {
  assert(sequences.size() > 0 && "Expected at least one sequence");
  const unsigned int sequence_lengths = sequences[0].second.length();
  for (auto &s : sequences) {
    assert(s.second.length() == sequence_lengths && "Sequence lengths do not match");
  }

  const unsigned int rate_category_count = 4;
  const double rate_categories[4] =
  { 0.13695378267140107,
    0.47675185617665189,
    0.99999999997958422,
    2.38629436117236260
  };

  const unsigned int subst_model_count = 1;
  const double subst_params[6] = { 1, 1, 1, 1, 1, 1 };

  const unsigned int nucleotide_states = 4;
  const double nucleotide_frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };

  pll_partition* partition = pll_partition_create(sequences.size(),
                                                  0, // Don't allocate any inner CLV's.
                                                  nucleotide_states,
                                                  sequence_lengths,
                                                  subst_model_count,
                                                  0, // Don't allocate any pmatrices.
                                                  rate_category_count,
                                                  0, // Don't allocate any scale buffers.
                                                  PLL_ATTRIB_ARCH_SSE);

  assert(partition);
  pll_set_frequencies(partition, 0, nucleotide_frequencies);
  pll_set_category_rates(partition, rate_categories);
  pll_set_subst_params(partition, 0, subst_params);

  for (unsigned int i = 0; i < sequences.size(); i++) {
    std::string sequence = sequences[i].second;

    pll_set_tip_states(partition, i, pll_map_nt, sequence.data());
  }

  // Once for each param index
  pll_update_eigen(partition, 0);

  return partition;
}

std::vector<Particle*> run_smc(const unsigned int particle_count,
                               const std::vector<std::pair<std::string, std::string>> sequences) {

  const pll_partition_t* reference_partition = create_reference_partition(sequences);
  PLLBufferManager* pll_buffer_manager = new PLLBufferManager;

  std::vector<Particle*> particles = create_particles(particle_count, sequences, reference_partition, pll_buffer_manager);

  const unsigned int sequence_count = sequences.size();
  const unsigned int iterations = sequence_count - 1;

  for (int i = 0; i < iterations; i++) {
    std::cerr << "Iteration " << i << std::endl;
    double rate = compute_rate(sequence_count, i);

    resample(particles, i);
    propose(particles, rate, i);
    normalize_weights(particles, i);
  }

  return particles;
}

void resample(std::vector<Particle*> &particles, const unsigned int iteration) {
  int offset = iteration % 2 == 0 ? 0 : particles.size() / 2;

  std::vector<double> normalized_weights;
  double ess_sum = 0.0;
  for (int i = offset; i < particles.size() / 2 + offset; i++) {
    double normalized_weight = particles[i]->normalized_weight;

    normalized_weights.push_back(normalized_weight);
    ess_sum += normalized_weight * normalized_weight;

    assert(normalized_weight <= 1.0);
    assert(ess_sum <= 1);
  }
  double ess = 1 / ess_sum;

  std::cerr << "ESS: " << ess << std::endl;

  std::random_device random;
  std::discrete_distribution<int> dist(normalized_weights.begin(), normalized_weights.end());

  for (int i = particles.size() / 2 - offset; i < particles.size() - offset; i++) {
    int index = dist(random);

    *particles[i] = *particles[index + offset];
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
