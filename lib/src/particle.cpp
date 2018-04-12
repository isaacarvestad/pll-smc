#include "particle.h"

Particle::Particle(double weight,
                   const std::vector<std::pair<std::string, std::string>> sequences,
                   const unsigned int sequence_lengths,
                   const pll_partition_t* reference_partition)
  : weight(weight),
    normalized_weight(1)
{
  forest = new PhyloForest(sequences, sequence_lengths, reference_partition);

  std::random_device random;
  mt_generator = std::mt19937(random());
}

Particle::Particle(const Particle &original)
  : weight(original.weight),
    normalized_weight(original.normalized_weight)
{
  forest = new PhyloForest(*original.forest);

  std::random_device random;
  mt_generator = std::mt19937(random());
}

Particle& Particle::operator=(const Particle& original) {
  if (this == &original) return *this;

  weight = original.weight;
  *forest = *original.forest;

  return *this;
}

Particle::~Particle() {
  delete(forest);
}

void Particle::propose(const double rate) {
  assert(forest->root_count() > 1 && "Cannot propose a continuation on a single root node");

  std::uniform_int_distribution<int> int_dist(0, forest->root_count() - 1);
  int i = int_dist(mt_generator);
  int j = -1;

  while (j == -1) {
    j = int_dist(mt_generator);
    if (j == i) j = -1;
  }

  std::exponential_distribution<double> exponential_dist(rate);
  double height = exponential_dist(mt_generator);

  phylo_tree_node* node = forest->connect(i, j, height);

  weight = forest->likelihood_factor(node);
  assert(!isnan(weight) && !isinf(weight));
}
