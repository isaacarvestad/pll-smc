#include "particle.h"

Particle::Particle(double weight,
                   const std::vector<std::pair<std::string, std::string>> sequences,
                   const unsigned int sequence_lengths)
  : weight(weight),
    normalized_weight(1)
{
  forest = new PhyloForest(sequences, sequence_lengths);

  std::random_device random;
  mt_generator = std::mt19937(random());
}

Particle::Particle(const Particle &original)
  : weight(original.weight),
    normalized_weight(original.normalized_weight)
{
  forest = new PhyloForest(*original.get_forest());

  std::random_device random;
  mt_generator = std::mt19937(random());
}

Particle::~Particle() {
  delete(forest);
}

void Particle::shallow_copy(const Particle &original) {
  weight = original.weight;
  normalized_weight = original.normalized_weight;

  forest->shallow_copy(*original.get_forest());
}

void Particle::propose() {
  assert(forest->root_count() > 1 &&
         "Cannot propose a continuation on a single root node");

  std::uniform_int_distribution<int> int_dist(0, forest->root_count() - 1);
  int i = int_dist(mt_generator);
  int j = -1;

  while (j == -1) {
    j = int_dist(mt_generator);
    if (j == i) j = -1;
  }

  double root_count = forest->root_count();
  double rate = root_count * (root_count - 1) / 2;

  std::exponential_distribution<double> exponential_dist(rate);
  double height = exponential_dist(mt_generator);
  double ln_height_prior = log(rate) * (-rate * height);

  phylo_tree_node* node = forest->connect(i, j, height);

  weight = forest->likelihood_factor(node) + ln_height_prior;

  assert(!isnan(weight) && !isinf(weight));
}
