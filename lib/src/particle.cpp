#include "particle.h"

Particle::Particle(double weight,
                   const std::vector<std::pair<std::string, std::string>> sequences,
                   const unsigned int sequence_lengths)
  : weight(weight),
    forest(sequences, sequence_lengths)
{
}

void Particle::propose() {
  std::random_device random;
  std::uniform_int_distribution<int> dist(0, forest.root_count() - 1);
  int i = dist(random);
  int j = -1;

  while (j == -1) {
    j = dist(random);
    if (j == i) j = -1;
  }

  pll_rnode_s* node = forest.connect(i, j, 1.0f, 1.0f);

  weight *= forest.likelihood_factor(node);
}
