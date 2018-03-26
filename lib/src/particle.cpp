#include "particle.h"

Particle::Particle(double weight, const std::vector<std::pair<std::string, std::string>> sequences)
  : weight(weight),
    forest(sequences, 0)
{
}
