#include "particle.h"

Particle::Particle(double weight, const std::vector<std::string> sequences)
  : weight(weight),
    forest(sequences, 0)
{
}
