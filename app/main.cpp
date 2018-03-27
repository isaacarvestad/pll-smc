#include <iostream>

#include "pll_smc.h"


void print_tree(pll_rnode_s* root) {
  if (root->left && root->right) {
    std::cout << "(* : ";
    print_tree(root->left);
    std::cout << ", ";
    print_tree(root->right);
    std::cout << ")";
  } else {
    std::string label(root->label);
    std::cout << label;
  }
}


int main() {
  std::cout << "Hello world!" << std::endl;

  std::vector<std::pair<std::string, std::string>> sequences =
    { {"H1", "TAAAAC"},
      {"H2", "CACACG"},
      {"H3", "AGGACA"},
      {"H4", "CGTAGT"},
      {"H5", "CGAATT"}
    };

  std::vector<Particle> particles = create_particles(10000, sequences);
  particles = run_smc(particles, 4);

  Particle* particle = nullptr;
  double max = 0.0f;
  for (auto &p : particles) {
    if (p.weight > max) {
      max = p.weight;
      particle = &p;
    }
  }
  std::cout << std::endl;

  std::cout << "Likelihood: " << particle->weight << std::endl;

  assert(particle->get_roots().size() == 1);
  print_tree(particle->get_roots()[0]);
  std::cout << std::endl;
}
