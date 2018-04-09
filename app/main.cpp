#include <iostream>
#include <algorithm>
#include <float.h>

#include "pll_smc.h"
#include "fasta_helper.h"

void print_tree(phylo_tree_node* root) {
  if (root->left && root->right) {
    std::cout << "(";
    print_tree(root->left);
    std::cout << ", ";
    print_tree(root->right);
    std::cout << "):" << root->length;
  } else {
    std::string label(root->label);
    std::cout << label << ":" << root->length;
  }
}

int main(int argc, char* argv[]) {
  unsigned int particle_count;
  if (argc < 2) {
    std::cerr << "Missing Fasta file path argument!";
    return 1;
  } else if (argc == 2) {
    particle_count = 1000;
  } else if (argc == 3) {
    particle_count = atoi(argv[2]);
  }

  std::vector<std::pair<std::string, std::string>> sequences = parse_sequences(argv[1]);

  std::cerr << "Creating 2 * " << particle_count << " particles" << std::endl;
  std::vector<Particle*> particles = create_particles(2*particle_count, sequences);

  std::cerr << "Running SMC for " << sequences.size() - 1 << " iterations" << std::endl;
  run_smc(particles, sequences.size());

  Particle* particle = nullptr;
  double max = __DBL_MIN__;

  for (auto &p : particles) {
    if (p->get_roots().size() > 1) continue;
    if (p->weight > max) {
      max = p->weight;
      particle = p;
    }
  }
  std::cerr << std::endl;

  std::cerr << "Weight: " << particle->weight << ", Normalized weight: " << particle->normalized_weight << std::endl;

  assert(particle->get_roots().size() == 1);
  print_tree(particle->get_roots()[0]);
  std::cout << ";" << std::endl;
}
