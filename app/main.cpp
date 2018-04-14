#include <iostream>
#include <algorithm>
#include <float.h>

#include "pll_smc.h"
#include "fasta_helper.h"

void print_tree(phylo_tree_node* root, std::ostream &stream) {
  if (root->left && root->right) {
    stream << "(";
    print_tree(root->left, stream);
    stream << ", ";
    print_tree(root->right, stream);
    stream << "):" << root->length;
  } else {
    std::string label(root->label);
    stream << label << ":" << root->length;
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
  double max = -DBL_MAX;

  for (auto &p : particles) {
    if (p->get_roots().size() > 1) continue;
    if (p->normalized_weight > max) {
      max = p->normalized_weight;
      particle = p;
    }

    std::cout << p->normalized_weight << " ";
    print_tree(p->get_roots().front(), std::cout);
    std::cout << ";" << std::endl;
  }
  std::cerr << std::endl;

  if (particle) {
    std::cerr << "Weight: " << particle->weight << ", Normalized weight: " << particle->normalized_weight << std::endl;

    assert(particle->get_roots().size() == 1);
    print_tree(particle->get_roots()[0], std::cerr);
    std::cout << ";" << std::endl;
  } else {
    std::cerr << "Couldn't find particle with largest normalized weight" << std::endl;
  }
}
