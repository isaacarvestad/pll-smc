#include <iostream>
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
  if (argc != 2) {
    std::cerr << "Missing Fasta file path argument!";
    return 1;
  }

  std::vector<std::pair<std::string, std::string>> sequences = parse_sequences(argv[1]);

  std::vector<Particle*> particles = create_particles(500, sequences);
  run_smc(particles, sequences.size());

  Particle* particle = nullptr;
  double max = __DBL_MIN__;
  for (auto &p : particles) {
    if (p->normalized_weight > max) {
      max = p->normalized_weight;
      particle = p;
    }
  }
  std::cout << std::endl;

  std::cout << "Weight: " << particle->weight << ", Normalized weight: " << particle->normalized_weight << std::endl;

  assert(particle->get_roots().size() == 1);
  print_tree(particle->get_roots()[0]);
  std::cout << ";" << std::endl;
}
