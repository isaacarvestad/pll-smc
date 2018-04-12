#include <iostream>
#include <algorithm>
#include <float.h>

#include "pll_smc.h"
#include "fasta_helper.h"

void print_tree(phylo_tree_node* root) {
  if (root->child_edge_l && root->child_edge_r) {
    phylo_tree_node* left = root->child_edge_l->child;
    phylo_tree_node* right = root->child_edge_r->child;

    std::cout << "(";
    print_tree(left);
    std::cout << ":" << root->child_edge_l->length;
    std::cout << ", ";
    print_tree(right);
    std::cout << ":" << root->child_edge_r->length;
    std::cout << ")";
  } else {
    std::string label(root->label);
    std::cout << label;
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

  std::cerr << "Running SMC for " << sequences.size() - 1 <<
    " iterations with 2 * " << particle_count << " particles" << std::endl;

  std::vector<Particle*> particles = run_smc(2*particle_count, sequences);

  Particle* particle = nullptr;
  double max = __DBL_MIN__;

  for (auto &p : particles) {
    if (p->get_roots().size() > 1) continue;
    if (p->weight > max) {
      max = p->weight;
      particle = p;
    }
  }

  std::cerr << "Weight: " << particle->weight << ", Normalized weight: " << particle->normalized_weight << std::endl;

  assert(particle->get_roots().size() == 1);
  print_tree(particle->get_roots()[0]);
  std::cout << ";" << std::endl;
}
