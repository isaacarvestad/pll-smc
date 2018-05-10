#include <algorithm>
#include <float.h>
#include <iostream>
#include <memory>

#include "fasta_helper.h"
#include "pll_smc.h"

void print_tree(std::shared_ptr<PhyloTreeNode> root, std::ostream &stream) {
  if (root->edge_l && root->edge_r) {
    stream << "(";
    print_tree(root->edge_l->child, stream);
    stream << ":" << root->edge_l->length;
    stream << ", ";
    print_tree(root->edge_r->child, stream);
    stream << ":" << root->edge_r->length;
    stream << ")";
  } else {
    std::string label(root->label);
    stream << label;
  }
}

int main(int argc, char *argv[]) {
  unsigned int particle_count = 1000;
  if (argc < 2) {
    std::cerr << "Missing Fasta file path argument!";
    return 1;
  } else if (argc == 3) {
    particle_count = atoi(argv[2]);
  }

  std::vector<std::pair<std::string, std::string>> sequences =
      parse_sequences(argv[1]);

  std::cerr << "Running SMC for " << sequences.size() - 1 << " iterations with "
            << particle_count << " particles" << std::endl;

  std::vector<Particle *> particles = run_smc(particle_count, sequences);

  Particle *particle = nullptr;
  double max = -DBL_MAX;

  for (auto &p : particles) {
    if (p->get_roots().size() > 1)
      continue;
    if (p->normalized_weight > max) {
      max = p->normalized_weight;
      particle = p;
    }

    std::cout << p->normalized_weight << " ";
    print_tree(p->get_roots().front(), std::cout);
    std::cout << ";" << std::endl;
  }

  if (particle) {
    std::cerr << "Weight: " << particle->weight
              << ", Normalized weight: " << particle->normalized_weight
              << std::endl;

    assert(particle->get_roots().size() == 1);
    print_tree(particle->get_roots()[0], std::cerr);
    std::cout << ";" << std::endl;
  } else {
    std::cerr << "Couldn't find particle with largest normalized weight"
              << std::endl;
  }
}
