#include <iostream>
#include <float.h>

#include "pll_smc.h"


void print_tree(pll_rnode_s* root) {
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


int main() {
  std::cout << "Hello world!" << std::endl;

  std::vector<std::pair<std::string, std::string>> sequences =
    { {"H0", "TATTGCTAATGCTACATCGCTCGGAACAAATACGGGGCACTTGTTTTCAGCACAAGCAAGCGGGAAGTCATGACCCCCTAGAGTCCCTTCTGGAAAGCCT"},
      {"H1", "GAGGCTACCCGCTGACAGTTCCCGAACAGGAATAATTATTGGGTGTTCTGATCGAACAGGGTATGGCAATGGCCATCACGAGCCGAGCGGCCTCCCAAAG"},
      {"H2", "GATACATCGGCGTTTACTACCCCTAGCTGAGTTCCTTTCAGTTGCGTTTGGGCGATTCAGGTTTACCCGGAGACGAACTCACTCAACTTTTAAACTAAGA"},
      {"H4", "TGAGACTGTTTGACGGTCTAGCGGATTCTCGTCAATGTAGTACCATATAGCGGGGCGTTCTCAACCATTTAGCGGAAATTAATGCGTTGACCCAGTCCCG"},
      {"H5", "CGACTGCGGGACCGCGGCTGCAGGATAGTTAATTAATTAGCGGCTGTTATGTCGCGCCACTCATTAGGGGTACGGCTTTAATTGAAGGTAGGTATTCACG"},
      {"H7", "CGTATCGAGTAACGCCCCCGATGGCTGGTGAGTAAATTAGGCTTGAATAGCGAGACGCACTACTCCGACTACAGATAAAAATTCAAGTTCCTGAGTCTCG"},
      {"H8", "TCTTACGGAAGTCTTATGATATACGCTGTCAGCCCGCCACAGTGTTTAAATAATCATGCTGCGGAAGTGTTCGTATCATGGGAGTTCTGCAGACGCCATA"},
      {"H12", "ACGCAACAGTCGTCTCTAGTTTTACGCAGGAGACTAGTGGCATGGGTTGGCTGTCACATCAGACTCCGGTGTCGTTCTGGGGAAATTGATCGAGGTTAGG"},
      {"H15", "ATTGTCCACACACTCGTTAACAGGCAAAAAGTCCGTCTAAAAGGGCTTAGGTGGAGCTCGATCCTGCAGTTCAACATCGAGTAACCACACTCGTACTGCG"},
      {"H16", "ATCCCTCTTACAAGCCTTTCCAGTCCATGGGGCTGAGTAAAATAGACAGTGCCTCGCGCGAGCACGGAGAATCTAATCATGAACCCCTCCGTCCTTCCGC"},
    };

  std::vector<Particle> particles = create_particles(10000, sequences);
  particles = run_smc(particles, sequences.size());

  Particle* particle = nullptr;
  double max = -DBL_MAX;
  for (auto &p : particles) {
    if (p.weight > max) {
      max = p.weight;
      particle = &p;
    }
  }
  std::cout << std::endl;

  std::cout << "Weight: " << particle->weight << ", normalized weight: " << particle->normalized_weight << std::endl;

  assert(particle->get_roots().size() == 1);
  print_tree(particle->get_roots()[0]);
  std::cout << ";" << std::endl;
}
