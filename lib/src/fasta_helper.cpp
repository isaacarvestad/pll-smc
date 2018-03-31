#include "fasta_helper.h"

std::vector<std::pair<std::string, std::string>> parse_sequences(std::string file_path) {
  pll_fasta_t* fasta = pll_fasta_open(file_path.c_str(), pll_map_fasta);

  std::vector<std::pair<std::string, std::string>> sequences;

  char* header;
  char* sequence;
  long sequence_length;
  long header_length;
  long sequence_number;

  while(pll_fasta_getnext(fasta, &header, &header_length, &sequence, &sequence_length, &sequence_number)) {
    std::string h(header);
    std::string s(sequence);

    sequences.push_back({h, s});
  }

  return sequences;
}
