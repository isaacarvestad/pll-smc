#ifndef LIB_PLL_SMC_FASTA_HELPER_H
#define LIB_PLL_SMC_FASTA_HELPER_H

#include <string>
#include <vector>

#include <libpll/pll.h>

/**
   Parses a Fasta file and returns a vector containing pairs of (name,
   sequence).

   Implementation based on PLL "newick-fasta-unrooted" example.
 */
std::vector<std::pair<std::string, std::string>> parse_sequences(std::string file_path);

#endif
