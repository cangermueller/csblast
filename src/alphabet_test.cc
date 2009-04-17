#include <gtest/gtest.h>

#include <iostream>

#include "amino_acid.h"
#include "nucleotide.h"
#include "utils-inl.h"

namespace cs {

TEST(AlphabetTest, AminoAcidAlphabet) {
  for (int i = 1; i < 100; ++i) {
    double x = 0.0 - 0.9 * i;
    std::cout << pow(2.0, x) << std::endl;
    std::cout << fast_pow2(x) << std::endl << std::endl;
  }

  EXPECT_EQ(20, AminoAcid::instance().size());
  EXPECT_EQ('R', AminoAcid::instance().itoc(AminoAcid::instance().ctoi('R')));
  EXPECT_EQ(1, AminoAcid::instance().ctoi(AminoAcid::instance().itoc(1)));
  EXPECT_EQ('A', AminoAcid::instance().itoc(0));
}

TEST(AlphabetTest, NucleotideAlphabet) {
  EXPECT_EQ(4, Nucleotide::instance().size());
  EXPECT_EQ('C', Nucleotide::instance().itoc(Nucleotide::instance().ctoi('C')));
  EXPECT_EQ(1, Nucleotide::instance().ctoi( Nucleotide::instance().itoc(1)));
  EXPECT_EQ('A', Nucleotide::instance().itoc(0));
}

}  // namespace cs
