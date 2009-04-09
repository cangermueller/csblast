#include <gtest/gtest.h>

#include "amino_acid.h"
#include "nucleotide.h"

namespace cs
{

TEST(AlphabetTest, AminoAcidAlphabet)
{
    EXPECT_EQ(20, AminoAcid::instance().size());
    EXPECT_EQ('R', AminoAcid::instance().itoc(AminoAcid::instance().ctoi('R')));
    EXPECT_EQ(1, AminoAcid::instance().ctoi(AminoAcid::instance().itoc(1)));
    EXPECT_EQ('A', AminoAcid::instance().itoc(0));
}

TEST(AlphabetTest, NucleotideAlphabet)
{
    EXPECT_EQ(4, Nucleotide::instance().size());
    EXPECT_EQ('C', Nucleotide::instance().itoc(Nucleotide::instance().ctoi('C')));
    EXPECT_EQ(1, Nucleotide::instance().ctoi( Nucleotide::instance().itoc(1)));
    EXPECT_EQ('A', Nucleotide::instance().itoc(0));
}

}  // namespace cs
