#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid.h"
// #include "blosum_matrix.h"
// #include "log.h"
// #include "matrix_pseudocounts.h"
#include "nucleotide.h"
// #include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"

#include "log.cc" // hack!!!

namespace cs
{

const float DELTA = 0.01f;

TEST(SequenceTest, ConstructionFromAlphabetVector)
{
    const Sequence<AminoAcid> sequence("header", "A RNDCQEGHILKMFPSTWYV\n");

    EXPECT_EQ(AminoAcid::instance().size(), sequence.length());
    EXPECT_EQ(AminoAcid::instance().ctoi('R'), sequence[1]);
}

TEST(SequenceTest, ConstructionFromInputStream)
{
    std::istringstream ss(">dummy header\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC");
    const Sequence<Nucleotide> sequence(ss);

    EXPECT_EQ(80, sequence.length());
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), sequence[1]);
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), sequence[79]);
}

TEST(SequenceTest, ConstructionOfMultipleSequencesFromInputStream)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    std::vector< shared_ptr<Sequence<Nucleotide> > > seqs(Sequence<Nucleotide>::readall(ss));

    EXPECT_EQ(2, static_cast<int>(seqs.size()));
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[0])[1]);
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[0])[79]);
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[1])[1]);
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[1])[79]);
}

// TEST(SequenceTest, AddMatrixPseudocountsToSequence)
// {
//     AminoAcidAlphabet* aa = AminoAcidAlphabet::instance();

//     const Sequence sequence("header", "ARNDCQEGHILKMFPSTWYV", aa);
//     Profile profile(sequence.length(), aa);

//     ASSERT_EQ(AminoAcidAlphabet::instance().size(), sequence.length());
//     ASSERT_EQ(AminoAcidAlphabet::instance().ctoi('R'), sequence[1]);
//     ASSERT_EQ(sequence.length(), profile.ncols());

//     BlosumMatrix m;
//     MatrixPseudocounts mpc(m);
//     mpc.add_to_sequence(sequence, DivergenceDependentAdmixture(1.0f, 10.0f), profile);

//     EXPECT_NEAR(0.06f, profile[0][AminoAcidAlphabet::instance().ctoi('V')], DELTA);
// }

}  // cs

