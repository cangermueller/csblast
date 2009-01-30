#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid_alphabet.h"
#include "blosum_matrix.h"
#include "log.h"
#include "matrix_pseudocounts.h"
#include "nucleotide_alphabet.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"

const float kDelta = 0.01f;

TEST(SequenceTest, ConstructionFromAlphabetVector)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    const cs::Sequence sequence("header", "A RNDCQEGHILKMFPSTWYV\n", aa);

    EXPECT_EQ(aa->size(), sequence.length());
    EXPECT_EQ(aa->ctoi('R'), sequence[1]);
}

TEST(SequenceTest, ConstructionFromInputStream)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::istringstream ss(">dummy header\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC");
    cs::Sequence sequence(ss, na);

    EXPECT_EQ(80, sequence.length());
    EXPECT_EQ(na->ctoi('C'), sequence[1]);
    EXPECT_EQ(na->ctoi('C'), sequence[79]);
}

TEST(SequenceTest, ConstructionOfMultipleSequencesFromInputStream)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    std::vector< shared_ptr<cs::Sequence> > seqs(cs::Sequence::readall(ss, na));

    EXPECT_EQ(2, static_cast<int>(seqs.size()));
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])[1]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])[79]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])[1]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])[79]);
}

TEST(SequenceTest, AddMatrixPseudocountsToSequence)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();

    const cs::Sequence sequence("header", "ARNDCQEGHILKMFPSTWYV", aa);
    cs::Profile profile(sequence.length(), aa);

    ASSERT_EQ(aa->size(), sequence.length());
    ASSERT_EQ(aa->ctoi('R'), sequence[1]);
    ASSERT_EQ(sequence.length(), profile.ncols());

    cs::BlosumMatrix m;
    cs::MatrixPseudocounts mpc(m);
    mpc.add_to_sequence(sequence, cs::DivergenceDependentAdmixture(1.0f, 10.0f), profile);

    EXPECT_NEAR(0.06f, profile[0][aa->ctoi('V')], kDelta);
}
