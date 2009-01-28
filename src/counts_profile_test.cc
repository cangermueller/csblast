#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "amino_acid_alphabet.h"
#include "matrix_pseudocounts.h"
#include "nucleotide_alphabet.h"
#include "counts_profile.h"
#include "shared_ptr.h"
#include "nucleotide_matrix.h"

const float kDelta = 0.01f;

TEST(CountsProfileTest, ConstructionFromInputStream)
{
    std::string data;
    data.append("CountsProfile\n");
    data.append("ncols\t\t4\n");
    data.append("nalph\t\t4\n");
    data.append("logspace\t0\n");
    data.append("has_counts\t0\n");
    data.append(" \tA\tC\tG\tT\tneff\n");
    data.append("1\t0\t*\t*\t*\t1000\n");
    data.append("2\t*\t0\t*\t*\t1000\n");
    data.append("3\t*\t*\t0\t*\t1000\n");
    data.append("4\t*\t*\t*\t0\t1000\n");
    data.append("//\n");
    std::istringstream ss(data);

    cs::CountsProfile profile(ss, cs::NucleotideAlphabet::instance());

    EXPECT_EQ(4, profile.ncols());
    EXPECT_EQ(4, profile.nalph());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
    EXPECT_FLOAT_EQ(1.0f, profile.neff(1));
    EXPECT_FALSE(profile.has_counts());
}

TEST(CountsProfileTest, ConstructionFromAlignment)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA, na);

    cs::CountsProfile profile(alignment, true); // use position-dependent weights

    EXPECT_FLOAT_EQ(1.0f, profile.neff(3));
    EXPECT_FLOAT_EQ(0.0f, profile[3][na->ctoi('A')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][na->ctoi('C')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][na->ctoi('G')]);
    EXPECT_FLOAT_EQ(1.0f, profile[3][na->ctoi('T')]);
}

TEST(CountsProfileTest, ConversionToCounts)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA\n");
    data.append(">seq3\nGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGT\n");
    data.append(">seq4\nCGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACG\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA, na);

    cs::CountsProfile profile(alignment, true);
    ASSERT_EQ(0.25f, profile[3][na->ctoi('A')]);

    profile.convert_to_counts();

    EXPECT_TRUE(profile.has_counts());
    EXPECT_FLOAT_EQ(0.25*profile.neff(3), profile[3][na->ctoi('A')]);
    EXPECT_FLOAT_EQ(profile.neff(3), std::accumulate(profile.col_begin(3), profile.col_end(3), 0.0f));
}

TEST(CountsProfileTest, AlignmentBpdS)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    std::ifstream fin("../data/BpdS.fas");
    cs::Alignment alignment(fin, cs::Alignment::FASTA, aa);
    fin.close();
    alignment.assign_match_columns_by_gap_rule();
    cs::CountsProfile profile(alignment, true);

    EXPECT_NEAR(0.0f, profile[10][aa->ctoi('H')], kDelta);
    EXPECT_FLOAT_EQ(1.0f, std::accumulate(profile.col_begin(122), profile.col_end(122), 0.0f));
}

TEST(CountsProfileTest, Alignment1Q7L)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    std::ifstream fin("../data/1Q7L.fas");
    cs::Alignment alignment(fin, cs::Alignment::FASTA, aa);
    fin.close();
    alignment.assign_match_columns_by_gap_rule();
    cs::CountsProfile profile(alignment, true);

    EXPECT_NEAR(0.05f, profile[10][aa->ctoi('G')], kDelta);
}

TEST(CountsProfileTest, AlignmentCelegansRefGene)
{
    cs::NucleotideAlphabet* aa = cs::NucleotideAlphabet::instance();
    std::ifstream fin("../data/ce_refgene.fas");
    cs::Alignment alignment(fin, cs::Alignment::FASTA, aa);
    fin.close();
    cs::CountsProfile profile(alignment, false);

    EXPECT_FLOAT_EQ(1.0f, profile[1][aa->ctoi('T')]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[1][aa->ctoi('T')]);
}

TEST(CountsProfileTest, MatrixPseudocounts)
{
    cs::NucleotideAlphabet* aa = cs::NucleotideAlphabet::instance();
    std::ifstream fin("../data/ce_refgene.fas");
    cs::Alignment alignment(fin, cs::Alignment::FASTA, aa);
    fin.close();
    cs::CountsProfile profile(alignment, false);

    ASSERT_FLOAT_EQ(1.0f, profile[1][aa->ctoi('T')]);

    cs::MatrixPseudocounts mpc(cs::NucleotideMatrix(1, -1));
    mpc.add_to_profile(profile, cs::ProfileSequenceAdmixture(1.0f, 10.0f));

    EXPECT_NEAR(0.25f, profile[0][aa->ctoi('T')], kDelta);
}
