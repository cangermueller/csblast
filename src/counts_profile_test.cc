#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "amino_acid.h"
#include "matrix_pseudocounts.h"
#include "nucleotide.h"
#include "counts_profile.h"
#include "shared_ptr.h"
#include "nucleotide_matrix.h"

namespace cs
{

const float DELTA = 0.01f;

TEST(CountsProfileTest, ConstructionFromInputStream)
{
    std::string data;
    data.append("CountsProfile\n");
    data.append("num_cols\t\t4\n");
    data.append("alphabet_size\t4\n");
    data.append("logspace\t0\n");
    data.append("has_counts\t0\n");
    data.append(" \tA\tC\tG\tT\tneff\n");
    data.append("1\t0\t*\t*\t*\t1000\n");
    data.append("2\t*\t0\t*\t*\t1000\n");
    data.append("3\t*\t*\t0\t*\t1000\n");
    data.append("4\t*\t*\t*\t0\t1000\n");
    data.append("//\n");
    std::istringstream ss(data);

    CountsProfile<Nucleotide> profile(ss);

    EXPECT_EQ(4, profile.num_cols());
    EXPECT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
    EXPECT_FLOAT_EQ(1.0f, profile.neff(1));
    EXPECT_FALSE(profile.has_counts());
}

TEST(CountsProfileTest, ConstructionFromAlignment)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    CountsProfile<Nucleotide> profile(alignment, true); // use position-dependent weights

    EXPECT_FLOAT_EQ(1.0f, profile.neff(3));
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('A')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('C')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('G')]);
    EXPECT_FLOAT_EQ(1.0f, profile[3][Nucleotide::instance().ctoi('T')]);
}

TEST(CountsProfileTest, ConversionToCounts)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA\n");
    data.append(">seq3\nGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGT\n");
    data.append(">seq4\nCGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACG\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    CountsProfile<Nucleotide> profile(alignment, true);
    ASSERT_EQ(0.25f, profile[3][Nucleotide::instance().ctoi('A')]);

    profile.convert_to_counts();

    EXPECT_TRUE(profile.has_counts());
    EXPECT_FLOAT_EQ(0.25*profile.neff(3), profile[3][Nucleotide::instance().ctoi('A')]);
    EXPECT_FLOAT_EQ(profile.neff(3), std::accumulate(profile.col_begin(3), profile.col_end(3), 0.0f));
}

TEST(CountsProfileTest, AlignmentBpdS)
{
    std::ifstream fin("../data/BpdS.fas");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::FASTA);
    fin.close();
    alignment.assign_match_columns_by_gap_rule();
    CountsProfile<AminoAcid> profile(alignment, true);

    EXPECT_NEAR(0.0f, profile[10][AminoAcid::instance().ctoi('H')], DELTA);
    EXPECT_FLOAT_EQ(1.0f, std::accumulate(profile.col_begin(122), profile.col_end(122), 0.0f));
}

TEST(CountsProfileTest, Alignment1Q7L)
{
    std::ifstream fin("../data/1Q7L.fas");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::FASTA);
    fin.close();
    alignment.assign_match_columns_by_gap_rule();
    CountsProfile<AminoAcid> profile(alignment, true);

    EXPECT_NEAR(0.05f, profile[10][AminoAcid::instance().ctoi('G')], DELTA);
}

TEST(CountsProfileTest, AlignmentCelegansRefGene)
{
    std::ifstream fin("../data/ce_refgene.fas");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fin.close();
    CountsProfile<Nucleotide> profile(alignment, false);

    EXPECT_FLOAT_EQ(1.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[1][Nucleotide::instance().ctoi('T')]);
}

TEST(CountsProfileTest, AddMatrixPseudocountsToProfile)
{
    std::ifstream fin("../data/ce_refgene.fas");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fin.close();
    CountsProfile<Nucleotide> profile(alignment, false);

    ASSERT_FLOAT_EQ(1.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    NucleotideMatrix m(1,-1);
    DivergenceDependentAdmixture pca(1.0f, 10.0f);
    MatrixPseudocounts<Nucleotide> mpc(&m, &pca);
    mpc.add_to_profile(profile);

    EXPECT_NEAR(0.25f, profile[0][Nucleotide::instance().ctoi('T')], DELTA);
}

TEST(CountsProfileTest, AddMatrixPseudocountsToLogProfile)
{
    std::ifstream fin("../data/ce_refgene.fas");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fin.close();
    CountsProfile<Nucleotide> profile(alignment, false);
    profile.transform_to_logspace();

    ASSERT_FLOAT_EQ(0.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    NucleotideMatrix m(1,-1);
    DivergenceDependentAdmixture pca(1.0f, 10.0f);
    MatrixPseudocounts<Nucleotide> mpc(&m, &pca);
    mpc.add_to_profile(profile);

    profile.transform_to_linspace();
    EXPECT_NEAR(0.25f, profile[0][Nucleotide::instance().ctoi('T')], DELTA);
}

};  // cs
