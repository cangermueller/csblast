#include <gtest/gtest.h>

#include <cstdio>

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
    FILE* fin = fopen("../data/nt_counts_profile.prf", "r");
    CountsProfile<Nucleotide> profile(fin);
    fclose(fin);

    EXPECT_EQ(4, profile.num_cols());
    EXPECT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
    EXPECT_FLOAT_EQ(1.0f, profile.neff(1));
    EXPECT_FALSE(profile.has_counts());
}

TEST(CountsProfileTest, ConstructionFromAlignment)
{
    FILE* fin = fopen("../data/nt_alignment7.fas", "r");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fclose(fin);

    CountsProfile<Nucleotide> profile(alignment, true);

    EXPECT_FLOAT_EQ(1.0f, profile.neff(3));
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('A')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('C')]);
    EXPECT_FLOAT_EQ(0.0f, profile[3][Nucleotide::instance().ctoi('G')]);
    EXPECT_FLOAT_EQ(1.0f, profile[3][Nucleotide::instance().ctoi('T')]);
}

TEST(CountsProfileTest, ConversionToCounts)
{
    FILE* fin = fopen("../data/nt_alignment8.fas", "r");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fclose(fin);

    CountsProfile<Nucleotide> profile(alignment, true);
    ASSERT_EQ(0.25f, profile[3][Nucleotide::instance().ctoi('A')]);

    profile.convert_to_counts();

    EXPECT_TRUE(profile.has_counts());
    EXPECT_FLOAT_EQ(0.25*profile.neff(3), profile[3][Nucleotide::instance().ctoi('A')]);
    EXPECT_FLOAT_EQ(profile.neff(3), std::accumulate(profile.col_begin(3),
                                                     profile.col_end(3), 0.0f));
}

TEST(CountsProfileTest, AlignmentBpdS)
{
    FILE* fin = fopen("../data/BpdS.fas", "r");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::FASTA);
    fclose(fin);

    alignment.assign_match_columns_by_gap_rule();
    CountsProfile<AminoAcid> profile(alignment, true);

    EXPECT_NEAR(0.0f, profile[10][AminoAcid::instance().ctoi('H')], DELTA);
    EXPECT_FLOAT_EQ(1.0f, std::accumulate(profile.col_begin(122), profile.col_end(122), 0.0f));
}

TEST(CountsProfileTest, Alignment1Q7L)
{
    FILE* fin = fopen("../data/1Q7L.fas", "r");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::FASTA);
    fclose(fin);

    alignment.assign_match_columns_by_gap_rule();
    CountsProfile<AminoAcid> profile(alignment, true);

    EXPECT_NEAR(0.05f, profile[10][AminoAcid::instance().ctoi('G')], DELTA);
}

TEST(CountsProfileTest, AlignmentCelegansRefGene)
{
    FILE* fin = fopen("../data/ce_refgene.fas", "r");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fclose(fin);

    CountsProfile<Nucleotide> profile(alignment, false);

    EXPECT_FLOAT_EQ(1.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[1][Nucleotide::instance().ctoi('T')]);
}

TEST(CountsProfileTest, AddMatrixPseudocountsToProfile)
{
    FILE* fin = fopen("../data/ce_refgene.fas", "r");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fclose(fin);

    CountsProfile<Nucleotide> profile(alignment, false);

    ASSERT_FLOAT_EQ(1.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    NucleotideMatrix m(1,-1);
    MatrixPseudocounts<Nucleotide> mpc(&m);
    mpc.add_to_profile(DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

    EXPECT_NEAR(0.25f, profile[0][Nucleotide::instance().ctoi('T')], DELTA);
}

TEST(CountsProfileTest, AddMatrixPseudocountsToLogProfile)
{
    FILE* fin = fopen("../data/ce_refgene.fas", "r");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fclose(fin);

    CountsProfile<Nucleotide> profile(alignment, false);
    profile.transform_to_logspace();

    ASSERT_FLOAT_EQ(0.0f, profile[1][Nucleotide::instance().ctoi('T')]);

    NucleotideMatrix m(1,-1);
    MatrixPseudocounts<Nucleotide> mpc(&m);
    mpc.add_to_profile(DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

    profile.transform_to_linspace();
    EXPECT_NEAR(0.25f, profile[0][Nucleotide::instance().ctoi('T')], DELTA);
}

};  // cs
