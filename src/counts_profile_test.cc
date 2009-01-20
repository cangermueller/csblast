#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "amino_acid_alphabet.h"
#include "nucleotide_alphabet.h"
#include "counts_profile.h"
#include "smart_ptr.h"

const float kDelta = 0.01f;

TEST(CountsProfileTest, ConstructionFromInputStream)
{
    std::string data;
    data.append("CountsProfile\n");
    data.append("ncols\t4\n");
    data.append("ndim\t4\n");
    data.append("has_counts\t0\n");
    data.append(" \tA\tC\tG\tT\tneff\n");
    data.append("1\t0\t*\t*\t*\t0\n");
    data.append("2\t*\t0\t*\t*\t0\n");
    data.append("3\t*\t*\t0\t*\t0\n");
    data.append("4\t*\t*\t*\t0\t0\n");
    data.append("//\n");
    std::istringstream ss(data);

    cs::CountsProfile profile(ss, cs::NucleotideAlphabet::instance());

    EXPECT_EQ(4, profile.ncols());
    EXPECT_EQ(4, profile.ndim());
    EXPECT_FLOAT_EQ(1.0f, profile(0,0));
    EXPECT_FLOAT_EQ(0.0f, profile(1,0));
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
    cs::Alignment alignment(ss, na);

    cs::CountsProfile profile(alignment, true); // use position-dependent weights

    EXPECT_FLOAT_EQ(1.0, profile.neff(3));
    EXPECT_FLOAT_EQ(0.0, profile(3, na->ctoi('A')));
    EXPECT_FLOAT_EQ(0.0, profile(3, na->ctoi('C')));
    EXPECT_FLOAT_EQ(0.0, profile(3, na->ctoi('G')));
    EXPECT_FLOAT_EQ(1.0, profile(3, na->ctoi('T')));
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
    cs::Alignment alignment(ss, na);

    cs::CountsProfile profile(alignment, true);
    ASSERT_EQ(0.25, profile(3, na->ctoi('A')));

    profile.convert_to_counts();

    EXPECT_TRUE(profile.has_counts());
    EXPECT_FLOAT_EQ(0.25*profile.neff(3), profile(3, na->ctoi('A')));
}

TEST(CountsProfileTest, DISABLED_AlignmentBpdS)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    std::ifstream fin("../data/BpdS.fas");
    cs::Alignment alignment(fin, aa);
    fin.close();
    cs::CountsProfile profile(alignment, true);

    EXPECT_NEAR(0.92, profile(122, aa->ctoi('H')), kDelta);
}

TEST(CountsProfileTest, Alignment1Q7L)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    std::ifstream fin("../data/1Q7L.fas");
    cs::Alignment alignment(fin, aa);
    fin.close();
    cs::CountsProfile profile(alignment, true);

    EXPECT_NEAR(0.61, profile(579, aa->ctoi('G')), kDelta);
}

TEST(CountsProfileTest, AlignmentCelegansRefGene)
{
    cs::NucleotideAlphabet* aa = cs::NucleotideAlphabet::instance();
    std::ifstream fin("../data/ce_refgene.fas");
    cs::Alignment alignment(fin, aa);
    fin.close();
    cs::CountsProfile profile(alignment, false);

    EXPECT_FLOAT_EQ(1.0, profile(1, aa->ctoi('T')));
}
