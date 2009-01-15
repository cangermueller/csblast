#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"
#include "alignment_profile.h"
#include "smart_ptr.h"

TEST(AlignmentProfileTest, ConstructionFromInputStream)
{
    std::string data;
    data.append("AlignmentProfile\n");
    data.append("Profile\n");
    data.append("ncols\t4\n");
    data.append("ndim\t4\n");
    data.append("\tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("//\n");
    data.append("has_counts\t0\n");
    data.append("\tneff\n");
    data.append("1\t0\n");
    data.append("2\t0\n");
    data.append("3\t0\n");
    data.append("4\t0\n");
    data.append("//\n");
    std::istringstream ss(data);

    cs::AlignmentProfile profile(ss, cs::NucleicAcidAlphabet::instance());

    EXPECT_EQ(4, profile.ncols());
    EXPECT_EQ(4, profile.ndim());
    EXPECT_EQ(1.0f, profile(0,0));
    EXPECT_EQ(0.0f, profile(1,0));
    EXPECT_FALSE(profile.has_counts());
}

TEST(AlignmentProfileTest, ConstructionFromAlignment)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, na);

    cs::AlignmentProfile profile(alignment, true); // use position-dependent weights

    EXPECT_EQ(1.0, profile.neff(3));
    EXPECT_EQ(0.0, profile(3, na->ctoi('A')));
    EXPECT_EQ(0.0, profile(3, na->ctoi('C')));
    EXPECT_EQ(0.0, profile(3, na->ctoi('G')));
    EXPECT_EQ(1.0, profile(3, na->ctoi('T')));
}

TEST(AlignmentProfileTest, ConversionToCounts)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTA\n");
    data.append(">seq3\nGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGT\n");
    data.append(">seq4\nCGTACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACGTACACGTACG\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, na);

    cs::AlignmentProfile profile(alignment, true); // use position-dependent weights
    ASSERT_EQ(0.25, profile(3, na->ctoi('A')));

    profile.convert_to_counts();

    EXPECT_TRUE(profile.has_counts());
    EXPECT_EQ(0.25*profile.neff(3), profile(3, na->ctoi('A')));
}

TEST(AlignmentProfileTest, RealAlignment1ena)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
    std::ifstream file;
    file.open("../data/1ena.fas");
    cs::Alignment alignment(file, aa);
    file.close();
    cs::AlignmentProfile profile(alignment, true); // use position-dependent weights

    ASSERT_TRUE(profile(0, aa->ctoi('L')) > 0.8);
}

