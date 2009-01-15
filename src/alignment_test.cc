#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "alignment.h"
#include "matrix.h"
#include "nucleic_acid_alphabet.h"

TEST(AlignmentTest, ConstructionFromInputStream)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, na);

    EXPECT_EQ(2, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());
    EXPECT_EQ(na->ctoi('A'), alignment(0,0));
    EXPECT_EQ(na->ctoi('C'), alignment(1,1));
    EXPECT_TRUE(alignment.gap(1,4));
    EXPECT_TRUE(alignment.endgap(1,78));
}

TEST(AlignmentTest, CalculationOfGlobalWeights)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, na);

    EXPECT_EQ(4, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());

    std::pair<std::vector<float>, float> wg_neff = cs::global_weights_and_diversity(alignment);

    EXPECT_EQ(4, static_cast<int>(wg_neff.first.size()));
    EXPECT_FLOAT_EQ(0.25, wg_neff.first[0]);
    EXPECT_FLOAT_EQ(1.0, wg_neff.second);
}

TEST(AlignmentTest, CalculationOfPositionDependentWeights)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTTACGTACACGTACGTACACGTACGTA\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\n----GTACGTACACGTACGTACACGTACGT\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\n----CGTACGTACACGTACGTACACGTACG\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, na);

    EXPECT_EQ(4, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());

    std::pair< Matrix<float>, std::vector<float> > wi_neff = cs::position_dependent_weights_and_diversity(alignment);

    EXPECT_FLOAT_EQ(0.5, wi_neff.first(0,0));
}

TEST(AlignmentTest, ConstructionFromCelegansRefGene)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::ifstream fin("../data/ce_refgene.fas");
    cs::Alignment alignment(fin, na);
    fin.close();

    EXPECT_EQ(na->ctoi('C'), alignment(0,0));
}


