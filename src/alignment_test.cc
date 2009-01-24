#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "alignment.h"
#include "matrix.h"
#include "nucleotide_alphabet.h"

TEST(AlignmentTest, ConstructionFromInputStream)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA_INPUT, na);

    EXPECT_EQ(2, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());
    EXPECT_EQ(na->ctoi('A'), alignment[0][0]);
    EXPECT_EQ(na->ctoi('C'), alignment[1][1]);
    EXPECT_EQ(na->gap(), alignment[4][1]);
    EXPECT_EQ(na->endgap(), alignment[78][1]);
}

TEST(AlignmentTest, CalculationOfGlobalWeights)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA_INPUT, na);

    EXPECT_EQ(4, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());

    std::pair<std::vector<float>, float> wg_neff = cs::global_weights_and_diversity(alignment);

    EXPECT_EQ(4, static_cast<int>(wg_neff.first.size()));
    EXPECT_FLOAT_EQ(0.25, wg_neff.first[0]);
    EXPECT_FLOAT_EQ(1.0, wg_neff.second);
}

TEST(AlignmentTest, CalculationOfPositionSpecificWeights)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTTACGTACACGTACGTACACGTACGTA\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\n----GTACGTACACGTACGTACACGTACGT\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\n----CGTACGTACACGTACGTACACGTACG\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA_INPUT, na);

    EXPECT_EQ(4, alignment.nseqs());
    EXPECT_EQ(80, alignment.ncols());

    std::pair< matrix<float>, std::vector<float> > wi_neff = cs::position_specific_weights_and_diversity(alignment);

    EXPECT_FLOAT_EQ(0.5, wi_neff.first[0][0]);
}

TEST(AlignmentTest, ConstructionFromCelegansRefGene)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::ifstream fin("../data/ce_refgene.fas");
    cs::Alignment alignment(fin, cs::Alignment::FASTA_INPUT, na);
    fin.close();

    EXPECT_EQ(na->ctoi('C'), alignment[0][0]);
}

TEST(AlignmentTest, RemoveColumnsWithGapInFirst)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nA-GTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA_INPUT, na);

    ASSERT_EQ(2, alignment.nseqs());
    ASSERT_EQ(80, alignment.ncols());

    alignment.assign_match_columns_by_first_sequence();

    EXPECT_EQ(76, alignment.nmatch());
    EXPECT_EQ(na->gap(), alignment.seq(0,1));
    EXPECT_EQ(na->ctoi('G'), alignment[1][1]);
}

TEST(AlignmentTest, RemoveColumnsByGapRule)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nA-GTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    data.append(">seq3\nACGT--GTACGTACGTACACGTACGTACAC\nACGTACCACGTACGTACACGGTATACGTAC\nACGTACGTAGTACGT-----\n");
    std::istringstream ss(data);
    cs::Alignment alignment(ss, cs::Alignment::FASTA_INPUT, na);

    ASSERT_EQ(3, alignment.nseqs());
    ASSERT_EQ(80, alignment.ncols());

    alignment.assign_match_columns_by_gap_rule();

    EXPECT_EQ(76, alignment.nmatch());
    EXPECT_EQ(na->ctoi('G'), alignment[4][0]);
}


