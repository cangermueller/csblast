#include <gtest/gtest.h>

#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "alignment.h"
#include "amino_acid.h"
#include "matrix.h"
#include "nucleotide.h"

namespace cs
{

TEST(AlignmentTest, ConstructionFromInputStream)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    EXPECT_EQ(2, alignment.num_seqs());
    EXPECT_EQ(80, alignment.num_cols());
    EXPECT_EQ(Nucleotide::instance().ctoi('A'), alignment[0][0]);
    EXPECT_EQ(Nucleotide::instance().ctoi('C'), alignment[1][1]);
    EXPECT_EQ(Nucleotide::instance().gap(), alignment[4][1]);
    EXPECT_EQ(Nucleotide::instance().endgap(), alignment[78][1]);
}

TEST(AlignmentTest, CalculationOfGlobalWeights)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    EXPECT_EQ(4, alignment.num_seqs());
    EXPECT_EQ(80, alignment.num_cols());

    std::vector<float> wg;
    float neff = global_weights_and_diversity(alignment, wg);

    EXPECT_EQ(4, static_cast<int>(wg.size()));
    EXPECT_FLOAT_EQ(0.25, wg[0]);
    EXPECT_FLOAT_EQ(1.0, neff);
}

TEST(AlignmentTest, CalculationOfPositionSpecificWeights)
{
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTTACGTACACGTACGTACACGTACGTA\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq3\n----GTACGTACACGTACGTACACGTACGT\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq4\n----CGTACGTACACGTACGTACACGTACG\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    EXPECT_EQ(4, alignment.num_seqs());
    EXPECT_EQ(80, alignment.num_cols());

    matrix<float> w;
    position_specific_weights_and_diversity(alignment, w);

    EXPECT_FLOAT_EQ(0.5, w[0][0]);
}

TEST(AlignmentTest, ConstructionFromCelegansRefGene)
{
    std::ifstream fin("../data/ce_refgene.fas");
    Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
    fin.close();

    EXPECT_EQ(Nucleotide::instance().ctoi('C'), alignment[0][0]);
}

TEST(AlignmentTest, RemoveColumnsWithGapInFirst)
{
    std::string data;
    data.append(">seq1\nA-GTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    ASSERT_EQ(2, alignment.num_seqs());
    ASSERT_EQ(80, alignment.num_cols());

    alignment.assign_match_columns_by_sequence(0);

    EXPECT_EQ(76, alignment.num_match_cols());
    EXPECT_EQ(Nucleotide::instance().gap(), alignment.seq(0,1));
    EXPECT_EQ(Nucleotide::instance().ctoi('G'), alignment[1][1]);
}

TEST(AlignmentTest, RemoveColumnsByGapRule)
{
    std::string data;
    data.append(">seq1\nA-GTACGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT--GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    data.append(">seq3\nACGT--GTACGTACGTACACGTACGTACAC\nACGTACCACGTACGTACACGGTATACGTAC\nACGTACGTAGTACGT-----\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::FASTA);

    ASSERT_EQ(3, alignment.num_seqs());
    ASSERT_EQ(80, alignment.num_cols());

    alignment.assign_match_columns_by_gap_rule();

    EXPECT_EQ(76, alignment.num_match_cols());
    EXPECT_EQ(Nucleotide::instance().ctoi('G'), alignment[4][0]);
}

TEST(AlignmentTest, ConstructionFromA2M)
{
    std::ifstream fin("../data/d1alx.a2m");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::A2M);
    fin.close();

    EXPECT_EQ(AminoAcid::instance().gap(), alignment.seq(0,27));
}

TEST(AlignmentTest, ConstructionFromA3M)
{
    std::ifstream fin("../data/d1alx.a3m");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::A3M);
    fin.close();

    EXPECT_EQ(AminoAcid::instance().gap(), alignment.seq(0,27));
}

TEST(AlignmentTest, RemoveInsertColumns)
{
    std::string data;
    data.append(">seq1\nA-GTacGTACACGTACGTACACGTACGTAC\nACGTACGTACA---ACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGT..GTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTA---GTACGT--\n");
    data.append(">seq3\nACGT..GTACGTACGTACACGTACGTACAC\nACGTACCACGTACGTACACGGTATACGTAC\nACGTACGTAGTACGT-----\n");
    std::istringstream ss(data);
    Alignment<Nucleotide> alignment(ss, Alignment<Nucleotide>::A2M);

    ASSERT_EQ(3, alignment.num_seqs());
    ASSERT_EQ(80, alignment.num_cols());
    ASSERT_TRUE(alignment.match_column(3));
    ASSERT_FALSE(alignment.match_column(4));
    ASSERT_EQ('A', alignment.chr(0,4));

    alignment.remove_insert_columns();

    ASSERT_EQ(3, alignment.num_seqs());
    ASSERT_EQ(78, alignment.num_cols());
    ASSERT_TRUE(alignment.match_column(3));
    ASSERT_TRUE(alignment.match_column(4));
    ASSERT_EQ('G', alignment.chr(0,4));
}

}  // cs
