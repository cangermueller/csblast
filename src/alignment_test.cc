#include <gtest/gtest.h>

#include <cstdio>

#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "matrix.h"
#include "nucleotide.h"

namespace cs {

TEST(AlignmentTest, ConstructionFromInputStream) {
  FILE* fin = fopen("../data/nt_alignment1.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  EXPECT_EQ(2, alignment.num_seqs());
  EXPECT_EQ(80, alignment.num_cols());
  EXPECT_EQ(Nucleotide::instance().ctoi('A'), alignment[0][0]);
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), alignment[1][1]);
  EXPECT_EQ(Nucleotide::instance().gap(), alignment[4][1]);
  EXPECT_EQ(Nucleotide::instance().endgap(), alignment[78][1]);
}

TEST(AlignmentTest, ConstructionFromBlastResults) {
  FILE* fin = fopen("../data/d1w6ga2.seq", "r");
  Sequence<AminoAcid> query(fin);
  fclose(fin);
  fin = fopen("../data/blast_results_broken_hitlist.txt", "r");
  BlastResults results(fin);
  fclose(fin);
  ASSERT_EQ(5, results.num_hits());

  Alignment<AminoAcid> ali(query, results, 2.0);
  EXPECT_EQ(3, ali.num_seqs());
}

TEST(AlignmentTest, CalculationOfGlobalWeights) {
  FILE* fin = fopen("../data/nt_alignment2.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  EXPECT_EQ(4, alignment.num_seqs());
  EXPECT_EQ(80, alignment.num_cols());

  std::vector<float> wg;
  float neff = global_weights_and_diversity(alignment, wg);

  EXPECT_EQ(4, static_cast<int>(wg.size()));
  EXPECT_FLOAT_EQ(0.25, wg[0]);
  EXPECT_FLOAT_EQ(1.0, neff);
}

TEST(AlignmentTest, CalculationOfPositionSpecificWeights) {
  FILE* fin = fopen("../data/nt_alignment3.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  EXPECT_EQ(4, alignment.num_seqs());
  EXPECT_EQ(80, alignment.num_cols());

  matrix<float> w;
  position_specific_weights_and_diversity(alignment, w);

  EXPECT_FLOAT_EQ(0.5, w[0][0]);
}

TEST(AlignmentTest, ConstructionFromCelegansRefGene) {
  FILE* fin = fopen("../data/ce_refgene.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  EXPECT_EQ(Nucleotide::instance().ctoi('C'), alignment[0][0]);
}

TEST(AlignmentTest, RemoveColumnsWithGapInFirst) {
  FILE* fin = fopen("../data/nt_alignment4.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  ASSERT_EQ(2, alignment.num_seqs());
  ASSERT_EQ(80, alignment.num_cols());

  alignment.assign_match_columns_by_sequence(0);

  EXPECT_EQ(76, alignment.num_match_cols());
  EXPECT_EQ(Nucleotide::instance().gap(), alignment.seq(0,1));
  EXPECT_EQ(Nucleotide::instance().ctoi('G'), alignment[1][1]);
}

TEST(AlignmentTest, RemoveColumnsByGapRule) {
  FILE* fin = fopen("../data/nt_alignment5.fas", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::FASTA);
  fclose(fin);

  ASSERT_EQ(3, alignment.num_seqs());
  ASSERT_EQ(80, alignment.num_cols());

  alignment.assign_match_columns_by_gap_rule();

  EXPECT_EQ(76, alignment.num_match_cols());
  EXPECT_EQ(Nucleotide::instance().ctoi('G'), alignment[4][0]);
}

TEST(AlignmentTest, ConstructionFromA2M) {
  FILE* fin = fopen("../data/d1alx.a2m", "r");
  Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::A2M);
  fclose(fin);

  EXPECT_EQ(AminoAcid::instance().gap(), alignment.seq(0,27));
}

TEST(AlignmentTest, ConstructionFromA3M) {
  FILE* fin = fopen("../data/d1alx.a3m", "r");
  Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::A3M);
  fclose(fin);

  EXPECT_EQ(AminoAcid::instance().gap(), alignment.seq(0,27));
}

TEST(AlignmentTest, RemoveInsertColumns) {
  FILE* fin = fopen("../data/nt_alignment6.a2m", "r");
  Alignment<Nucleotide> alignment(fin, Alignment<Nucleotide>::A2M);
  fclose(fin);

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

}  // namespace cs
