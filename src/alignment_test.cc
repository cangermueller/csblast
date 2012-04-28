#include <gtest/gtest.h>

#include "cs.h"
#include "alignment-inl.h"

namespace cs {
using std::string;
using std::vector;

const string test_dir = PathCat(getenv("CS_DATA") ? getenv("CS_DATA") : "../data", "test");

TEST(AlignmentTest, DISABLED_ConstructionFromPsiInput) {
  FILE* fin = fopen("../data/26SPS9.psi", "r");
  Alignment<AA> alignment(fin, PSI_ALIGNMENT);
  fclose(fin);

  EXPECT_EQ(14, static_cast<int>(alignment.nseqs()));
  EXPECT_EQ(263, static_cast<int>(alignment.ncols()));

  alignment.AssignMatchColumnsBySequence();

  EXPECT_EQ(176, static_cast<int>(alignment.ncols()));
  EXPECT_EQ(0, alignment.ninsert());
  EXPECT_EQ(AA::kCharToInt['D'], alignment[23][0]);
  EXPECT_EQ(AA::kGap, alignment[23][1]);
}

TEST(AlignmentTest, DISABLED_ConstructionFromSequence) {
  FILE* fin = fopen("../data/d1w6ga2.seq", "r");
  Sequence<AA> query(fin);
  fclose(fin);

  Alignment<AA> ali(query);
  EXPECT_EQ(1, static_cast<int>(ali.nseqs()));
  EXPECT_EQ(query.ToString(), ali.GetSequence(0).ToString());
}

TEST(AlignmentTest, DISABLED_ReadFastaDnaAlignment) {
  FILE* fin = fopen("../data/align00.fas", "r");
  Alignment<Dna> alignment(fin, FASTA_ALIGNMENT);
  fclose(fin);

  EXPECT_EQ(Dna::kCharToInt['T'], alignment[0][0]);
  EXPECT_EQ(Dna::kCharToInt['A'], alignment[1][0]);
  EXPECT_EQ(Dna::kCharToInt['A'], alignment[2][0]);
  EXPECT_EQ(Dna::kCharToInt['T'], alignment[3][0]);
  EXPECT_EQ(Dna::kGap, alignment[10][0]);
}

TEST(AlignmentTest, DISABLED_ConstructionFromBlastHits) {
  FILE* fin = fopen("../data/blast_results_broken_hitlist.txt", "r");
  BlastHits blast_hits(fin);
  fclose(fin);

  ASSERT_EQ(5, static_cast<int>(blast_hits.nhits()));
  blast_hits.Filter(2.0);
  ASSERT_EQ(2, static_cast<int>(blast_hits.nhits()));

  Alignment<AA> ali(blast_hits);
  EXPECT_EQ(2, static_cast<int>(ali.nseqs()));
}

TEST(AlignmentTest, DISABLED_ConstructionFromBlastHitsWithMultipleHSPs) {
  FILE* fin = fopen("../data/blast_results.txt", "r");
  BlastHits blast_hits(fin);
  fclose(fin);

  ASSERT_EQ((size_t)500, blast_hits.nhits());
  blast_hits.Filter(0.2);
  ASSERT_EQ((size_t)4, blast_hits.nhits());

  Alignment<AA> ali(blast_hits, false);
  EXPECT_EQ((size_t)5, ali.nseqs());

  Alignment<AA> ali_best(blast_hits, true);
  EXPECT_EQ((size_t)4, ali_best.nseqs());
}

TEST(AlignmentTest, DISABLED_ConstructionFromA2M) {
  FILE* fin = fopen("../data/d1alx.a2m", "r");
  Alignment<AA> alignment(fin, A2M_ALIGNMENT);
  fclose(fin);

  EXPECT_EQ(AA::kGap, alignment.seq(0,27));
}

TEST(AlignmentTest, ConstructionFromA3M) {
  vector<string> alis;
  alis.push_back("101mA.a3m");
  alis.push_back("101mA_ss.a3m");
  alis.push_back("101mA_ss_whitespace.a3m");
  for (vector<string>::iterator it = alis.begin(); it != alis.end(); ++it) {
    FILE* fin = fopen(PathCat(test_dir, *it), "r");
    Alignment<AA> ali(fin, A3M_ALIGNMENT);
    fclose(fin);
    EXPECT_EQ(ali.nseqs(), 1066);
    EXPECT_EQ(ali.ncols(), 858);
  }
  alis.clear();
  alis.push_back("3nkuB.a3m");
  alis.push_back("3nkuB_ss.a3m");
  alis.push_back("3nkuB_ss_whitespace.a3m");
  for (vector<string>::iterator it = alis.begin(); it != alis.end(); ++it) {
    FILE* fin = fopen(PathCat(test_dir, *it), "r");
    Alignment<AA> ali(fin, A3M_ALIGNMENT);
    fclose(fin);
    EXPECT_EQ(ali.nseqs(), 1);
    EXPECT_EQ(ali.ncols(), 200);
  }
}

TEST(AlignmentTest, DISABLED_AssignMatchColumnsByGapRule) {
  FILE* fin = fopen("../data/MalT_diverse.fas", "r");
  Alignment<AA> ali(fin, FASTA_ALIGNMENT);
  fclose(fin);

  EXPECT_EQ(ali.ncols(), ali.nmatch() + ali.ninsert());
  EXPECT_TRUE(ali.is_match(51));
  EXPECT_TRUE(ali.is_match(52));
  EXPECT_TRUE(ali.is_match(52));

  ali.AssignMatchColumnsByGapRule(20);

  EXPECT_EQ(ali.ncols(), ali.nmatch() + ali.ninsert());
  EXPECT_TRUE(ali.is_match(51));
  EXPECT_FALSE(ali.is_match(52));
  EXPECT_FALSE(ali.is_match(52));
}

TEST(AlignmentTest, DISABLED_Merging) {
  FILE* fin = fopen("../data/MalT_slim.fas", "r");
  Alignment<AA> ali_slim(fin, FASTA_ALIGNMENT);
  fclose(fin);
  ASSERT_EQ(5, static_cast<int>(ali_slim.nseqs()));

  fin = fopen("../data/MalT_diverse.fas", "r");
  Alignment<AA> ali_diverse(fin, FASTA_ALIGNMENT);
  fclose(fin);
  ASSERT_EQ(21, static_cast<int>(ali_diverse.nseqs()));

  ali_slim.Merge(ali_diverse);
  EXPECT_EQ(23, static_cast<int>(ali_slim.nseqs()));
  EXPECT_EQ(ali_slim[180][22], ali_diverse[180][20]);
  EXPECT_EQ(ali_slim[181][22], ali_diverse[181][20]);
  EXPECT_EQ(ali_slim[182][22], ali_diverse[182][20]);
  EXPECT_EQ(ali_slim[183][22], ali_diverse[183][20]);
  EXPECT_EQ(ali_slim[184][22], ali_diverse[184][20]);
}

}  // namespace cs
