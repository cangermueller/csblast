#include <gtest/gtest.h>

#include <cstdio>

#include <string>

#include "blast_hits.h"

using std::string;

namespace cs {

TEST(BlastHitsTest, SimpleResults) {
  FILE* fin = fopen("../data/blast_results.txt", "r");
  BlastHits hits(fin);
  fclose(fin);

  EXPECT_EQ(500, hits.num_hits());
  EXPECT_EQ(151, hits.query_length());
  EXPECT_DOUBLE_EQ(264.0, hits[0].bit_score);
  EXPECT_DOUBLE_EQ(264.0, hits[0].hsps[0].bit_score);
  EXPECT_DOUBLE_EQ(9e-73, hits[0].evalue);
  EXPECT_DOUBLE_EQ(9e-73, hits[0].hsps[0].evalue);
  EXPECT_EQ(1, static_cast<int>(hits.hit(2).hsps.size()));
  EXPECT_EQ(3, static_cast<int>(hits.hit(3).hsps.size()));
  EXPECT_EQ(1, static_cast<int>(hits.hit(4).hsps.size()));
  EXPECT_EQ(string("HVLHSRHP"), string(hits.hit(3).hsps[2].query_seq.begin(),
                                       hits.hit(3).hsps[2].query_seq.end()));

  hits.Filter(0.2);

  EXPECT_EQ(4, hits.num_hits());
  EXPECT_EQ(2, static_cast<int>(hits.hit(3).hsps.size()));
}

TEST(BlastHitsTest, ResultsWithBrokenHitlist) {
  FILE* fin = fopen("../data/blast_results_broken_hitlist.txt", "r");
  BlastHits hits(fin);
  fclose(fin);

  EXPECT_EQ(5, hits.num_hits());
  EXPECT_EQ(88, hits.query_length());
}

}  // namespace cs
