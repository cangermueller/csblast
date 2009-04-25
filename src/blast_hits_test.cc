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
  EXPECT_EQ(3, hits.hit(3).hsps.size());
  EXPECT_EQ(string("HVLHSRHP"), string(hits.hit(3).hsps[2].query_seq.begin(),
                                       hits.hit(3).hsps[2].query_seq.end()));
}

TEST(BlastHitsTest, ResultsWithBrokenHitlist) {
  FILE* fin = fopen("../data/blast_results_broken_hitlist.txt", "r");
  BlastHits hits(fin);
  fclose(fin);

  EXPECT_EQ(5, hits.num_hits());
}

}  // namespace cs
