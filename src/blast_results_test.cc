#include <gtest/gtest.h>

#include <cstdio>

#include <string>

#include "blast_results.h"

using std::string;

namespace cs {

TEST(BlastResultsTest, SimpleResults) {
  FILE* fin = fopen("../data/blast_results.txt", "r");
  BlastResults results(fin);
  fclose(fin);

  EXPECT_EQ(500, results.num_hits());
  EXPECT_EQ(3, results.hit(3).hsps.size());
  EXPECT_EQ(string("HVLHSRHP"), string(results.hit(3).hsps[2].query_seq.begin(),
                                       results.hit(3).hsps[2].query_seq.end()));
}

TEST(BlastResultsTest, ResultsWithBrokenHitlist) {
  FILE* fin = fopen("../data/blast_results_broken_hitlist.txt", "r");
  BlastResults results(fin);
  fclose(fin);

  EXPECT_EQ(5, results.num_hits());
}

}  // namespace cs
