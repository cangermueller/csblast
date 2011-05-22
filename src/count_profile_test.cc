#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"

namespace cs {

const double kDelta = 0.01;

TEST(CountProfileTest, DISABLED_AlignmentBpdSWithConsensus) {
  FILE* fin = fopen("../data/BpdS.fas", "r");
  Alignment<AA> alignment(fin, FASTA_ALIGNMENT);
  fclose(fin);

  alignment.AssignMatchColumnsByGapRule();
  CountProfile<AA> profile(alignment, true);

  EXPECT_NEAR(0.0, profile.counts[10][AA::kCharToInt['H']], kDelta);

  BlosumMatrix sm;
  Sequence<AA> cons = ConsensusSequence(profile, sm);

  LOG(ERROR) << cons;
  LOG(ERROR) << profile;
}

TEST(CountProfileTest, DISABLED_Alignment1Q7L) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AA> alignment(fin, FASTA_ALIGNMENT);
  fclose(fin);

  alignment.AssignMatchColumnsByGapRule();
  CountProfile<AA> profile(alignment, true);

  Normalize(profile.counts, 1.0);  // convert to relative frequencies
  EXPECT_NEAR(0.05, profile.counts[10][AA::kCharToInt['G']], kDelta);
}

TEST(CountProfileTest, DISABLED_ConstructionFromSerializedCountProfile) {
  FILE* fin = fopen("../data/dna_count_profile.prf", "r");
  CountProfile<Dna> cp(fin);
  fclose(fin);

  EXPECT_EQ("Count-Profile 1", cp.name);
  EXPECT_EQ(4, static_cast<int>(cp.length()));

  LOG(ERROR) << cp;
}

TEST(ProfileTest, Insert) {
  const size_t len = 100;
  Profile<AA> dest(len);
  for (size_t r = 0; r < 10; ++r) {
    Profile<AA> src(len);
    for (size_t i = 0; i < len; ++i)
      for (size_t a = 0; a < AA::kSize; ++a) 
        src[i][a] = rand();
    size_t idx = static_cast<size_t>(rand() % len);
    dest.Insert(idx, src);
    for (size_t i = idx; i < len; ++i) 
      for (size_t a = 0; a < AA::kSize; ++a)
        EXPECT_EQ(src[i - idx][a], dest[i][a]);
  }
}



};  // namespace cs
