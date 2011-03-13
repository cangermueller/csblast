#include <gtest/gtest.h>

#include "cs.h"
#include "context_profile-inl.h"

namespace cs {

const double kDelta = 0.001;

TEST(ContextProfileTest, ConstructionAndLinLogTransform) {
  FILE* fin = fopen("../data/context_profile.prf", "r");
  ContextProfile<AA> cp(fin);
  fclose(fin);

  EXPECT_EQ("1", cp.name);

  TransformToLin(cp);
  EXPECT_FALSE(cp.is_log);
  EXPECT_NEAR(0.0136, cp.probs[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(0.0483, cp.probs[2][AA::kCharToInt['N']], kDelta);
  EXPECT_NEAR(1.0000, cp.probs[0][AA::kCharToInt['X']], kDelta);

  TransformToLog(cp);
  EXPECT_TRUE(cp.is_log);
  EXPECT_NEAR(-4.2975, cp.probs[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-3.0333, cp.probs[2][AA::kCharToInt['N']], kDelta);
  EXPECT_NEAR( 0.0, cp.probs[0][AA::kCharToInt['X']], kDelta);

  FILE* fout = fopen("/tmp/context_profile.prf", "w");
  cp.Write(fout);
  fclose(fout);
}

};  // namespace cs
