#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "context_profile-inl.h"
#include "crf_state-inl.h"

namespace cs {

const double kDelta = 0.01;

TEST(CrfStateTest, ConstructionFromContextProfile) {
  BlosumMatrix m;
  FILE* fin = fopen("../data/context_profile.prf", "r");
  ContextProfile<AA> cp(fin);
  ASSERT_TRUE(cp.is_log);
  fclose(fin);

  CrfState<AA> state1(cp, m);
  EXPECT_NEAR(-1.70, state1.context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.47, state1.context_weights[1][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.09, state1.context_weights[2][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-0.74, state1.context_weights[3][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-1.02, state1.context_weights[4][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.00, state1.context_weights[0][AA::kCharToInt['X']], kDelta);

  TransformToLin(cp);

  ASSERT_FALSE(cp.is_log);
  ASSERT_NEAR(0.0136, cp.probs[0][AA::kCharToInt['A']], kDelta);
  ASSERT_NEAR(0.0483, cp.probs[2][AA::kCharToInt['N']], kDelta);
  ASSERT_NEAR(1.0000, cp.probs[0][AA::kCharToInt['X']], kDelta);

  CrfState<AA> state2(cp, m);
  EXPECT_NEAR(-1.70, state2.context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.47, state2.context_weights[1][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.09, state2.context_weights[2][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-0.74, state2.context_weights[3][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-1.02, state2.context_weights[4][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.00, state2.context_weights[0][AA::kCharToInt['X']], kDelta);
}

};  // namespace cs
