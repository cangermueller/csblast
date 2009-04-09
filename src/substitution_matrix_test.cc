#include <gtest/gtest.h>

#include "blosum_matrix.h"
#include "nucleotide_matrix.h"

namespace cs {

const float kDelta = 0.005f;

TEST(BlosumMatrixTest, Blosum45) {
  BlosumMatrix matrix(BlosumMatrix::BLOSUM45);
  EXPECT_EQ(20, matrix.size());
  EXPECT_NEAR(1.5609, matrix.s(0,0), kDelta);
}

TEST(BlosumMatrixTest, Blosum62) {
  BlosumMatrix matrix;
  EXPECT_EQ(20, matrix.size());
  EXPECT_NEAR(1.9646, matrix.s(0,0), kDelta);
}

TEST(BlosumMatrixTest, Blosum80) {
  BlosumMatrix matrix(BlosumMatrix::BLOSUM80);
  EXPECT_EQ(20, matrix.size());
  EXPECT_NEAR(2.2549, matrix.s(0,0), kDelta);
}

TEST(NucleotideMatrixTest, Match1Mismatch1) {
  NucleotideMatrix matrix(1, -1);
  EXPECT_EQ(4, matrix.size());
  EXPECT_FLOAT_EQ(0.25, matrix.f(0));
}

TEST(NucleotideMatrixTest, Match1Mismatch2) {
  NucleotideMatrix matrix(1, -2);
  EXPECT_EQ(4, matrix.size());
  EXPECT_FLOAT_EQ(0.25, matrix.f(0));
}

}  // namespace cs
