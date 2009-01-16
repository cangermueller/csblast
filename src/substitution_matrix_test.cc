#include <gtest/gtest.h>

#include "blosum_matrix.h"

const float kDelta = 0.005f;

TEST(BlosumMatrixTest, Blosum45)
{
    cs::BlosumMatrix matrix(cs::BlosumMatrix::BLOSUM45);
    EXPECT_EQ(20, matrix.size());
    EXPECT_NEAR(1.5609, matrix.s(0,0), kDelta);
}

TEST(BlosumMatrixTest, Blosum62)
{
    cs::BlosumMatrix matrix;
    EXPECT_EQ(20, matrix.size());
    EXPECT_NEAR(1.9646, matrix.s(0,0), kDelta);
}

TEST(BlosumMatrixTest, Blosum80)
{
    cs::BlosumMatrix matrix(cs::BlosumMatrix::BLOSUM80);
    EXPECT_EQ(20, matrix.size());
    EXPECT_NEAR(2.2549, matrix.s(0,0), kDelta);
}
