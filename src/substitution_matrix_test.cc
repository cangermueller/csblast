#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "gonnet_matrix.h"
#include "profile-inl.h"
#include "pdf_writer-inl.h"
#include "tamura_nei_matrix.h"

namespace cs {

const double kDelta = 0.01;

TEST(BlosumMatrixTest, Blosum45) {
    BlosumMatrix matrix(BLOSUM45);
    EXPECT_NEAR(4.68806, matrix.s(0,0), kDelta);
    EXPECT_NEAR(14.2328, Neff<AA>(matrix), kDelta);
}

TEST(GonnetMatrixTest, GonnetMatrix) {
    GonnetMatrix matrix;
    EXPECT_NEAR( 0.5544, matrix.s(0,0), kDelta);
    EXPECT_NEAR(15.4168, Neff<AA>(matrix), kDelta);

    Profile<AA> profile(AA::kSize);
    for (size_t a = 0; a < AA::kSize; ++a)
        for (size_t b = 0; b < AA::kSize; ++b)
            profile[a][b] = matrix.r(b,a);
}

TEST(BlosumMatrixTest, Blosum62) {
    BlosumMatrix matrix(BLOSUM45);
    EXPECT_NEAR(3.93474, matrix.s(0,0), kDelta);
    EXPECT_NEAR(11.4868, Neff<AA>(matrix), kDelta);

    Profile<AA> profile(AA::kSize);
    for (size_t a = 0; a < AA::kSize; ++a)
        for (size_t b = 0; b < AA::kSize; ++b)
            profile[a][b] = matrix.r(b,a);
}

TEST(BlosumMatrixTest, Blosum80) {
    BlosumMatrix matrix(BLOSUM80);
    EXPECT_NEAR(4.51468, matrix.s(0,0), kDelta);
    EXPECT_NEAR(9.5318, Neff<AA>(matrix), kDelta);
}

TEST(NucleotideMatrixTest, TamuraNeiDefault) {
    TamuraNeiMatrix matrix;

    EXPECT_NEAR(0.90, matrix.s(0,0), kDelta);
    EXPECT_NEAR(-1.11, matrix.s(0,1), kDelta);
    EXPECT_NEAR(-0.14, matrix.s(0,2), kDelta);
    EXPECT_NEAR(0.25, matrix.p(0), kDelta);
    EXPECT_NEAR(0.25, matrix.p(1), kDelta);
    EXPECT_NEAR(0.25, matrix.p(2), kDelta);
    EXPECT_NEAR(0.25, matrix.p(3), kDelta);
}

TEST(NucleotideMatrixTest, TamuraNeiNonDefault) {
    double my_bg_freqs[] = { 0.4, 0.3, 0.2, 0.1 };
    TamuraNeiMatrix matrix(0.4, my_bg_freqs, 1.3, 1.3, 1.0);

    EXPECT_NEAR(0.57, matrix.s(0,0), kDelta);
    EXPECT_NEAR(-1.10, matrix.s(0,1), kDelta);
    EXPECT_NEAR(-0.25, matrix.s(0,2), kDelta);
    EXPECT_NEAR(0.4, matrix.p(0), kDelta);
    EXPECT_NEAR(0.3, matrix.p(1), kDelta);
    EXPECT_NEAR(0.2, matrix.p(2), kDelta);
    EXPECT_NEAR(0.1, matrix.p(3), kDelta);
}

}  // namespace cs
