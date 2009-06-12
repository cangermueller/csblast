#include <gtest/gtest.h>

#include <cstdio>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile-inl.h"
#include "shared_ptr.h"
#include "sum_product_algorithm-inl.h"
#include "utils-inl.h"

namespace cs {

const double kDelta = 0.0001;

TEST(SumProductAlgorithmTest, 1Q7L) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
  fclose(fin);
  ali.AssignMatchColumnsByGapRule();
  CountProfile<AminoAcid> profile(ali, true);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.AddPseudocountsToProfile(ConstantAdmixture(0.1f), &profile);

  const int num_cols = 5;
  CRF<AminoAcid> crf(profile.length() - num_cols + 1, num_cols);
  for (int i = 0; i < profile.length() - num_cols + 1; ++i) {
    Profile<AminoAcid> p(profile, i, num_cols);
    crf.AddState(p);
  }
  crf.InitTransitions(HomogeneousTransitionInitializerCRF<AminoAcid>());

  SumProductMatrices mat(profile.length(), crf.num_states());
  SumProductAlgorithm(crf, profile, &mat);

  EXPECT_NEAR(1.0, mat.alpha[2][0] * mat.beta[2][0], kDelta);
  EXPECT_NEAR(1.0, mat.alpha[3][1] * mat.beta[3][1], kDelta);
  EXPECT_NEAR(1.0, mat.alpha[4][2] * mat.beta[4][2], kDelta);

  EXPECT_NEAR(1.0, mat.alpha_pc[2][0] * mat.beta_pc[2][0], kDelta);
  EXPECT_NEAR(1.0, mat.alpha_pc[3][1] * mat.beta_pc[3][1], kDelta);
  EXPECT_NEAR(1.0, mat.alpha_pc[4][2] * mat.beta_pc[4][2], kDelta);

  LOG(ERROR) << mat;
}

}  // namespace cs
