#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"

namespace cs {

const double kDelta = 0.01;

class CrfTest : public testing::Test {
 protected:

  virtual void SetUp() {
    FILE* fin = fopen("../data/crf_state.cst", "r");
    s1_.Read(fin);
    rewind(fin);
    s2_.Read(fin);
    rewind(fin);
    s3_.Read(fin);
    fclose(fin);
  }

  CrfState<AA> s1_;
  CrfState<AA> s2_;
  CrfState<AA> s3_;
};

TEST_F(CrfTest, SimpleConstruction) {
  Crf<AA> crf(3, 13);

  EXPECT_EQ((size_t)3, crf.size());
  EXPECT_EQ((size_t)13, crf.wlen());
  crf.SetState(0, s1_);
  crf.SetState(1, s2_);
  crf.SetState(2, s3_);

  EXPECT_NEAR(-1.696, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.135, crf[0].context_weights[0][AA::kCharToInt['R']], kDelta);
  EXPECT_NEAR( 0.000, crf[0].context_weights[0][AA::kCharToInt['X']], kDelta);
}

TEST_F(CrfTest, ConstructionFromSerializedCRF) {
  FILE* fin = fopen("../data/test_crf.crf", "r");
  Crf<AA> crf(fin);
  fclose(fin);

  EXPECT_NEAR(-1.696, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.135, crf[0].context_weights[0][AA::kCharToInt['R']], kDelta);
  EXPECT_NEAR( 0.000, crf[0].context_weights[0][AA::kCharToInt['X']], kDelta);
}

TEST(CrfTestInit, SamplingCrfInit) {
  BlosumMatrix m;
  MatrixPseudocounts<AA> pc(m);

  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AA> ali(fin, FASTA_ALIGNMENT);
  fclose(fin);
  ali.AssignMatchColumnsBySequence();

  CountProfile<AA> p_full(ali, true);
  p_full.counts = pc.AddTo(p_full, ConstantAdmix(0.01));
  Normalize(p_full.counts, 1.0);
  CountProfile<AA> p_win(p_full, 100, 13);

  Sequence<AA> seq_win(ali.GetSequence(0), 100, 13);
  ProfileColumn<AA> pc_col(&p_full.counts[106][0]);
  TrainingSequence<AA> train_pair(seq_win, pc_col);
  std::vector<TrainingSequence<AA> > trainset(100, train_pair);

  ConstantAdmix admix(1.0);
  SamplingCrfInit<AA, TrainingSequence<AA> > init(trainset, pc, admix, m);
  Crf<AA> crf(10, 13, init);

  EXPECT_EQ((size_t)10, crf.size());
  EXPECT_EQ((size_t)13, crf.wlen());
  EXPECT_EQ((size_t)13, crf[0].context_weights.length());

  EXPECT_NEAR(-0.59, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.39, crf[0].context_weights[1][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-0.29, crf[0].context_weights[2][AA::kCharToInt['A']], kDelta);
}

TEST(CrfTestInit, LibraryBasedCrfInit) {
  BlosumMatrix m;
  FILE* fin = fopen("../data/context_library.lib", "r");
  ContextLibrary<AA> lib(fin);
  fclose(fin);

  LibraryBasedCrfInit<AA> init(lib, m);
  Crf<AA> crf(3, 13, init);

  EXPECT_EQ((size_t)3, crf.size());
  EXPECT_EQ((size_t)13, crf.wlen());
  EXPECT_EQ((size_t)13, crf[0].context_weights.length());
  EXPECT_NEAR(-0.51, crf[0].bias_weight, kDelta);
  EXPECT_NEAR(-1.70, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.47, crf[0].context_weights[1][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.09, crf[0].context_weights[2][AA::kCharToInt['A']], kDelta);
}

TEST(CrfTestInit, GaussianCrfInit) {
  BlosumMatrix m;
  GaussianCrfInit<AA> init(0.1, m, 0);
  Crf<AA> crf(10, 13, init);

  EXPECT_EQ((size_t)10, crf.size());
  EXPECT_EQ((size_t)13, crf.wlen());
  EXPECT_EQ((size_t)13, crf[0].context_weights.length());
  EXPECT_NEAR( 0.12, crf[0].bias_weight, kDelta);
  EXPECT_NEAR(-0.03, crf[1].bias_weight, kDelta);
  EXPECT_NEAR(-0.01, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-0.02, crf[0].context_weights[1][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR(-0.20, crf[0].context_weights[2][AA::kCharToInt['A']], kDelta);
}

}  // namespace cs
