#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "hmc.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"

namespace cs {

const double kDelta = 0.1;

TEST(HmcTest, Likelihood) {
  FILE* fp = fopen("../data/nr30_neff2.5_1psi_W13_N100.tsq", "r");
  std::vector<TrainingSequence<AA> > trainset;
  ReadAll(fp, trainset);
  fclose(fp);

  BlosumMatrix m;
  GaussianCrfInit<AA> init(0.0001, m, 0);
  Crf<AA> crf(10, 13, init);

  CrfFunc<AA, TrainingSequence<AA> > llfunc(trainset, m);
  double loglike = llfunc(crf);
  EXPECT_NEAR(0.0, loglike, 0.01);
}

class HmcTestBlosum : public testing::Test {
 protected:
  typedef std::vector<TrainingSequence<AA> > TrainingSet;

  static const size_t kNumStates    = 20;
  static const size_t kWindowLength = 13;
  static const size_t kTrainSetSize = 10000;

  HmcTestBlosum() {}

  virtual void SetUp() {
    srand(0);
    const size_t c = (kWindowLength - 1) / 2;
    for (size_t n = 0; n < kTrainSetSize; ++n) {
      // Generate input sequence by sampling from Blosum background freqs
      Sequence<AA> seq(kWindowLength);
      seq.set_header(strprintf("seq%u", n + 1));
      for (size_t i = 0; i < kWindowLength; ++i) {
        double r = rand() / double(RAND_MAX);
        double sum = 0.0f;
        for (size_t a = 0; a < AA::kSize; ++a) {
          sum += m_.p(a);
          if (r <= sum) { seq[i] = a; break; }
        }
      }
      // Set pseudocounts to conditional probabilities P(a|b) given aa 'b' in
      // central position of input sequence window.
      ProfileColumn<AA> col;
      for (size_t a = 0; a < AA::kSize; ++a)
        col[a] = m_.r(a, seq[c]);

      trainset_.push_back(TrainingSequence<AA>(seq, col));
    }

    // HMC parameter settings
    params_.nblocks    = 1;
    params_.lfsteps    = 100;
    params_.epsilon    = 0.01;
    params_.temp       = 1.0;
    params_.sgd_epochs = 30;
  }

  TrainingSet trainset_;
  BlosumMatrix m_;
  HmcParams params_;
};

TEST_F(HmcTestBlosum, GaussianInitLeapfrog) {
  GaussianCrfInit<AA> init(0.01, m_);
  Crf<AA> crf(kNumStates, kWindowLength, init);

  HmcState<AA> state(crf);
  CrfFunc<AA, TrainingSequence<AA> > llfunc(trainset_, m_);
  DerivCrfFunc<AA, TrainingSequence<AA> > gradfunc(trainset_, m_, 1., 1., 10.);
  LeapfrogProposal<AA, TrainingSequence<AA> > propose(gradfunc, params_);
  HmcStep(3, state, propose, llfunc, params_.seed);
  EXPECT_NEAR(3411.12, state.loglike, kDelta);
}

TEST_F(HmcTestBlosum, GaussianInitBasinHopping) {
  GaussianCrfInit<AA> init(0.01, m_);
  Crf<AA> crf(kNumStates, kWindowLength, init);

  HmcState<AA> state(crf);
  CrfFunc<AA, TrainingSequence<AA> > llfunc(trainset_, m_);
  DerivCrfFunc<AA, TrainingSequence<AA> > gradfunc(trainset_, m_, 1., 1., 10.);
  SgdParams sgd_params;
  BasinHoppingProposal<AA, TrainingSequence<AA> > propose(gradfunc,
                                                          params_,
                                                          sgd_params);
  HmcStep(3, state, propose, llfunc, params_.seed);
  EXPECT_NEAR(4598.56, state.loglike, kDelta);
}

TEST_F(HmcTestBlosum, InitWithOptimalSolution) {
  Crf<AA> crf(kNumStates, kWindowLength);
  for (size_t k = 0; k < kNumStates; ++k) {
    CrfState<AA> state(kWindowLength);
    state.bias_weight = 0.0;
    for (size_t j = 0; j < kWindowLength; ++j) {
      for (size_t a = 0; a < AA::kSize; ++a)
        if (j == crf.center() && k % AA::kSize == a)
          state.context_weights[j][a] = 10.0;
        else
          state.context_weights[j][a] = 0.0;
      state.context_weights[j][AA::kAny] = 0.0;
    }
    for (size_t a = 0; a < AA::kSize; ++a)
      state.pc_weights[a] = log(m_.r(a,k % AA::kSize));
    UpdatePseudocounts(state);
    crf.SetState(k, state);
  }

  CrfFunc<AA, TrainingSequence<AA> > loglike(trainset_, m_);
  EXPECT_NEAR(4838.52, loglike(crf), kDelta);
}

TEST_F(HmcTestBlosum, InitWithAlternativeSolution) {
  Crf<AA> crf(kNumStates, kWindowLength);
  for (size_t k = 0; k < kNumStates; ++k) {
    CrfState<AA> state(kWindowLength);
    state.bias_weight = log(m_.p(k));;
    for (size_t j = 0; j < kWindowLength; ++j) {
      for (size_t a = 0; a < AA::kSize; ++a)
        if (j == crf.center())
          state.context_weights[j][a] = log(m_.r(a,k % AA::kSize) / m_.p(a));
        else
          state.context_weights[j][a] = 0.0;
      state.context_weights[j][AA::kAny] = 0.0;
    }
    for (size_t a = 0; a < AA::kSize; ++a)
      state.pc_weights[a] = k % AA::kSize == a ? 10.0 : 0;
    UpdatePseudocounts(state);
    crf.SetState(k, state);
  }

  CrfFunc<AA, TrainingSequence<AA> > loglike(trainset_, m_);
  EXPECT_NEAR(4838.52, loglike(crf), kDelta);
}

}  // namespace cs
