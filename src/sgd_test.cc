#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "sgd.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"

namespace cs {

class SgdTestBlosum : public testing::Test {
 protected:
  typedef std::vector<TrainingSequence<AA> > TrainingSet;

  static const size_t kNumStates    = 20;
  static const size_t kWindowLength = 13;
  static const size_t kTrainSetSize = 10000;

  SgdTestBlosum() {}

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
  }

  TrainingSet trainset_;
  BlosumMatrix m_;
  SgdParams params_;
};

TEST_F(SgdTestBlosum, InitWithGaussian) {
  GaussianCrfInit<AA> init(0.1, m_);
  Crf<AA> crf(kNumStates, kWindowLength, init);

  DerivCrfFunc<AA, TrainingSequence<AA> > func(trainset_, m_, 1.0, 1.0, 10.0);
  SgdOptimizer<AA, TrainingSequence<AA> > sgd_opt(func, func, params_);
  double loglike = sgd_opt.Optimize(crf, stdout);

  EXPECT_NEAR(0.4308, loglike, 0.0001);
}

TEST_F(SgdTestBlosum, InitWithOptimalSolution) {
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

  DerivCrfFunc<AA, TrainingSequence<AA> > func(trainset_, m_, 1.0, 1.0, 10.0);
  SgdOptimizer<AA, TrainingSequence<AA> > sgd(func, func, params_);
  double loglike = sgd.Optimize(crf, stdout);

  EXPECT_NEAR(0.4839, loglike, 0.0001);
}

}  // namespace cs
