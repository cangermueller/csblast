#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "sgd.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"

namespace cs {

const double kDelta = 1e-6;

template<typename TrainingPair>
class SgdTestBlosum : public testing::Test {
 protected:
  typedef std::vector<TrainingPair> TrainingSet;

  static const size_t kNumStates    = 20;
  static const size_t kWindowLength = 13;
  static const size_t kTrainSetSize = 10000;

  SgdTestBlosum() : prior_(10.0, 10.0, 1.0, 10.0) {}

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

      trainset_.push_back(TrainingPair(seq, col));
    }
  }

  TrainingSet trainset_;
  BlosumMatrix m_;
  SgdParams params_;
  LassoDerivCrfFuncPrior<AA> prior_;
};

typedef testing::Types<TrainingSequence<AA>, TrainingProfile<AA> > MyTypes;
TYPED_TEST_CASE(SgdTestBlosum, TrainingProfile<AA>);


TYPED_TEST(SgdTestBlosum, DISABLED_InitWithGaussian) {
  GaussianCrfInit<AA> init(0.1, this->m_);
  Crf<AA> crf(this->kNumStates, this->kWindowLength, init);

  DerivCrfFunc<AA, TypeParam> func(this->trainset_, this->m_, this->prior_);
  SgdOptimizer<AA, TypeParam, TypeParam> sgd_opt(func, func, this->params_);
  double loglike = sgd_opt.Optimize(crf, stdout);

  EXPECT_NEAR(0.4308, loglike, 0.0001);
}

TYPED_TEST(SgdTestBlosum, InitWithOptimalSolution) {
  Crf<AA> crf(this->kNumStates, this->kWindowLength);
  for (size_t k = 0; k < this->kNumStates; ++k) {
    CrfState<AA> state(this->kWindowLength);
    state.bias_weight = 0.0;
    for (size_t j = 0; j < this->kWindowLength; ++j) {
      for (size_t a = 0; a < AA::kSize; ++a)
        if (j == crf.center() && k % AA::kSize == a)
          state.context_weights[j][a] = 10.0;
        else
          state.context_weights[j][a] = 0.0;
      state.context_weights[j][AA::kAny] = 0.0;
    }
    for (size_t a = 0; a < AA::kSize; ++a)
      state.pc_weights[a] = log(this->m_.r(a,k % AA::kSize));
    UpdatePseudocounts(state);
    crf.SetState(k, state);
  }

  DerivCrfFunc<AA, TypeParam> func(this->trainset_, this->m_, this->prior_);
  SgdOptimizer<AA, TypeParam, TypeParam> sgd(func, func, this->params_);
  double loglike = sgd.Optimize(crf, stdout);

  EXPECT_NEAR(0.4839, loglike, 0.0001);
}


TYPED_TEST(SgdTestBlosum, PcPrior) {
  MatrixPseudocounts<AA> pc(this->m_);
  ConstantAdmix admix(1.0);
  LassoDerivCrfFuncPrior<AA> prior = this->prior_;
  prior.sigma_pc = 0.01;
  SamplingCrfInit<AA, TypeParam> init(this->trainset_, pc, admix, this->m_); 
  Crf<AA> crf(this->kNumStates, this->kWindowLength, init);
  DerivCrfFunc<AA, TypeParam> func(this->trainset_, this->m_, prior);
  Sgd<AA, TypeParam> sgd(func, this->params_);
  SgdState<AA> s(crf);
  const size_t c = crf.center();

  // pc_weights must be equal to context_weights[center] after initialization
  for (size_t k = 0; k < crf.size(); ++k) {
    for (size_t a = 0; a < AA::kSize; ++a)
      ASSERT_EQ(s.crf[k].context_weights[c][a], s.crf[k].pc_weights[a]);
  }

  // pc_weights must be near to context_weights[center] due to the strict prior
  for (size_t r = 0; r < 10; ++r) {sgd(s); }
  for (size_t k = 0; k < crf.size(); ++k) {
    for (size_t a = 0; a < AA::kSize; ++a)
      EXPECT_NEAR(s.crf[k].context_weights[c][a], s.crf[k].pc_weights[a], 0.1);
  }
}

}  // namespace cs
