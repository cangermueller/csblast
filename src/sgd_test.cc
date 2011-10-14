#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "sgd.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"

namespace cs {

const double kDelta = 1e-6;
const std::string test_dir = PathCat(getenv("CS_DATA") ? getenv("CS_DATA") : "../data", "test");

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
TYPED_TEST_CASE(SgdTestBlosum, TrainingSequence<AA>);


TYPED_TEST(SgdTestBlosum, DISABLED_InitWithGaussian) {
  GaussianCrfInit<AA> init(0.1, this->m_);
  Crf<AA> crf(this->kNumStates, this->kWindowLength, init);

  DerivCrfFunc<AA, TypeParam> func(this->trainset_, this->m_, this->prior_);
  SgdOptimizer<AA, TypeParam, TypeParam> sgd_opt(func, func, this->params_);
  double loglike = sgd_opt.Optimize(crf, stdout);

  EXPECT_NEAR(0.4308, loglike, 0.0001);
}

TYPED_TEST(SgdTestBlosum, DISABLED_InitWithOptimalSolution) {
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


TYPED_TEST(SgdTestBlosum, DISABLED_PcPrior) {
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

class SgdTest : public testing::Test {

  virtual void SetUp() {
    FILE* fin = fopen(PathCat(test_dir, "sgd_test.tsq"), "r");
    for (size_t n = 0; n < kSamples; ++n)
      tset_.push_back(TrainingSequence<AA>(fin));
    fclose(fin);
    vset_ = tset_;
    printf("%zu training samples read!\n", tset_.size());

    MatrixPseudocounts<AA> pc(sm_);
    ConstantAdmix admix(0.75);
    SamplingCrfInit<AA, TrainingSequence<AA> > init(tset_, pc, admix, sm_, 0, 1.6, 0.85);
    crf_.reset(new Crf<AA>(kStates, kWlen, init));
  }

  protected:
  typedef vector<TrainingSequence<AA> > TrainingSet;

  static const size_t kSamples = 10000;
  static const size_t kStates  = 20;
  static const size_t kWlen    = 13;

  TrainingSet tset_;
  TrainingSet vset_;
  BlosumMatrix sm_;
  scoped_ptr<Crf<AA> > crf_;
};

TEST_F(SgdTest, SigmaPenalty) {
  const double kSigmaContext   = 0.6;
  const double kSigmaDecay     = 0.8;
  const double kContextPenalty = 0.0;

  Vector<double> column_weights(kWlen);
  size_t c = crf_->center();
  column_weights[c] = kSigmaContext;
  for (size_t i = 1; i < c; ++i) {
    double weight = column_weights[c + i - 1] * kSigmaDecay;
    column_weights[c - i] = weight;
    column_weights[c + i] = weight;
  }

  LassoDerivCrfFuncPrior<AA> prior(1.0, kSigmaContext, kSigmaDecay, 1.0, kContextPenalty);
  DerivCrfFunc<AA, TrainingSequence<AA> > tf(tset_, sm_, prior);
  CrfFunc<AA, TrainingSequence<AA> > vf(tset_, sm_);
  
  SgdParams params;
  params.eta_mode = SgdParams::ETA_MODE_FUNC;
  params.eta_init = 0.1;
  params.eta_decay = 2.0;
  params.nblocks = 1000;

  Sgd<AA, TrainingSequence<AA> > sgd(tf, params);
  SgdState<AA> s(*crf_);

  printf("#    %5s %5s %5s %9s %9s\n", "MIN", "AVG", "MAX", "PRIOR", "LLTRAIN");
  for (size_t i = 0; i < 100; ++i) {
    Crf<AA>& crf = s.crf;
    double min = DBL_MAX;
    double max = -DBL_MAX;
    double avg = 0.0;
    size_t n   = 0;
    size_t c   = crf.center();
    for (size_t k = 0; k < crf_->size(); ++k) {
      for (size_t i = c; i <=c; ++i) {
        double sum = 0.0;
        for (size_t a = 0; a < AA::kSize; ++a)
          sum += exp(crf[k].context_weights[i][a] / column_weights[i]);
        min = MIN(min, sum);
        max = MAX(max, sum);
        if (sum < DBL_MAX) {
          avg += sum;
          n++;
        }
      }
    }
    avg /= n;
    double loglike = s.loglike / tset_.size();
    double p = prior(crf) / crf.nweights();

    printf("%3zu: %5.2f %5.2f %5.2f %9.4f %9.4f\n", i, min, avg, max, p, i > 0 ? loglike : 0.0);
    sgd(s);
  }
}

}  // namespace cs
