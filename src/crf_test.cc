#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "crf-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "training_sequence.h"
#include "library_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"
#include "context_library-inl.h"

namespace cs {

const double kDelta = 0.01;
const double kDeltaTiny = 1e-6;
const std::string test_dir = PathCat(getenv("CS_DATA") ? getenv("CS_DATA") : "../data", "test");

class CrfTest : public testing::Test {
 protected:
       

  virtual void SetUp() {
    /*
    FILE* fin = fopen(PathCat(test_dir, "crf_state.cst"), "r");
    s1_.Read(fin);
    rewind(fin);
    s2_.Read(fin);
    rewind(fin);
    s3_.Read(fin);
    fclose(fin);
    */
  }

  CrfState<AA> s1_;
  CrfState<AA> s2_;
  CrfState<AA> s3_;
};

TEST_F(CrfTest, DISABLED_SimpleConstruction) {
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

TEST_F(CrfTest, DISABLED_ConstructionFromSerializedCRF) {
  FILE* fin = fopen(PathCat(test_dir, "test_crf.crf"), "r");
  Crf<AA> crf(fin);
  fclose(fin);

  EXPECT_NEAR(-1.696, crf[0].context_weights[0][AA::kCharToInt['A']], kDelta);
  EXPECT_NEAR( 0.135, crf[0].context_weights[0][AA::kCharToInt['R']], kDelta);
  EXPECT_NEAR( 0.000, crf[0].context_weights[0][AA::kCharToInt['X']], kDelta);
}


class CrfTestInit : public testing::Test {
  public:

    void IsEqual(const ContextLibrary<AA>& lib, const Crf<AA>& crf, 
                 double wcenter, double wdecay) {

      const size_t wlen = lib.wlen();
      const size_t c = 0.5 * (wlen - 1);
      double weights[wlen];
      for (size_t i = 1; i <= c; ++i) {
        double w = wcenter * pow(wdecay, i);
        weights[c - i] = w;
        weights[c + i] = w;
      }
      weights[c] = wcenter;

      EXPECT_EQ(lib.size(), crf.size());
      for (size_t k = 0; k < lib.size(); ++k) {
        EXPECT_NEAR(lib[k].prior, crf[k].bias_weight, kDelta);
        for (size_t i = 0; i < wlen; ++i)
          for (size_t a = 0; a < AA::kSize; ++a) 
            EXPECT_NEAR(weights[i] * lib[k].probs[i][a], crf[k].context_weights[i][a], kDelta);
        for (size_t a = 0; a < AA::kSize; ++a) 
          EXPECT_NEAR(lib[k].pc[a], crf[k].pc[a], kDelta);
      }
    }
};


TEST_F(CrfTestInit, DISABLED_SamplingCrfInit) {
  const Sequence<AA> seq("AHWTEXKLFADDF");
  const size_t nstates = 10;
  const size_t wlen    = seq.length();
  TrainingSequence<AA> tseq(seq, ProfileColumn<AA>(1.0));
  vector<TrainingSequence<AA> > tset(100, tseq);

  for (double tau = 0.0; tau <= 1.0; ++tau) {
    BlosumMatrix sm;
    ConstantAdmix admix(tau);
    MatrixPseudocounts<AA> pseudocounts(sm);
    SamplingCrfInit<AA, TrainingSequence<AA> > init(tset, pseudocounts, admix, sm, 0, 1.0, 1.0);
    Crf<AA> crf(nstates, wlen, init);

    const double bias_weight = log(1.0 / nstates);
    const double pc_weight = 0.0;
    const double pc = 1.0 / AA::kSize;
    for (size_t k = 0; k < crf.size(); ++k) {
      EXPECT_NEAR(bias_weight, crf[k].bias_weight, kDeltaTiny);
      for (size_t i = 0; i < wlen; ++i) {
        if (tau == 0.0) {
          double context_weight = log(DBL_MIN);
          for (size_t a = 0; a < AA::kSize; ++a)
            EXPECT_NEAR(a == seq[i] ? 0.0 : context_weight, crf[k].context_weights[i][a], kDeltaTiny);
        } else {
          ProfileColumn<AA> col_pc;
          for (size_t a = 0; a < AA::kSize; ++a) col_pc[a] = sm.r(a, seq[i]);
          Normalize(col_pc, 1.0);
          for (size_t a = 0; a < AA::kSize; ++a)
            EXPECT_NEAR(log(MAX(col_pc[a], DBL_MIN)), crf[k].context_weights[i][a], kDeltaTiny);
        }
      }
      for (size_t a = 0; a < AA::kSize; ++a) {
        EXPECT_NEAR(pc_weight, crf[k].pc_weights[a], kDeltaTiny);
        EXPECT_NEAR(pc, crf[k].pc[a], kDeltaTiny);
      }
    }
  }
}


TEST_F(CrfTestInit, DISABLED_GaussianCrfInit) {
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

TEST_F(CrfTestInit, DISABLED_LibraryCrfInit) {
  const size_t states  = 50;
  const size_t wlen    = 13;
  const double wcenter = 1.6;
  const double wdecay  = 0.85;
  const size_t c       = 0.5 * (wlen - 1);

  BlosumMatrix m;
  GaussianLibraryInit<AA> init(1.0, m);
  ContextLibrary<AA> lib(states, wlen, init);
  for (size_t k = 0; k < lib.size(); ++k) {
    lib[k].prior = MAX(lib[k].prior, DBL_MIN);
    lib[k].pc = ProfileColumn<AA>(lib[k].probs[c]);
    for (size_t a = 0; a < AA::kSize; ++a)
      lib[k].pc[a] = MAX(lib[k].pc[a], DBL_MIN);
  }
  TransformToLog(lib);
  Crf<AA> crf(lib, wcenter, wdecay);
  EXPECT_EQ(lib.size(), crf.size());
  IsEqual(lib, crf, wcenter, wdecay);

  LibraryBasedCrfInit<AA> lib_init(lib, wcenter, wdecay);
  crf = Crf<AA>(lib.size(), wlen, lib_init);
  IsEqual(lib, crf, wcenter, wdecay);

  FILE* fout = fopen(PathCat(test_dir, "equality.lib"), "w");
  lib.Write(fout);
  fclose(fout);
  fout = fopen(PathCat(test_dir, "equality.crf"), "w");
  crf.Write(fout);
  fclose(fout);
  FILE* fin = fopen(PathCat(test_dir, "equality.lib"), "r");
  ContextLibrary<AA> lib_s(fin);
  TransformToLog(lib_s);
  fclose(fin);
  fin = fopen(PathCat(test_dir, "equality.crf"), "r");
  Crf<AA> crf_s(fin);
  fclose(fin);
  IsEqual(lib_s, crf_s, wcenter, wdecay);
  remove(PathCat(test_dir, "equality.crf"));
  remove(PathCat(test_dir, "equality.lib"));


  Ran ran(0);
  ConstantAdmix admix(1.0);
  LibraryPseudocounts<AA> lib_pc(lib, wcenter, wdecay);
  CrfPseudocounts<AA> crf_pc(crf);
  for (size_t s = 0; s < 10; ++s) {
    Sequence<AA> seq(50);
    for (size_t i = 0; i < seq.length(); ++i) 
      seq[i] = static_cast<size_t>(ran(AA::kSize));
    CountProfile<AA> lib_profile(seq.length());
    lib_profile.counts = lib_pc.AddTo(seq, admix);
    CountProfile<AA> crf_profile(seq.length());
    crf_profile.counts = crf_pc.AddTo(seq, admix);
    for (size_t i = 0; i < seq.length(); ++i)
      for (size_t a = 0; a < AA::kSize; ++a) {
        EXPECT_TRUE(lib_profile.counts[i][a] > 0);
        EXPECT_NEAR(lib_profile.counts[i][a], crf_profile.counts[i][a], kDelta);
      }
  }
}

TEST_F(CrfTestInit, LibToCrf) {
  const double wcenter = 1.6;
  const double wdecay  = 0.85;
  std::string data_dir = getenv("CS_DATA") ? getenv("CS_DATA") : "../data";

  FILE* fin = fopen(PathCat(data_dir, "K4000.lib"), "r");
  ContextLibrary<AA> lib(fin);
  fclose(fin);
  Crf<AA> crf(lib, wcenter, wdecay);
  TransformToLog(lib);
  IsEqual(lib, crf, wcenter, wdecay);
  FILE* fout = fopen(PathCat(data_dir, "K4000.crf"), "w");
  crf.Write(fout);
  fclose(fout);
}


}  // namespace cs
