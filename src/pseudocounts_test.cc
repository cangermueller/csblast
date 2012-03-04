#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "po_hmm-inl.h"
#include "pseudocounts-inl.h"
#include "crf-inl.h"
#include "context_library-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "crf_pseudocounts-inl.h"
#include "library_pseudocounts-inl.h"


namespace cs {

using std::string;
using std::vector;

const double kDelta = 1e-5;
const string test_dir = PathCat(getenv("CS_DATA") ? getenv("CS_DATA") : "../data", "test");

class PseudocountsTest : public testing::Test {

  protected:
    typedef vector<Pseudocounts<AA>* > PcEngines;


    PseudocountsTest() : ran(Ran(0)) {}

    void SetUp() {
      FILE* fin = fopen(PathCat(test_dir, "K4000.lib"), "r");
      lib.reset(new ContextLibrary<AA>(fin));
      TransformToLog(*lib);
      fclose(fin);
      crf.reset(new Crf<AA>(*lib, 1.6, 0.85));

      pc_engines.push_back(new MatrixPseudocounts<AA>(sm));
      pc_engines.push_back(new LibraryPseudocounts<AA>(*lib, 1.6, 0.85));
      pc_engines.push_back(new CrfPseudocounts<AA>(*crf));
    }

    void TearDown() {
      for (size_t i = 0; i < pc_engines.size(); ++i) 
        delete pc_engines[i];
    }


    Sequence<AA> GetRndSeq(size_t len = kLen) {
      Sequence<AA> seq(len);
      for (size_t i = 0; i < len; ++i)
        seq[i] = ran(AA::kSize);
      return seq;
    }

    bool IsEqual(Profile<AA>& p, Profile<AA>& q, double delta = 0.0) {
      if (p.length() != q.length()) return false;
      for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < AA::kSize; ++a)
          if (fabs(p[i][a] - q[i][a]) > delta) return false;
      }
      return true;
    }

    Profile<AA> ToProfile(POHmm<AA>& hmm) {
      Profile<AA> p(hmm.size());
      for (size_t i = 1; i <= hmm.size(); ++i) {
        for (size_t a = 0; a < AA::kSize; ++a)
          p[i - 1][a] = hmm.g[i].probs[a];
      }
      return p;
    }

          

 protected:
    Ran ran;
    BlosumMatrix sm;
    scoped_ptr<ContextLibrary<AA> > lib;
    scoped_ptr<Crf<AA> > crf;
    PcEngines pc_engines;

    static const size_t kLen    = 50;
    static const size_t kRounds = 2;
};

TEST_F(PseudocountsTest, ConstantAdmix) {
  ConstantAdmix admix(0.0);
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    for (PcEngines::iterator it = pc_engines.begin(); it != pc_engines.end(); ++it) {
      Profile<AA> p = (*it)->AddTo(seq, admix);
      ASSERT_TRUE(IsEqual(p, cp.counts));
      p = (*it)->AddTo(cp, admix);
      ASSERT_TRUE(IsEqual(p, cp.counts));
    }
  }
}

TEST_F(PseudocountsTest, AdjustNeffSeq) {
  CSBlastAdmix admix(0.5, 0.5);
  const double kDeltaNeff = 0.01;
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    const double neff[] = {1.0, 4.0, 8.0};
    for (PcEngines::iterator it = pc_engines.begin(); it != pc_engines.end(); ++it) {
      for (size_t i = 0; i < 3; ++i) {
        // printf("%.3f   %.3f\n", 
        //    Neff((*it)->AddTo(seq, neff[i], admix, kDeltaNeff)), Neff((*it)->AddTo(cp, admix, neff[i], kDeltaNeff)));        
        ASSERT_NEAR(Neff((*it)->AddTo(seq, admix, neff[i], kDeltaNeff)), neff[i], kDeltaNeff);
        ASSERT_NEAR(Neff((*it)->AddTo(cp, admix, neff[i], kDeltaNeff)), neff[i], kDeltaNeff);
      }
      Profile<AA> p, q;

      // seq: a -> 0.0
      q = (*it)->AddTo(seq, admix, 1.0, 0.0);
      ASSERT_TRUE(IsEqual(cp.counts, q));

      // seq: a -> 1.0
      admix.pca = 1.0;
      p = (*it)->AddTo(seq, admix);
      q = (*it)->AddTo(seq, admix, 30, 0.0);
      ASSERT_TRUE(IsEqual(p, q));

      // cp: a -> 0.0
      q = (*it)->AddTo(cp, admix, 1.0, 0.0);
      ASSERT_TRUE(IsEqual(cp.counts, q));

      // cp: a -> 1.0
      admix.pca = 1.0;
      p = (*it)->AddTo(seq, admix);
      q = (*it)->AddTo(cp, admix, 30, 0.0);
      ASSERT_TRUE(IsEqual(p, q));
    }
  }
}

TEST_F(PseudocountsTest, AdjustNeffProf) {
  const double kDeltaNeff = 0.01;
  const double kNeffProf = 3.0;
  CSBlastAdmix admix(0.5, 12.0);

  for (size_t r = 0; r < kRounds; ++r) {
    // Prepare count profile cp with Neff = kNeffProf
    Sequence<AA> cp_seq = GetRndSeq();
    Profile<AA> cp_p = pc_engines[1]->AddTo(cp_seq, admix, kNeffProf, kDeltaNeff);
    ASSERT_NEAR(Neff(cp_p), kNeffProf, kDeltaNeff);
    CountProfile<AA> cp(cp_p);
    for (size_t i = 0; i < cp.length(); ++i) 
      cp.neff[i] = 2 * kNeffProf - ran(2 * kNeffProf - 1);
    Normalize(cp.counts, cp.neff);
    Profile<AA> pp(cp);

    const double neff[] = {4.0, 6.0, 8.0};
    for (PcEngines::iterator it = pc_engines.begin(); it != pc_engines.end(); ++it) {
      for (size_t i = 0; i < 3; ++i) {
        // printf("%.3f   %.3f\n", 
        //    Neff((*it)->AddTo(seq, neff[i], admix, kDeltaNeff)), Neff((*it)->AddTo(cp, admix, neff[i], kDeltaNeff)));        
        ASSERT_NEAR(Neff((*it)->AddTo(cp, admix, neff[i], kDeltaNeff)), neff[i], kDeltaNeff);
      }
      Profile<AA> p, q;

      // seq: a -> 0.0
      q = (*it)->AddTo(cp, admix, 1.5, 0.0);
      ASSERT_NEAR(Neff(cp_p), Neff(q), kDeltaNeff);
      ASSERT_TRUE(IsEqual(cp_p, q, 1e-5));

      // seq: a -> 1.0
      admix.pca = 1.0;
      p = (*it)->AddTo(cp, admix);
      q = (*it)->AddTo(cp, admix, 30, 0.0);
      ASSERT_TRUE(IsEqual(p, q));
    }
  }
}

TEST_F(PseudocountsTest, CrfLibEquality) {
  const double kNeffDelta = 0.001;
  const double tau[]      = {0.1, 0.5, 0.9};
  const double neff[]     = {3.0, 5.0, 8.0};
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    Profile<AA> p, q;
    for (size_t i = 0; i < 3; ++i) {
      ConstantAdmix admix(tau[i]);
      p = pc_engines[1]->AddTo(seq, admix);
      q = pc_engines[2]->AddTo(seq, admix);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
      p = pc_engines[1]->AddTo(cp, admix);
      q = pc_engines[2]->AddTo(cp, admix);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
    for (size_t i = 0; i < 3; ++i) {
      CSBlastAdmix admix(1.0, 12.0);
      p = pc_engines[1]->AddTo(seq, admix, neff[i], kNeffDelta);
      q = pc_engines[2]->AddTo(seq, admix, neff[i], kNeffDelta);
      ASSERT_NEAR(Neff(p), neff[i], kNeffDelta);
      ASSERT_NEAR(Neff(q), neff[i], kNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
      p = pc_engines[1]->AddTo(cp, admix, neff[i], kNeffDelta);
      q = pc_engines[2]->AddTo(cp, admix, neff[i], kNeffDelta);
      ASSERT_NEAR(Neff(p), neff[i], kNeffDelta);
      ASSERT_NEAR(Neff(q), neff[i], kNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
  }
}

TEST_F(PseudocountsTest, POHmm) {
  const double kNeffDelta = 0.001;
  const double tau[]      = {0.1, 0.5, 0.9};
  const double neff[]     = {3.0, 5.0, 8.0};
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    POHmm<AA> hmm(seq);
    Profile<AA> p, q;
    for (size_t i = 0; i < 3; ++i) {
      ConstantAdmix admix(tau[i]);
      p = pc_engines[1]->AddTo(seq, admix);
      pc_engines[1]->AddTo(&hmm, admix);
      q = ToProfile(hmm);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
    for (size_t i = 0; i < 3; ++i) {
      CSBlastAdmix admix(1.0, 12.0);
      p = pc_engines[1]->AddTo(seq, admix, neff[i], kNeffDelta);
      pc_engines[1]->AddTo(&hmm, admix, neff[i], kNeffDelta);
      q = ToProfile(hmm);
      ASSERT_NEAR(Neff(p), neff[i], kNeffDelta);
      ASSERT_NEAR(Neff(q), neff[i], kNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
  }
}



};  // namespace cs
