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

      pc_admix.push_back(new ConstantAdmix(0.5));
      pc_admix.push_back(new CSBlastAdmix(0.5, 12.0));
      pc_admix.push_back(new HHsearchAdmix(0.5, 1000));
    }

    void TearDown() {
      for (size_t i = 0; i < pc_engines.size(); ++i) 
        delete pc_engines[i];
      for (size_t i = 0; i < pc_admix.size(); ++i) 
        delete pc_admix[i];
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
    vector<Admix*> pc_admix;

    static const size_t kLen    = 50;
    static const size_t kRounds = 2;
};

TEST_F(PseudocountsTest, ConstantAdmix) {
  ConstantAdmix admix(0.0);
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    for (PcEngines::iterator it = pc_engines.begin(); it != pc_engines.end(); ++it) {
      (*it)->SetTargetNeff(0.0);
      Profile<AA> p = (*it)->AddTo(seq, admix);
      ASSERT_TRUE(IsEqual(p, cp.counts));
      p = (*it)->AddTo(cp, admix);
      ASSERT_TRUE(IsEqual(p, cp.counts));
    }
  }
}

TEST_F(PseudocountsTest, Sequence) {
  const double kTargetNeffDelta = 0.01;
  const double kTargetNeff[] = {1.0, 4.0, 8.0};
  
  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    for (PcEngines::iterator it_pc = pc_engines.begin(); it_pc != pc_engines.end(); ++it_pc) {
      Pseudocounts<AA>& pc = **it_pc;

      for (vector<Admix*>::iterator it_admix = pc_admix.begin(); it_admix != pc_admix.end(); ++it_admix) {
        Admix& admix = **it_admix;

        pc.SetTargetNeffDelta(kTargetNeffDelta);
        for (size_t i = 0; i < 3; ++i) {
          pc.SetTargetNeff(kTargetNeff[i]);
          ASSERT_NEAR(Neff(pc.AddTo(seq, admix)), kTargetNeff[i], kTargetNeffDelta);
          ASSERT_NEAR(Neff(pc.AddTo(cp, admix)), kTargetNeff[i], kTargetNeffDelta);
        }

        Profile<AA> p, q;
        pc.SetTargetNeffDelta(0.0);

        // seq: a -> 0.0
        pc.SetTargetNeff(1.0);
        q = pc.AddTo(seq, admix);
        ASSERT_TRUE(IsEqual(cp.counts, q));

        // cp: a -> 0.0
        pc.SetTargetNeff(1.0);
        q = pc.AddTo(cp, admix);
        ASSERT_TRUE(IsEqual(cp.counts, q));

        // seq: a -> 1.0
        admix.SetTargetNeffParam(1.0);
        pc.SetTargetNeff(0.0);
        p = pc.AddTo(seq, admix);
        pc.SetTargetNeff(30.0);
        q = pc.AddTo(seq, admix);
        ASSERT_TRUE(IsEqual(p, q));

        // cp: a -> 1.0
        admix.SetTargetNeffParam(1.0);
        pc.SetTargetNeff(0.0);
        p = pc.AddTo(cp, admix);
        pc.SetTargetNeff(30.0);
        q = pc.AddTo(cp, admix);
        ASSERT_TRUE(IsEqual(p, q));
      }
    }
  }
}

TEST_F(PseudocountsTest, CountProfile) {
  const double kTargetNeffDelta = 0.01;
  const double kTargetNeff[] = {4.0, 6.0, 8.0};
  const double kTargetNeffProf = 3.0;
  
  for (size_t r = 0; r < kRounds; ++r) {
    // Create random count profile cp with kTargetNeffProf
    Sequence<AA> cp_seq = GetRndSeq();
    pc_engines[1]->SetTargetNeff(kTargetNeffProf);
    pc_engines[1]->SetTargetNeffDelta(kTargetNeffDelta);
    Profile<AA> cp_p = pc_engines[1]->AddTo(cp_seq, *pc_admix[1]);
    ASSERT_NEAR(Neff(cp_p), kTargetNeffProf, kTargetNeffDelta);
    CountProfile<AA> cp(cp_p);
    for (size_t i = 0; i < cp.length(); ++i) 
      cp.neff[i] = 2 * kTargetNeffProf - ran(2 * kTargetNeffProf - 1);
    Normalize(cp.counts, cp.neff);

    for (PcEngines::iterator it_pc = pc_engines.begin(); it_pc != pc_engines.end(); ++it_pc) {
      Pseudocounts<AA>& pc = **it_pc;

      for (vector<Admix*>::iterator it_admix = pc_admix.begin(); it_admix != pc_admix.end(); ++it_admix) {
        Admix& admix = **it_admix;

        pc.SetTargetNeffDelta(kTargetNeffDelta);
        for (size_t i = 0; i < 3; ++i) {
          pc.SetTargetNeff(kTargetNeff[i]);
          ASSERT_NEAR(Neff(pc.AddTo(cp, admix)), kTargetNeff[i], kTargetNeffDelta);
        }

        Profile<AA> q;
        pc.SetTargetNeffDelta(0.0);

        // pca -> 0.0
        pc.SetTargetNeff(kTargetNeffProf - 1);
        q = pc.AddTo(cp, admix);
        ASSERT_NEAR(Neff(q), Neff(cp_p), kTargetNeffDelta);
        ASSERT_TRUE(IsEqual(cp_p, q, 1e-5));
      }
    }
  }
}

TEST_F(PseudocountsTest, POHmm) {
  const double kTargetNeffDelta = 0.001;
  const double kTargetNeff[]    = {3.0, 5.0, 8.0};
  const double tau[]            = {0.1, 0.5, 0.9};
  pc_engines[1]->SetTargetNeffDelta(kTargetNeffDelta);

  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    POHmm<AA> hmm(seq);
    Profile<AA> p, q;
    pc_engines[1]->SetTargetNeff(0.0);
    for (size_t i = 0; i < 3; ++i) {
      ConstantAdmix admix(tau[i]);
      p = pc_engines[1]->AddTo(seq, admix);
      pc_engines[1]->AddTo(&hmm, admix);
      q = ToProfile(hmm);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
    for (size_t i = 0; i < 3; ++i) {
      pc_engines[1]->SetTargetNeff(kTargetNeff[i]);
      CSBlastAdmix admix(1.0, 12.0);

      p = pc_engines[1]->AddTo(seq, admix);
      pc_engines[1]->AddTo(&hmm, admix);
      q = ToProfile(hmm);
      ASSERT_NEAR(Neff(p), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_NEAR(Neff(q), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
  }
}

TEST_F(PseudocountsTest, CrfLibEquality) {
  const double kTargetNeffDelta = 0.001;
  const double kTargetNeff[]    = {3.0, 5.0, 8.0};
  const double tau[]            = {0.1, 0.5, 0.9};

  for (size_t r = 0; r < kRounds; ++r) {
    Sequence<AA> seq = GetRndSeq();
    CountProfile<AA> cp(seq);
    Profile<AA> p, q;

    for (size_t i = 1; i <= 2; ++i) {
      pc_engines[i]->SetTargetNeffDelta(kTargetNeffDelta);
      pc_engines[i]->SetTargetNeff(0.0);
    }

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
      pc_engines[1]->SetTargetNeff(kTargetNeff[i]);
      pc_engines[2]->SetTargetNeff(kTargetNeff[i]);
      CSBlastAdmix admix(1.0, 12.0);

      p = pc_engines[1]->AddTo(seq, admix);
      q = pc_engines[2]->AddTo(seq, admix);
      ASSERT_NEAR(Neff(p), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_NEAR(Neff(q), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));

      p = pc_engines[1]->AddTo(cp, admix);
      q = pc_engines[2]->AddTo(cp, admix);
      ASSERT_NEAR(Neff(p), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_NEAR(Neff(q), kTargetNeff[i], kTargetNeffDelta);
      ASSERT_TRUE(IsEqual(p, q, 0.001));
    }
  }
}




};  // namespace cs
