#include <gtest/gtest.h>

#include <cstdio>

#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "clustering-inl.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

const float kFloatDelta = 0.01;

class ClusteringTest : public testing::Test {
 protected:
  static const int kWindowLength = 5;
  static const int kNumStates    = 10;

  ClusteringTest()
      : lib_(kNumStates, kWindowLength) {}

  virtual void SetUp() {
    // Build generic zinc finger HMM
    Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
    Profile<AminoAcid> profile(seq.length());

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(&m);
    mpc.AddPseudocountsToSequence(seq, ConstantAdmixture(1.0f), &profile);

    // Initialize profile library
    for (int i = 0; i < kNumStates; ++i) {
      Profile<AminoAcid> p(profile, i, kWindowLength);
      lib_.AddProfile(p);
    }
    lib_.TransformToLogSpace();

    // Read zinc finger alignments
    std::vector< shared_ptr< Alignment<AminoAcid> > > alis;
    FILE* ali_in = fopen("../data/zinc_finger_alignments.fas", "r");
    Alignment<AminoAcid>::ReadAll(ali_in, Alignment<AminoAcid>::FASTA, &alis);
    fclose(ali_in);

    // Convert alignments to counts and add pseudocounts
    typedef std::vector< shared_ptr< Alignment<AminoAcid> > > ali_vector;
    typedef ali_vector::iterator alignment_iterator;
    for (alignment_iterator ai = alis.begin(); ai != alis.end(); ++ai) {
      CountProfile<AminoAcid> p_full(**ai, true);
      for (int i = 0; i < p_full.num_cols() - kWindowLength + 1; ++i) {
        shared_ptr< CountProfile<AminoAcid> > p_ptr(
            new CountProfile<AminoAcid>(p_full, i, kWindowLength));
        mpc.AddPseudocountsToProfile(ConstantAdmixture(0.01f), p_ptr.get());
        counts_.push_back(p_ptr);
      }
    }
  }

  ProfileLibrary<AminoAcid> lib_;
  std::vector< shared_ptr< CountProfile<AminoAcid> > > counts_;
};

TEST_F(ClusteringTest, ZincFingerAlisClustering) {
  ClusteringOptions p;
  p.num_blocks   = 1;
  p.epsilon_null = 1.0f;
  p.beta         = 0.0f;
  p.min_scans    = 5;
  Clustering<AminoAcid, CountProfile> em_clust(p, counts_, lib_, stdout);
  em_clust.Run();

  lib_.TransformToLinSpace();
  EXPECT_NEAR(0.9948, lib_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.4975, lib_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(ClusteringTest, ZincFingerAlisOnlineClustering) {
  ClusteringOptions p;
  p.num_blocks   = 2;
  p.epsilon_null = 0.05f;
  p.beta         = 0.5f;
  p.min_scans    = 5;
  Clustering<AminoAcid, CountProfile> em_clust(p, counts_, lib_, stdout);
  em_clust.Run();

  lib_.TransformToLinSpace();
  EXPECT_NEAR(0.9948, lib_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.4975, lib_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

}  // namespace cs
