#include <gtest/gtest.h>

#include <cstdio>

#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "baum_welch_training-inl.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "hmm-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

const float kFloatDelta = 0.01;

class BaumWelchTrainingTest : public testing::Test {
 protected:
  BaumWelchTrainingTest()
      : hmm_(23, 1) {}

  virtual void SetUp() {
    // Build generic zinc finger HMM
    Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
    Profile<AminoAcid> profile(seq.length());

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(&m);
    mpc.add_to_sequence(seq, ConstantAdmixture(1.0f), &profile);

    for (int i = 0; i < seq.length(); ++i) {
      Profile<AminoAcid> p(profile, i, 1);
      hmm_.AddState(p);
    }
    hmm_.InitTransitions(HomogeneousTransitionInitializerHMM<AminoAcid>());
    hmm_.TransformStatesToLogSpace();

    // Read zinc finger sequences
    FILE* seq_in = fopen("../data/zinc_finger_proteins.fas", "r");
    Sequence<AminoAcid>::ReadAll(seq_in, &seqs_);
    fclose(seq_in);

    // Read zinc finger alignments and construct count profiles
    std::vector< shared_ptr< Alignment<AminoAcid> > > alis;
    FILE* ali_in = fopen("../data/zinc_finger_alignments.fas", "r");
    Alignment<AminoAcid>::ReadAll(ali_in, Alignment<AminoAcid>::FASTA, &alis);
    fclose(ali_in);

    // Convert alignments to counts and add pseudocounts
    typedef std::vector< shared_ptr< Alignment<AminoAcid> > > ali_vector;
    typedef ali_vector::iterator alignment_iterator;
    for (alignment_iterator ai = alis.begin(); ai != alis.end(); ++ai) {
      shared_ptr< CountProfile<AminoAcid> > p_ptr(
          new CountProfile<AminoAcid>(**ai, true));
      mpc.add_to_profile(ConstantAdmixture(0.01f), p_ptr.get());
      counts_.push_back(p_ptr);
    }
  }

  HMM<AminoAcid> hmm_;
  std::vector< shared_ptr< Sequence<AminoAcid> > > seqs_;
  std::vector< shared_ptr< CountProfile<AminoAcid> > > counts_;
};

TEST_F(BaumWelchTrainingTest, ZincFingerSeqsTraining) {
  BaumWelchOptions p;
  p.num_blocks         = 1;
  p.epsilon_null       = 1.0;
  p.beta               = 0.0;
  p.epsilon_batch      = 0.0f;
  p.min_scans          = 3;
  p.weight_center      = 1.0f;
  BaumWelchTraining<AminoAcid, Sequence> bwt(p, seqs_, hmm_, stdout);
  bwt.Run();

  hmm_.TransformStatesToLinSpace();
  EXPECT_NEAR(0.8977, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.8977, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisTraining) {
  BaumWelchOptions p;
  p.transition_pc  = 0.8f;
  p.log_likelihood_change    = 0.02f;
  p.max_connectivity         = 5;
  p.num_blocks               = 1;
  p.epsilon_null             = 1.0f;
  p.beta                     = 0.0f;
  p.min_scans                = 3;
  p.weight_center            = 1.0f;
  BaumWelchTraining<AminoAcid, CountProfile> bwt(p, counts_, hmm_, stdout);
  bwt.Run();

  hmm_.TransformStatesToLinSpace();
  EXPECT_NEAR(0.9948, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.9948, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisOnlineTraining) {
  BaumWelchOptions p;
  p.transition_pc  = 0.8f;
  p.log_likelihood_change    = 0.02f;
  p.max_connectivity         = 5;
  p.num_blocks               = 2;
  p.epsilon_null             = 0.05f;
  p.beta                     = 0.5f;
  p.min_scans                = 3;
  p.weight_center            = 1.0f;
  BaumWelchTraining<AminoAcid, CountProfile> bwt(p, counts_, hmm_, stdout);
  bwt.Run();

  hmm_.TransformStatesToLinSpace();
  EXPECT_NEAR(0.9948, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.9948, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

}  // namespace cs
