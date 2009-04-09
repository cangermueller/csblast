#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
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
      : hmm_(23, 1) { }

  virtual void SetUp() {
    // Build generic zinc finger HMM
    Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
    Profile<AminoAcid> profile(seq.length());

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(&m);
    mpc.add_to_sequence(seq, ConstantAdmixture(1.0f), &profile);

    for (int i = 0; i < seq.length(); ++i) {
      Profile<AminoAcid> p(profile, i, 1);
      hmm_.add_profile(p);
    }
    hmm_.init_transitions(HomogeneousTransitionInitializer<AminoAcid>());
    hmm_.transform_states_to_logspace();

    // Read zinc finger sequences
    std::ifstream seq_in("../data/zinc_finger_proteins.fas");
    Sequence<AminoAcid>::readall(seq_in, seqs_);
    seq_in.close();

    // Read zinc finger alignments and construct count profiles
    std::ifstream ali_in("../data/zinc_finger_alignments.fas");
    std::vector< shared_ptr< Alignment<AminoAcid> > > alis;
    Alignment<AminoAcid>::readall(ali_in, Alignment<AminoAcid>::FASTA, alis);
    ali_in.close();

    // Convert alignments to counts and add pseudocounts
    typedef std::vector< shared_ptr< Alignment<AminoAcid> > > ali_vector;
    typedef ali_vector::iterator alignment_iterator;
    for (alignment_iterator ai = alis.begin(); ai != alis.end(); ++ai) {
      shared_ptr< CountProfile<AminoAcid> > p_ptr(
          new CountProfile<AminoAcid>(**ai, true));
      mpc.add_to_profile(ConstantAdmixture(0.01f), p_ptr.get());
      p_ptr->convert_to_counts();
      counts_.push_back(p_ptr);
    }
  }

  HMM<AminoAcid> hmm_;
  std::vector< shared_ptr< Sequence<AminoAcid> > > seqs_;
  std::vector< shared_ptr< CountProfile<AminoAcid> > > counts_;
};

TEST_F(BaumWelchTrainingTest, ZincFingerSeqsTraining) {
  BaumWelchParams p;
  p.num_blocks         = 1;
  p.epsilon_null       = 1.0;
  p.beta               = 0.0;
  p.epsilon_batch      = 0.0f;
  p.min_scans          = 3;
  p.weight_center      = 1.0f;
  BaumWelchTraining<AminoAcid, Sequence> bwt(p, seqs_, hmm_, std::cout);
  bwt.run();

  hmm_.transform_states_to_linspace();
  EXPECT_NEAR(0.8977, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.8977, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisTraining) {
  BaumWelchParams p;
  p.transition_pseudocounts  = 0.8f;
  p.log_likelihood_change    = 0.02f;
  p.max_connectivity         = 5;
  p.num_blocks               = 1;
  p.epsilon_null             = 1.0f;
  p.beta                     = 0.0f;
  p.min_scans                = 3;
  p.weight_center            = 1.0f;
  BaumWelchTraining<AminoAcid, CountProfile> bwt(p, counts_, hmm_, std::cout);
  bwt.run();

  hmm_.transform_states_to_linspace();
  EXPECT_NEAR(0.9948, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.9948, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisOnlineTraining) {
  BaumWelchParams p;
  p.transition_pseudocounts  = 0.8f;
  p.log_likelihood_change    = 0.02f;
  p.max_connectivity         = 5;
  p.num_blocks               = 2;
  p.epsilon_null             = 0.05f;
  p.beta                     = 0.5f;
  p.min_scans                = 3;
  p.weight_center            = 1.0f;
  BaumWelchTraining<AminoAcid, CountProfile> bwt(p, counts_, hmm_, std::cout);
  bwt.run();

  hmm_.transform_states_to_linspace();
  EXPECT_NEAR(0.9948, hmm_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.9948, hmm_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

}  // namespace cs
