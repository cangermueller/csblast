#include <gtest/gtest.h>

#include <cstdio>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "mult_emission-inl.h"
#include "forward_backward_algorithm-inl.h"
#include "hmm-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

const float kFloatDelta = 0.0001;

class ForwardBackwardAlgorithmTest : public testing::Test {
 protected:
  ForwardBackwardAlgorithmTest()
      : hmm_(13, 1) {}

  virtual void SetUp() {
    Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
    Profile<AminoAcid> profile(seq.length());

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(&m);
    mpc.add_to_sequence(seq, ConstantAdmixture(0.1f), &profile);

    for (int i = 0; i < seq.length(); ++i) {
      Profile<AminoAcid> p(profile, i, 1);
      hmm_.AddState(p);
    }
    hmm_.InitTransitions(HomogeneousTransitionInitializerHMM<AminoAcid>());
    hmm_.TransformStatesToLogSpace();
  }

  HMM<AminoAcid> hmm_;
};

TEST_F(ForwardBackwardAlgorithmTest, ZincFingerMotif) {
  Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
  MultEmission<AminoAcid> emission(1);
  ForwardBackwardMatrices m(seq.length(), hmm_.num_states());
  ForwardBackwardAlgorithm(hmm_, seq, emission, &m);

  EXPECT_NEAR(0.9566, m.f[0][0] * m.b[0][0], kFloatDelta);
  EXPECT_NEAR(0.4920, m.f[1][1] * m.b[1][1], kFloatDelta);
  EXPECT_NEAR(0.9326, m.f[2][2] * m.b[2][2], kFloatDelta);
}

TEST_F(ForwardBackwardAlgorithmTest, 1Q7L) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
  fclose(fin);
  ali.AssignMatchColumnsByGapRule();
  CountProfile<AminoAcid> profile(ali, true);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.add_to_profile(ConstantAdmixture(0.01f), &profile);

  HMM<AminoAcid> hmm(profile.length(), 1);
  for (int i = 0; i < profile.length(); ++i) {
    Profile<AminoAcid> p(profile, i, 1);
    hmm.AddState(p);
  }
  hmm.InitTransitions(HomogeneousTransitionInitializerHMM<AminoAcid>());
  hmm.TransformStatesToLogSpace();

  MultEmission<AminoAcid> emission(1);
  ForwardBackwardMatrices mat(profile.length(), hmm.num_states());
  ForwardBackwardAlgorithm(hmm, profile, emission, &mat);
}

}  // namespace cs
