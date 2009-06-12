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

const double kDelta = 0.0001;

TEST(ForwardBackwardAlgorithmTest, ZincFingerMotif) {
  Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
  Profile<AminoAcid> profile(seq.length());

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.AddPseudocountsToSequence(seq, ConstantAdmixture(0.1f), &profile);

  HMM<AminoAcid> hmm(profile.length(), 1);
  for (int i = 0; i < seq.length(); ++i) {
    Profile<AminoAcid> p(profile, i, 1);
    hmm.AddState(p);
  }
  hmm.InitTransitions(HomogeneousTransitionInitializerHMM<AminoAcid>());
  hmm.TransformStatesToLogSpace();

  MultEmission<AminoAcid> emission(1);
  ForwardBackwardMatrices mat(seq.length(), hmm.num_states());
  ForwardBackwardAlgorithm(hmm, seq, emission, &mat);

  EXPECT_NEAR(0.9566, mat.f[0][0] * mat.b[0][0], kDelta);
  EXPECT_NEAR(0.4920, mat.f[1][1] * mat.b[1][1], kDelta);
  EXPECT_NEAR(0.9326, mat.f[2][2] * mat.b[2][2], kDelta);
}

TEST(ForwardBackwardAlgorithmTest, 1Q7L) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
  fclose(fin);
  ali.AssignMatchColumnsByGapRule();
  CountProfile<AminoAcid> profile(ali, true);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.AddPseudocountsToProfile(ConstantAdmixture(0.01f), &profile);

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

  EXPECT_NEAR(0.6116, mat.f[0][0] * mat.b[0][0], kDelta);
  EXPECT_NEAR(0.5003, mat.f[1][1] * mat.b[1][1], kDelta);
  EXPECT_NEAR(0.1685, mat.f[2][2] * mat.b[2][2], kDelta);
}

}  // namespace cs
