#include <gtest/gtest.h>

#include <cstdio>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "blosum_matrix.h"
#include "hmm-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "nucleotide.h"
#include "profile-inl.h"
#include "shared_ptr.h"

namespace cs {

class HMMTest : public testing::Test {
 protected:

  virtual void SetUp() {
    FILE* fin = fopen("../data/hmm_profile.prf", "r");
    p1_.read(fin);
    rewind(fin);
    p2_.read(fin);
    rewind(fin);
    p3_.read(fin);
    fclose(fin);
  }

  Profile<Nucleotide> p1_;
  Profile<Nucleotide> p2_;
  Profile<Nucleotide> p3_;
};

TEST_F(HMMTest, SimpleConstruction) {
  HMM<Nucleotide> hmm(3, 5);

  EXPECT_FALSE(hmm.full());
  EXPECT_EQ(0, hmm.num_transitions());
  EXPECT_EQ(3, hmm.num_states());

  int index_p1 = hmm.add_profile(p1_);
  EXPECT_EQ(0, index_p1);
  int index_p2 = hmm.add_profile(p2_);
  EXPECT_EQ(1, index_p2);
  int index_p3 = hmm.add_profile(p3_);
  EXPECT_EQ(2, index_p3);

  EXPECT_FLOAT_EQ(1.0f, hmm[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[0][1][0]);
  EXPECT_FLOAT_EQ(1.0f, hmm[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[1][1][0]);
  EXPECT_FLOAT_EQ(1.0f, hmm[2][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[2][1][0]);

  hmm[2][0][0] = 0.0f;
  hmm[2][0][1] = 1.0f;
  EXPECT_FLOAT_EQ(0.0f, hmm[2][0][0]);
  EXPECT_FLOAT_EQ(1.0f, hmm[2][0][1]);

  EXPECT_FLOAT_EQ(0.0f, hmm(0,1));
  hmm(0,1) = 0.5f;
  EXPECT_EQ(1, hmm[0].num_out_transitions());
  EXPECT_EQ(0, hmm[0].num_in_transitions());
  EXPECT_EQ(1, hmm[1].num_in_transitions());
  EXPECT_EQ(0, hmm[1].num_out_transitions());
  EXPECT_FLOAT_EQ(0.5f, hmm(0,1));

  EXPECT_FLOAT_EQ(0.0f, hmm(0,2));
  hmm(0,2) = 0.5f;
  EXPECT_FLOAT_EQ(0.5f, hmm(0,2));

  hmm(1,2) = 1.0f;
  EXPECT_FLOAT_EQ(1.0f, hmm(1,2));

  EXPECT_EQ(3, hmm.num_transitions());
}

TEST_F(HMMTest, ConstructionFromSerializedHMM) {
  FILE* fin = fopen("../data/scop20_K100.hmm", "r");
  HMM<AminoAcid> hmm(fin);
  fclose(fin);

  EXPECT_EQ(869, hmm.num_transitions());
  EXPECT_EQ(100, hmm.num_states());
  EXPECT_EQ(13, hmm.num_cols());
  EXPECT_EQ(34, hmm.iterations());
  EXPECT_TRUE(hmm.states_logspace());
  EXPECT_FALSE(hmm.transitions_logspace());
  EXPECT_TRUE(hmm[0].logspace());
}

TEST_F(HMMTest, DISABLED_ConstructionFromSerializedHMM2) {
  std::ifstream fin("../data/scop20_K100.hmm");
  HMM<AminoAcid> hmm(fin);
  fin.close();

  EXPECT_EQ(869, hmm.num_transitions());
  EXPECT_EQ(100, hmm.num_states());
  EXPECT_EQ(13, hmm.num_cols());
  EXPECT_EQ(34, hmm.iterations());
  EXPECT_TRUE(hmm.states_logspace());
  EXPECT_FALSE(hmm.transitions_logspace());
}

TEST_F(HMMTest, NormalizeTransitions) {
  HMM<Nucleotide> hmm(3, 5);
  hmm.add_profile(p1_);
  hmm.add_profile(p2_);
  hmm.add_profile(p3_);
  hmm(0,1) = 0.3f;
  hmm(0,2) = 0.3f;
  hmm(1,2) = 0.9f;
  hmm(2,2) = 0.6f;

  normalize_transitions(hmm);

  EXPECT_FLOAT_EQ(0.5f, hmm(0,1));
  EXPECT_FLOAT_EQ(0.5f, hmm(0,2));
  EXPECT_FLOAT_EQ(1.0f, hmm(1,2));
  EXPECT_FLOAT_EQ(1.0f, hmm(2,2));
}

TEST_F(HMMTest, LinLogTransformation) {
  HMM<Nucleotide> hmm(3, 5);
  hmm.add_profile(p1_);
  hmm.add_profile(p2_);
  hmm.add_profile(p3_);
  hmm(0,1) = 0.5f;
  hmm(0,2) = 0.5f;
  hmm(1,2) = 1.0f;

  hmm.transform_transitions_to_logspace();
  hmm.transform_states_to_logspace();

  EXPECT_FLOAT_EQ(0.0f, hmm[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[2][0][0]);

  EXPECT_FLOAT_EQ(0.0f, hmm(1,2));

  hmm.transform_transitions_to_linspace();
  hmm.transform_states_to_linspace();

  EXPECT_FLOAT_EQ(1.0f, hmm[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[0][1][0]);
  EXPECT_FLOAT_EQ(1.0f, hmm[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[1][1][0]);
  EXPECT_FLOAT_EQ(1.0f, hmm[2][0][0]);
  EXPECT_FLOAT_EQ(0.0f, hmm[2][1][0]);

  EXPECT_FLOAT_EQ(1.0f, hmm(1,2));
}

TEST(HMMTestInitialization, RandomSampleInitializer) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
  fclose(fin);
  ali.assign_match_columns_by_gap_rule();

  typedef shared_ptr< CountProfile<AminoAcid> > profile_ptr;
  profile_ptr p(new CountProfile<AminoAcid>(ali, true));
  std::vector<profile_ptr> profiles(10, p);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> pc(&m);
  SamplingStateInitializer<AminoAcid> st_init(profiles, 0.2f, &pc, 0.2f);
  HomogeneousTransitionInitializer<AminoAcid> tr_init;
  HMM<AminoAcid> hmm(10, 5, st_init, tr_init);

  EXPECT_EQ(10, hmm.num_states());
  EXPECT_EQ(100, hmm.num_transitions());
  EXPECT_EQ(5, hmm.num_cols());
  EXPECT_EQ(5, hmm[0].num_cols());
}

}  // namespace cs
