#include <gtest/gtest.h>

#include <cstdio>

#include <vector>

#include "amino_acid.h"
#include "crf-inl.h"
#include "log.h"
#include "nucleotide.h"
#include "profile-inl.h"
#include "shared_ptr.h"

namespace cs {

class CRFTest : public testing::Test {
 protected:

  virtual void SetUp() {
    FILE* fin = fopen("../data/profile.prf", "r");
    p1_.Read(fin);
    rewind(fin);
    p2_.Read(fin);
    rewind(fin);
    p3_.Read(fin);
    fclose(fin);
  }

  Profile<Nucleotide> p1_;
  Profile<Nucleotide> p2_;
  Profile<Nucleotide> p3_;
};

TEST_F(CRFTest, SimpleConstruction) {
  CRF<Nucleotide> crf(3, 5);

  EXPECT_FALSE(crf.full());
  EXPECT_EQ(0, crf.num_transitions());
  EXPECT_EQ(3, crf.num_states());

  int index_p1 = crf.AddState(p1_);
  EXPECT_EQ(0, index_p1);
  int index_p2 = crf.AddState(p2_);
  EXPECT_EQ(1, index_p2);
  int index_p3 = crf.AddState(p3_);
  EXPECT_EQ(2, index_p3);

  EXPECT_FLOAT_EQ(0.0f, crf[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, crf[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, crf[2][0][0]);

  crf[2][0][0] = 0.0f;
  crf[2][0][1] = 1.0f;
  EXPECT_FLOAT_EQ(0.0f, crf[2][0][0]);
  EXPECT_FLOAT_EQ(1.0f, crf[2][0][1]);

  EXPECT_FLOAT_EQ(0.0f, crf(0,1));
  crf(0,1) = 0.5f;
  EXPECT_EQ(1, crf[0].num_out_transitions());
  EXPECT_EQ(0, crf[0].num_in_transitions());
  EXPECT_EQ(1, crf[1].num_in_transitions());
  EXPECT_EQ(0, crf[1].num_out_transitions());
  EXPECT_FLOAT_EQ(0.5f, crf(0,1));

  EXPECT_FLOAT_EQ(0.0f, crf(0,2));
  crf(0,2) = 0.5f;
  EXPECT_FLOAT_EQ(0.5f, crf(0,2));

  crf(1,2) = 1.0f;
  EXPECT_FLOAT_EQ(1.0f, crf(1,2));

  EXPECT_EQ(3, crf.num_transitions());
}

TEST_F(CRFTest, ConstructionFromSerializedCRF) {
  FILE* fin = fopen("../data/scop20_K100.crf", "r");
  CRF<AminoAcid> crf(fin);
  fclose(fin);

  EXPECT_EQ(864, crf.num_transitions());
  EXPECT_EQ(100, crf.num_states());
  EXPECT_EQ(11, crf.num_cols());
  EXPECT_EQ(34, crf.iterations());
  EXPECT_FALSE(crf.transitions_logspace());
}

TEST_F(CRFTest, NormalizeTransitions) {
  CRF<Nucleotide> crf(3, 5);
  crf.AddState(p1_);
  crf.AddState(p2_);
  crf.AddState(p3_);
  crf(0,1) = 0.3f;
  crf(0,2) = 0.3f;
  crf(1,2) = 0.9f;
  crf(2,2) = 0.6f;

  NormalizeTransitions(&crf);

  EXPECT_FLOAT_EQ(0.5f, crf(0,1));
  EXPECT_FLOAT_EQ(0.5f, crf(0,2));
  EXPECT_FLOAT_EQ(1.0f, crf(1,2));
  EXPECT_FLOAT_EQ(1.0f, crf(2,2));
}

}  // namespace cs
