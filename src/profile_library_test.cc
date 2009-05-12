#include <gtest/gtest.h>

#include <cstdio>

#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "blosum_matrix.h"
#include "profile_library-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "nucleotide.h"
#include "profile-inl.h"
#include "shared_ptr.h"

namespace cs {

class ProfileLibraryTest : public testing::Test {
 protected:

  virtual void SetUp() {
    FILE* fin = fopen("../data/hmm_profile.prf", "r");
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

TEST_F(ProfileLibraryTest, SimpleConstruction) {
  ProfileLibrary<Nucleotide> lib(3, 5);

  EXPECT_FALSE(lib.full());
  EXPECT_EQ(3, lib.num_profiles());
  EXPECT_EQ(5, lib.num_cols());

  int index_p1 = lib.AddState(p1_);
  EXPECT_EQ(0, index_p1);
  int index_p2 = lib.AddState(p2_);
  EXPECT_EQ(1, index_p2);
  int index_p3 = lib.AddState(p3_);
  EXPECT_EQ(2, index_p3);

  EXPECT_FLOAT_EQ(1.0f, lib[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[0][1][0]);
  EXPECT_FLOAT_EQ(1.0f, lib[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[1][1][0]);
  EXPECT_FLOAT_EQ(1.0f, lib[2][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[2][1][0]);

  lib[2][0][0] = 0.0f;
  lib[2][0][1] = 1.0f;
  EXPECT_FLOAT_EQ(0.0f, lib[2][0][0]);
  EXPECT_FLOAT_EQ(1.0f, lib[2][0][1]);
}

TEST_F(ProfileLibraryTest, ConstructionFromSerializedLibrary) {
  FILE* fin = fopen("../data/scop20_1.73_opt_N100000_W13.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  EXPECT_EQ(50, lib.num_profiles());
  EXPECT_EQ(13, lib.num_cols());
  EXPECT_TRUE(lib.logspace());
  EXPECT_TRUE(lib[0].logspace());
}

TEST_F(ProfileLibraryTest, LinLogTransformation) {
  ProfileLibrary<Nucleotide> lib(3, 5);
  lib.AddState(p1_);
  lib.AddState(p2_);
  lib.AddState(p3_);

  lib.transform_to_logspace();

  EXPECT_FLOAT_EQ(0.0f, lib[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[2][0][0]);

  lib.transform_to_linspace();

  EXPECT_FLOAT_EQ(1.0f, lib[0][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[0][1][0]);
  EXPECT_FLOAT_EQ(1.0f, lib[1][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[1][1][0]);
  EXPECT_FLOAT_EQ(1.0f, lib[2][0][0]);
  EXPECT_FLOAT_EQ(0.0f, lib[2][1][0]);
}

TEST(ProfileLibraryTestInitialization, RandomSampleInitializer) {
  FILE* fin = fopen("../data/1Q7L.fas", "r");
  Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
  fclose(fin);
  ali.AssignMatchColumnsByGapRule();

  typedef shared_ptr< CountProfile<AminoAcid> > profile_ptr;
  CountProfile<AminoAcid> p_full(ali, true);
  profile_ptr p_window(new CountProfile<AminoAcid>(p_full, 0, 13));
  std::vector<profile_ptr> profiles(100, p_window);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> pc(&m);
  SamplingProfileInitializer<AminoAcid> profile_init(profiles, &pc, 0.2f);
  ProfileLibrary<AminoAcid> lib(10, 13, profile_init);

  EXPECT_EQ(10, lib.num_profiles());
  EXPECT_EQ(13, lib.num_cols());
  EXPECT_EQ(13, lib[0].num_cols());
}

}  // namespace cs
