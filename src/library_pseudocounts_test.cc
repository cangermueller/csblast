#include <gtest/gtest.h>

#include <iostream>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "emitter-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile_library-inl.h"
#include "library_pseudocounts-inl.h"

namespace cs {

const float kFloatDelta = 0.01f;

TEST(LibraryPseudocountsTest, AddToDummySequence) {
  Sequence<AminoAcid> seq("header", "ARNDCQEGHILKMFPSTWYV");
  Profile<AminoAcid> profile(seq.length());

  ASSERT_EQ(AminoAcid::instance().size(), seq.length());
  ASSERT_EQ(AminoAcid::instance().ctoi('R'), seq[1]);
  ASSERT_EQ(seq.length(), profile.num_cols());

  FILE* fin = fopen("../data/scop20_1.73_opt_N100000_W13.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  ASSERT_EQ(50, lib.num_profiles());

  EmissionOptions opts;
  LibraryPseudocounts<AminoAcid> pc(&lib, opts);
  pc.add_to_sequence(seq, DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

  EXPECT_NEAR(0.0736f, profile[0][AminoAcid::instance().ctoi('V')], kFloatDelta);
}

TEST(LibraryPseudocountsTest, AddToSmallProfile) {
  FILE* ali_in = fopen("../data/zinc_finger_alignments.fas", "r");
  Alignment<AminoAcid> ali(ali_in, Alignment<AminoAcid>::FASTA);
  fclose(ali_in);
  CountProfile<AminoAcid> profile(ali, false);

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.add_to_profile(ConstantAdmixture(0.1f), &profile);

  FILE* fin = fopen("../data/scop20_1.73_opt_N100000_W13.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  ASSERT_EQ(50, lib.num_profiles());

  EmissionOptions opts;
  LibraryPseudocounts<AminoAcid> pc(&lib, opts);
  pc.add_to_profile(DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

  EXPECT_NEAR(0.80f, profile[0][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.81f, profile[5][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST(LibraryPseudocountsTest, AddToZnFingerSequenceK50) {
  FILE* seq_in = fopen("../data/zinc_finger.seq", "r");
  Sequence<AminoAcid> seq(seq_in);
  fclose(seq_in);

  Profile<AminoAcid> profile(seq.length());

  FILE* fin = fopen("../data/scop20_1.73_opt_N100000_W13.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  EmissionOptions opts;
  LibraryPseudocounts<AminoAcid> pc(&lib, opts);
  for (int n = 0; n < 80; ++n)
    pc.add_to_sequence(seq, ConstantAdmixture(1.0f), &profile);

  EXPECT_NEAR(0.8259, profile[53][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.8214, profile[56][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST(LibraryPseudocountsTest, AddToZnFingerSequenceK4000) {
  FILE* seq_in = fopen("../data/zinc_finger.seq", "r");
  Sequence<AminoAcid> seq(seq_in);
  fclose(seq_in);
  Profile<AminoAcid> profile(seq.length());

  FILE* fin = fopen("../data/K4000.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  ASSERT_EQ(4000, lib.num_profiles());
  ASSERT_EQ(13, lib.num_cols());
  ASSERT_EQ(20, lib.alphabet_size());
  ASSERT_TRUE(lib.logspace());

  EmissionOptions opts;
  LibraryPseudocounts<AminoAcid> pc(&lib, opts);
  pc.add_to_sequence(seq, ConstantAdmixture(1.0f), &profile);

  EXPECT_NEAR(0.78696, profile[53][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.84759, profile[56][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST(LibraryPseudocountsTest, AddToZnFingerAlignmentK4000) {
  FILE* ali_in = fopen("../data/zinc_finger.fas", "r");
  Alignment<AminoAcid> ali(ali_in, Alignment<AminoAcid>::FASTA);
  fclose(ali_in);
  CountProfile<AminoAcid> profile(ali, false);

  FILE* fin = fopen("../data/K4000.lib", "r");
  ProfileLibrary<AminoAcid> lib(fin);
  fclose(fin);

  ASSERT_EQ(4000, lib.num_profiles());
  ASSERT_EQ(13, lib.num_cols());
  ASSERT_EQ(20, lib.alphabet_size());
  ASSERT_TRUE(lib.logspace());

  EmissionOptions opts;
  LibraryPseudocounts<AminoAcid> pc(&lib, opts);
  pc.add_to_profile(DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

  EXPECT_NEAR(0.7612, profile[53][AminoAcid::instance().ctoi('C')], kFloatDelta);
  EXPECT_NEAR(0.7925, profile[56][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

}  // namespace cs
