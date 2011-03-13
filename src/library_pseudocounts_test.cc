#include <gtest/gtest.h>

#include "cs.h"
#include "context_library.h"
#include "library_pseudocounts-inl.h"
#include "sequence-inl.h"

namespace cs {

TEST(LibraryPseudocountsTest, K4000ZincFingerSequence) {
  FILE* seq_in = fopen("../data/zinc_finger.seq", "r");
  Sequence<AA> seq(seq_in);
  fclose(seq_in);

  FILE* fin = fopen("../data/nr30_neff2.5_1psi_W13_N1000000_K4000.lib", "r");
  ContextLibrary<AA> lib(fin);
  fclose(fin);
  TransformToLog(lib);

  LibraryPseudocounts<AA> pc(lib, 1.6, 0.85);
  Profile<AA> profile(pc.AddTo(seq, ConstantAdmix(1.0)));

  EXPECT_NEAR(0.85495, profile[53][AA::kCharToInt['C']], 0.0001);
  EXPECT_NEAR(0.85437, profile[56][AA::kCharToInt['C']], 0.0001);
}

TEST(LibraryPseudocountsTest, K4000ZincFingerAlignment) {
  FILE* ali_in = fopen("../data/zinc_finger.fas", "r");
  Alignment<AA> ali(ali_in, FASTA_ALIGNMENT);
  fclose(ali_in);
  CountProfile<AA> ali_profile(ali, false);

  FILE* fin = fopen("../data/nr30_neff2.5_1psi_W13_N1000000_K4000.lib", "r");
  ContextLibrary<AA> lib(fin);
  fclose(fin);
  TransformToLog(lib);

  LibraryPseudocounts<AA> pc(lib, 1.6, 0.85);
  Profile<AA> profile(pc.AddTo(ali_profile, CSBlastAdmix(1.0, 10.0)));

  EXPECT_NEAR(0.8001, profile[53][AA::kCharToInt['C']], 0.0001);
  EXPECT_NEAR(0.7961, profile[56][AA::kCharToInt['C']], 0.0001);
}

}  // namespace cs
