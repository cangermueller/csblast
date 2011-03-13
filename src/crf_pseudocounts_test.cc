#include <gtest/gtest.h>

#include "cs.h"
#include "crf-inl.h"
#include "crf_pseudocounts-inl.h"
#include "sequence-inl.h"

namespace cs {

TEST(CrfPseudocountsTest, K50ZincFingerSequence) {
  FILE* seq_in = fopen("../data/zinc_finger.seq", "r");
  Sequence<AA> seq(seq_in);
  fclose(seq_in);

  FILE* fin = fopen("../data/K50.crf", "r");
  Crf<AA> crf(fin);
  fclose(fin);

  CrfPseudocounts<AA> pc(crf);
  Profile<AA> profile(pc.AddTo(seq, ConstantAdmix(1.0)));

  EXPECT_NEAR(0.744567, profile[53][AA::kCharToInt['C']], 0.0001);
  EXPECT_NEAR(0.757459, profile[56][AA::kCharToInt['C']], 0.0001);
}

TEST(CrfPseudocountsTest, K50ZincFingerAlignment) {
  FILE* ali_in = fopen("../data/zinc_finger.fas", "r");
  Alignment<AA> ali(ali_in, FASTA_ALIGNMENT);
  fclose(ali_in);
  CountProfile<AA> ali_profile(ali, false);

  FILE* fin = fopen("../data/K50.crf", "r");
  Crf<AA> crf(fin);
  fclose(fin);

  CrfPseudocounts<AA> pc(crf);
  Profile<AA> profile(pc.AddTo(ali_profile, CSBlastAdmix(1.0, 10.0)));

  EXPECT_NEAR(0.744527, profile[53][AA::kCharToInt['C']], 0.0001);
  EXPECT_NEAR(0.740898, profile[56][AA::kCharToInt['C']], 0.0001);
}

}  // namespace cs
