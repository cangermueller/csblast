#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "matrix_pseudocounts-inl.h"
#include "sequence-inl.h"
#include "training_sequence.h"

namespace cs {

TEST(SequenceTest, ConstructionFromAlphabet) {
  Sequence<AA> sequence("A --RNDCQEGHILKMFPSTWYV\n", "header");

  EXPECT_EQ(AA::kSize, sequence.length());
  EXPECT_EQ(AA::kCharToInt['R'], sequence[1]);
  EXPECT_EQ(std::string("header"), sequence.header());
}

TEST(SequenceTest, ReadDnaSeqFile) {
  FILE* fin = fopen("../data/seq00.fas", "r");
  Sequence<Dna> sequence(fin);
  fclose(fin);

  EXPECT_EQ(Dna::kCharToInt['T'], sequence[0]);
  EXPECT_EQ(Dna::kCharToInt['A'], sequence[1]);
  EXPECT_EQ(Dna::kCharToInt['A'], sequence[2]);
  EXPECT_EQ(Dna::kCharToInt['T'], sequence[3]);
}

TEST(SequenceTest, ReadAS62) {
  FILE* fin = fopen("../data/AS62.seq", "r");
  Sequence<AS62> sequence(fin);
  fclose(fin);

  EXPECT_EQ(AS62::kCharToInt['X'], sequence[0]);
  EXPECT_EQ(AS62::kCharToInt['Z'], sequence[1]);
  EXPECT_EQ(AS62::kCharToInt['L'], sequence[2]);
  EXPECT_EQ(AS62::kCharToInt['5'], sequence[3]);
}

TEST(SequenceTest, ReadAS90) {
  FILE* fin = fopen("../data/AS90.seq", "r");
  Sequence<AS90> sequence(fin);
  fclose(fin);

  EXPECT_EQ(AS90::kCharToInt['%'], sequence[0]);
  EXPECT_EQ(AS90::kCharToInt['i'], sequence[1]);
  EXPECT_EQ(AS90::kCharToInt['O'], sequence[2]);
  EXPECT_EQ(AS90::kCharToInt[';'], sequence[3]);
}

TEST(SequenceTest, AddMatrixPseudocountsToSequence) {
  const Sequence<AA> sequence("ARNDCQEGHILKMFPSTWYV");

  ASSERT_EQ(AA::kSize, sequence.length());
  ASSERT_EQ(AA::kCharToInt['R'], sequence[1]);

  BlosumMatrix m;
  MatrixPseudocounts<AA> mpc(m);
  Profile<AA> p(mpc.AddTo(sequence, CSBlastAdmix(1.0, 10.0)));

  EXPECT_NEAR(0.06f, p[0][AA::kCharToInt['V']], 0.01);
}

}  // namespace cs
