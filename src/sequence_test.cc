#include <gtest/gtest.h>

#include <cstdio>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "nucleotide.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "shared_ptr.h"

namespace cs {

const float kFloatDelta = 0.01f;

TEST(SequenceTest, ConstructionFromAlphabetVector) {
  Sequence<AminoAcid> sequence("header", "A RNDCQEGHILKMFPSTWYV\n");

  EXPECT_EQ(AminoAcid::instance().size(), sequence.length());
  EXPECT_EQ(AminoAcid::instance().ctoi('R'), sequence[1]);
}

TEST(SequenceTest, ConstructionFromInputStream) {
  FILE* fin = fopen("../data/nt_seqs.fas", "r");
  Sequence<Nucleotide> sequence(fin);
  fclose(fin);

  EXPECT_EQ(80, sequence.length());
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), sequence[1]);
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), sequence[79]);
}

TEST(SequenceTest, ConstructionOfMultipleSequencesFromInputStream) {
  FILE* fin = fopen("../data/nt_seqs.fas", "r");
  std::vector< shared_ptr<Sequence<Nucleotide> > > seqs;
  Sequence<Nucleotide>::readall(fin, &seqs);
  fclose(fin);

  EXPECT_EQ(2, static_cast<int>(seqs.size()));
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[0])[1]);
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[0])[79]);
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[1])[1]);
  EXPECT_EQ(Nucleotide::instance().ctoi('C'), (*seqs[1])[79]);
}

TEST(SequenceTest, AddMatrixPseudocountsToSequence) {
  const Sequence<AminoAcid> sequence("header", "ARNDCQEGHILKMFPSTWYV");
  Profile<AminoAcid> profile(sequence.length());

  ASSERT_EQ(AminoAcid::instance().size(), sequence.length());
  ASSERT_EQ(AminoAcid::instance().ctoi('R'), sequence[1]);
  ASSERT_EQ(sequence.length(), profile.num_cols());

  BlosumMatrix m;
  MatrixPseudocounts<AminoAcid> mpc(&m);
  mpc.add_to_sequence(sequence,
                      DivergenceDependentAdmixture(1.0f, 10.0f),
                      &profile);

  EXPECT_NEAR(0.06f, profile[0][AminoAcid::instance().ctoi('V')], kFloatDelta);
}

TEST(DISABLED_SequenceTest, ReadSwissprotDatabaseCppStyle) {
  std::ifstream seq_in("../../../databases/uniprot_sprot.fasta");
  for (int i = 0; i < 200000; ++i) {
    const Sequence<AminoAcid> sequence(seq_in);
    EXPECT_LT(0, sequence.length());
  }
  seq_in.close();
}

TEST(DISABLED_SequenceTest, ReadSwissprotDatabaseCStyle) {
  FILE* fin = fopen("../../../databases/uniprot_sprot.fasta", "r");
  for (int i = 0; i < 200000; ++i) {
    const Sequence<AminoAcid> sequence(fin);
    EXPECT_LT(0, sequence.length());
  }
  fclose(fin);
}

}  // namespace cs
