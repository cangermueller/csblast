#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid_alphabet.h"
#include "nucleotide_alphabet.h"
#include "sequence.h"
#include "shared_ptr.h"

TEST(SequenceTest, ConstructionFromAlphabetVector)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();

    std::string header("dummy sequence header");
    std::string seq("ARNDCQEGHILKMFPSTWYV");
    seq.insert(seq.begin()+1, ' ');
    seq.push_back('\n');

    const cs::Sequence sequence(header, seq, aa);

    EXPECT_EQ(aa->size(), sequence.length());
    EXPECT_EQ(aa->ctoi('R'), sequence[1]);
}

TEST(SequenceTest, ConstructionFromInputStream)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::istringstream ss(">dummy header\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC");
    cs::Sequence sequence(ss, na);

    EXPECT_EQ(80, sequence.length());
    EXPECT_EQ(na->ctoi('C'), sequence[1]);
    EXPECT_EQ(na->ctoi('C'), sequence[79]);
}

TEST(SequenceTest, ConstructionOfMultipleSequencesFromInputStream)
{
    cs::NucleotideAlphabet* na = cs::NucleotideAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    std::vector< shared_ptr<cs::Sequence> > seqs(cs::Sequence::readall(ss, na));

    EXPECT_EQ(2, static_cast<int>(seqs.size()));
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])[1]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])[79]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])[1]);
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])[79]);
}
