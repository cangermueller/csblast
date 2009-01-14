#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"
#include "sequence.h"
#include "smart_ptr.h"

TEST(SequenceTest, ConstructionFromAlphabetVector)
{
    cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();

    std::string header("dummy sequence header");
    std::string seq(aa->begin(), aa->end());
    seq.insert(seq.begin()+1, ' ');
    seq.push_back('\n');
    std::vector<char> intvec;
    for (cs::AminoAcidAlphabet::const_iterator iter = aa->begin(); iter != aa->end(); ++iter)
        intvec.push_back(aa->ctoi(*iter));

    const cs::Sequence sequence(header, seq, aa);

    EXPECT_EQ(aa->size(), sequence.length());
    EXPECT_EQ(aa->ctoi('R'), sequence(1));
    EXPECT_TRUE(equal(sequence.begin(), sequence.end(), intvec.begin()));
}

TEST(SequenceTest, ConstructionFromInputStream)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::istringstream ss(">dummy header\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC");
    cs::Sequence sequence(ss, na);

    EXPECT_EQ(80, sequence.length());
    EXPECT_EQ(na->ctoi('C'), sequence(1));
    EXPECT_EQ(na->ctoi('C'), sequence(79));
}

TEST(SequenceTest, ConstructionOfMultipleSequencesFromInputStream)
{
    cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
    std::string data;
    data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
    std::istringstream ss(data);
    std::vector< SmartPtr<cs::Sequence> > seqs(cs::Sequence::read(ss, na));

    EXPECT_EQ(2, static_cast<int>(seqs.size()));
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])(1));
    EXPECT_EQ(na->ctoi('C'), (*seqs[0])(79));
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])(1));
    EXPECT_EQ(na->ctoi('C'), (*seqs[1])(79));
}
