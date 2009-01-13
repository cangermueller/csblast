#include <cxxtest/TestSuite.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"
#include "sequence.h"
#include "smart_ptr.h"

class SequenceTestSuite : public CxxTest::TestSuite
{
  public:
    void test_construction_from_character_vector( void )
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

        TS_ASSERT_EQUALS( sequence.length(), aa->size() );
        TS_ASSERT_EQUALS( sequence(1), aa->ctoi('R') );
        TS_ASSERT( equal(sequence.begin(), sequence.end(), intvec.begin()) );
    }

    void test_construction_from_input_stream( void )
    {
        cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
        std::istringstream ss(">dummy header\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC");
        cs::Sequence sequence(ss, na);

        TS_ASSERT_EQUALS( sequence.length(), 80 );
        TS_ASSERT_EQUALS( sequence(1), na->ctoi('C') );
        TS_ASSERT_EQUALS( sequence(79), na->ctoi('C') );
    }

    void test_construction_of_multiple_sequences_from_input_stream( void )
    {
        cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
        std::string data;
        data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
        data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
        std::istringstream ss(data);
        std::vector< SmartPtr<cs::Sequence> > seqs(cs::Sequence::read(ss, na));

        TS_ASSERT_EQUALS( static_cast<int>(seqs.size()), 2 );
        TS_ASSERT_EQUALS( (*seqs[0])(1), na->ctoi('C') );
        TS_ASSERT_EQUALS( (*seqs[0])(79), na->ctoi('C') );
        TS_ASSERT_EQUALS( (*seqs[1])(1), na->ctoi('C') );
        TS_ASSERT_EQUALS( (*seqs[1])(79), na->ctoi('C') );
    }

    void test_construction_from_invalid_character_vector( void )
    {
        cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
        std::string header("dummy sequence header");
        std::string seq(na->begin(), na->end());
        seq.insert(seq.begin()+1, ' ');
        seq.push_back('F'); //invalid character
        seq.push_back('\n');

        TS_ASSERT_THROWS_ANYTHING( cs::Sequence(header, seq, na) );
    }
};
