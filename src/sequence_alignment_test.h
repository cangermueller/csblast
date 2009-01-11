#include <cxxtest/TestSuite.h>

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#include "nucleic_acid_alphabet.h"
#include "sequence_alignment.h"

class SequenceAlignmentTestSuite : public CxxTest::TestSuite
{
  public:
    void test_construction_from_input_stream( void )
    {
        cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
        std::string data;
        data.append(">seq1\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
        data.append(">seq2\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTACACGTACGTAC\nACGTACGTACACGTACGTAC\n");
        std::istringstream ss(data);
        cs::SequenceAlignment alignment(ss, na);

        TS_ASSERT_EQUALS( alignment.nseqs(), 2 );
        TS_ASSERT_EQUALS( alignment.ncols(), 80 );
        TS_ASSERT_EQUALS( alignment(0,0), na->ctoi('A') );
        TS_ASSERT_EQUALS( alignment(1,1), na->ctoi('C') );
    }
};
