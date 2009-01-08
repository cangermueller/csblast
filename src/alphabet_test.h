#include <cxxtest/TestSuite.h>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"

class AlphabetTestSuite : public CxxTest::TestSuite
{
  public:
    void test_amino_acid_alphabet( void )
    {
        cs::AminoAcidAlphabet* aa = cs::AminoAcidAlphabet::instance();
        TS_ASSERT_EQUALS( aa->size(), 21 );
        TS_ASSERT_EQUALS( aa->itoc( aa->ctoi('R') ), 'R' );
        TS_ASSERT_EQUALS( aa->ctoi( aa->itoc(1) ), 1 );
        TS_ASSERT_EQUALS( *(aa->begin()), 'A' );
    }

    void test_nucleic_acid_alphabet( void )
    {
        cs::NucleicAcidAlphabet* na = cs::NucleicAcidAlphabet::instance();
        TS_ASSERT_EQUALS( na->size(), 5 );
        TS_ASSERT_EQUALS( na->itoc( na->ctoi('C') ), 'C' );
        TS_ASSERT_EQUALS( na->ctoi( na->itoc(1) ), 1 );
        TS_ASSERT_EQUALS( *(na->begin()), 'A' );
    }
};
