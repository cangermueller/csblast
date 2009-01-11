#include <cxxtest/TestSuite.h>

#include <string>
#include <sstream>

#include "amino_acid_alphabet.h"
#include "nucleic_acid_alphabet.h"
#include "sequence_profile.h"

class SequenceProfileTestSuite : public CxxTest::TestSuite
{
  public:
    void test_construction_from_input_stream( void )
    {
        std::string data;
        data.append("#\tA\tC\tG\tT\n");
        data.append("1\t0\t*\t*\t*\n");
        data.append("2\t*\t0\t*\t*\n");
        data.append("3\t*\t*\t0\t*\n");
        data.append("4\t*\t*\t*\t0\n");
        data.append("5\t0\t*\t*\t*\n");
        data.append("6\t*\t0\t*\t*\n");
        std::istringstream ss(data);

        cs::SequenceProfile profile(ss, cs::NucleicAcidAlphabet::instance());

        TS_ASSERT_EQUALS( profile.ncols(), 6 );
        TS_ASSERT_EQUALS( profile.ndim(), 4 );
        TS_ASSERT_EQUALS( profile(0,0), 1.0f );
        TS_ASSERT_EQUALS( profile(1,0), 0.0f );
    }

    void test_construction_multiple_profiles_from_input_stream( void )
    {
        std::string data;
        data.append("#\tA\tC\tG\tT\n");
        data.append("1\t0\t*\t*\t*\n");
        data.append("2\t*\t0\t*\t*\n");
        data.append("3\t*\t*\t0\t*\n");
        data.append("4\t*\t*\t*\t0\n");
        data.append("5\t0\t*\t*\t*\n");
        data.append("6\t*\t0\t*\t*\n");
        data.append("//\n");
        std::istringstream ss(data+data);

        std::vector<cs::SequenceProfile*> profiles(cs::SequenceProfile::read(ss, cs::NucleicAcidAlphabet::instance()));

        TS_ASSERT_EQUALS( static_cast<int>(profiles.size()), 2 );
        TS_ASSERT_EQUALS( profiles[0]->ncols(), 6 );
        TS_ASSERT_EQUALS( profiles[0]->ndim(), 4 );
        TS_ASSERT_EQUALS( profiles[1]->ncols(), 6 );
        TS_ASSERT_EQUALS( profiles[1]->ndim(), 4 );
        TS_ASSERT_EQUALS( (*profiles[0])(0,0), 1.0f );
        TS_ASSERT_EQUALS( (*profiles[0])(1,0), 0.0f );
        TS_ASSERT_EQUALS( (*profiles[1])(0,0), 1.0f );
        TS_ASSERT_EQUALS( (*profiles[1])(1,0), 0.0f );

        for (std::vector<cs::SequenceProfile*>::iterator iter = profiles.begin(); iter != profiles.end(); ++iter)
            delete *iter;
    }

    void test_construction_from_input_stream_with_improper_alphabet( void )
    {
        std::string data;
        data.append("#\tA\tC\tG\tT\n");
        data.append("1\t0\t*\t*\t*\n");
        data.append("2\t*\t0\t*\t*\n");
        data.append("3\t*\t*\t0\t*\n");
        data.append("4\t*\t*\t*\t0\n");
        data.append("5\t0\t*\t*\t*\n");
        data.append("6\t*\t0\t*\t*\n");
        std::istringstream ss(data);

        TS_ASSERT_THROWS_ANYTHING( cs::SequenceProfile(ss, cs::AminoAcidAlphabet::instance()) );
    }

    void test_construction_of_subprofile( void )
    {
        std::string data;
        data.append("#\tA\tC\tG\tT\n");
        data.append("1\t0\t*\t*\t*\n");
        data.append("2\t*\t0\t*\t*\n");
        data.append("3\t*\t*\t0\t*\n");
        std::istringstream ss(data);

        cs::SequenceProfile profile(ss, cs::NucleicAcidAlphabet::instance());
        cs::SequenceProfile subprofile(profile, 1, 2);

        TS_ASSERT_EQUALS( subprofile.ncols(), 2 );
        TS_ASSERT_EQUALS( subprofile.ndim(), 4 );
        TS_ASSERT_EQUALS( subprofile(0,0), 0.0f );
        TS_ASSERT_EQUALS( subprofile(0,1), 1.0f );
    }
};
