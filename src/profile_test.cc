#include <gtest/gtest.h>

#include <string>
#include <sstream>

#include "amino_acid_alphabet.h"
#include "nucleotide_alphabet.h"
#include "profile.h"
#include "shared_ptr.h"

TEST(ProfileTest, ConstructionFromInputStream)
{
    std::string data;
    data.append("Profile\n");
    data.append("ncols\t6\n");
    data.append("nalph\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    std::istringstream ss(data);

    cs::Profile profile(ss, cs::NucleotideAlphabet::instance());

    EXPECT_EQ(6, profile.ncols());
    EXPECT_EQ(4, profile.nalph());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
}

TEST(ProfileTest, ConstructionOfMultipleProfilesFromInputStream)
{
    std::string data;
    data.append("Profile\n");
    data.append("ncols\t6\n");
    data.append("nalph\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    data.append("//\n");
    std::istringstream ss(data+data);

    std::vector< shared_ptr<cs::Profile> > profiles(cs::Profile::readall(ss, cs::NucleotideAlphabet::instance()));

    EXPECT_EQ(2, static_cast<int>(profiles.size()));
    EXPECT_EQ(6, (*profiles[0]).ncols());
    EXPECT_EQ(4, (*profiles[0]).nalph());
    EXPECT_EQ(6, (*profiles[1]).ncols());
    EXPECT_EQ(4, (*profiles[1]).nalph());
    EXPECT_FLOAT_EQ(1.0f, (*profiles[0])[0][0]);
    EXPECT_FLOAT_EQ(0.0f, (*profiles[0])[1][0]);
    EXPECT_FLOAT_EQ(1.0f, (*profiles[1])[0][0]);
    EXPECT_FLOAT_EQ(0.0f, (*profiles[1])[1][0]);
}

TEST(ProfileTest, ConstructionOfSubprofile)
{
    std::string data;
    data.append("Profile\n");
    data.append("ncols\t\t3\n");
    data.append("nalph\t\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    std::istringstream ss(data);

    cs::Profile profile(ss, cs::NucleotideAlphabet::instance());
    cs::Profile subprofile(profile, 1, 2);

    EXPECT_EQ(2, subprofile.ncols());
    EXPECT_EQ(4, subprofile.nalph());
    EXPECT_FLOAT_EQ(0.0f, subprofile[0][0]);
    EXPECT_FLOAT_EQ(1.0f, subprofile[0][1]);
}

TEST(ProfileTest, Log2LinSpace)
{
    std::string data;
    data.append("Profile\n");
    data.append("ncols\t\t6\n");
    data.append("nalph\t\t4\n");
    data.append("logspace\t1\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    std::istringstream ss(data);

    cs::Profile profile(ss, cs::NucleotideAlphabet::instance());

    ASSERT_EQ(6, profile.ncols());
    ASSERT_EQ(4, profile.nalph());
    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);

    profile.transform_to_linspace();

    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);

    //    profile.write(std::cout);
}

TEST(ProfileTest, Lin2LogSpace)
{
    std::string data;
    data.append("Profile\n");
    data.append("ncols\t\t6\n");
    data.append("nalph\t\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    std::istringstream ss(data);

    cs::Profile profile(ss, cs::NucleotideAlphabet::instance());

    ASSERT_EQ(6, profile.ncols());
    ASSERT_EQ(4, profile.nalph());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);
}
