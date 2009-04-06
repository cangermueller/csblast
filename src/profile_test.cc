#include <gtest/gtest.h>

#include <cstdio>

#include <string>
#include <sstream>

#include "amino_acid.h"
#include "nucleotide.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs
{

TEST(ProfileTest, ConstructionFromInputStream)
{
    FILE* fin = fopen("../data/nt_profile.prf", "r");
    Profile<Nucleotide> profile(fin);
    fclose(fin);

    // EXPECT_EQ(6, profile.num_cols());
    // EXPECT_EQ(4, profile.alphabet_size());
    // EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    // EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
}

TEST(ProfileTest, ConstructionOfMultipleProfilesFromInputStream)
{
    std::string data;
    data.append("Profile\n");
    data.append("num_cols\t\t6\n");
    data.append("alphabet_size\t4\n");
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

    std::vector< shared_ptr<Profile<Nucleotide> > > profiles;
    Profile<Nucleotide>::readall(ss, profiles);

    EXPECT_EQ(2, static_cast<int>(profiles.size()));
    EXPECT_EQ(6, (*profiles[0]).num_cols());
    EXPECT_EQ(4, (*profiles[0]).alphabet_size());
    EXPECT_EQ(6, (*profiles[1]).num_cols());
    EXPECT_EQ(4, (*profiles[1]).alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, (*profiles[0])[0][0]);
    EXPECT_FLOAT_EQ(0.0f, (*profiles[0])[1][0]);
    EXPECT_FLOAT_EQ(1.0f, (*profiles[1])[0][0]);
    EXPECT_FLOAT_EQ(0.0f, (*profiles[1])[1][0]);
}

TEST(ProfileTest, ConstructionOfSubprofile)
{
    std::string data;
    data.append("Profile\n");
    data.append("num_cols\t\t3\n");
    data.append("alphabet_size\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    std::istringstream ss(data);

    Profile<Nucleotide> profile(ss);
    Profile<Nucleotide> subprofile(profile, 1, 2);

    EXPECT_EQ(2, subprofile.num_cols());
    EXPECT_EQ(4, subprofile.alphabet_size());
    EXPECT_FLOAT_EQ(0.0f, subprofile[0][0]);
    EXPECT_FLOAT_EQ(1.0f, subprofile[0][1]);
}

TEST(ProfileTest, Log2LinSpace)
{
    std::string data;
    data.append("Profile\n");
    data.append("num_cols\t\t6\n");
    data.append("alphabet_size\t4\n");
    data.append("logspace\t1\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    std::istringstream ss(data);

    Profile<Nucleotide> profile(ss);

    ASSERT_EQ(6, profile.num_cols());
    ASSERT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);

    profile.transform_to_linspace();

    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
}

TEST(ProfileTest, Lin2LogSpace)
{
    std::string data;
    data.append("Profile\n");
    data.append("num_cols\t\t6\n");
    data.append("alphabet_size\t4\n");
    data.append("logspace\t0\n");
    data.append(" \tA\tC\tG\tT\n");
    data.append("1\t0\t*\t*\t*\n");
    data.append("2\t*\t0\t*\t*\n");
    data.append("3\t*\t*\t0\t*\n");
    data.append("4\t*\t*\t*\t0\n");
    data.append("5\t0\t*\t*\t*\n");
    data.append("6\t*\t0\t*\t*\n");
    std::istringstream ss(data);

    Profile<Nucleotide> profile(ss);

    ASSERT_EQ(6, profile.num_cols());
    ASSERT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);
}

}  // cs
