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

    EXPECT_EQ(6, profile.num_cols());
    EXPECT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
    EXPECT_FLOAT_EQ(0.0f, profile[1][0]);
}

TEST(ProfileTest, ConstructionOfMultipleProfilesFromInputStream)
{
    FILE* fin = fopen("../data/nt_profiles.prf", "r");
    std::vector< shared_ptr<Profile<Nucleotide> > > profiles;
    Profile<Nucleotide>::readall(fin, &profiles);
    fclose(fin);

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
    FILE* fin = fopen("../data/nt_profile.prf", "r");
    Profile<Nucleotide> profile(fin);
    fclose(fin);
    Profile<Nucleotide> subprofile(profile, 1, 2);

    EXPECT_EQ(2, subprofile.num_cols());
    EXPECT_EQ(4, subprofile.alphabet_size());
    EXPECT_FLOAT_EQ(0.0f, subprofile[0][0]);
    EXPECT_FLOAT_EQ(1.0f, subprofile[0][1]);
}

TEST(ProfileTest, Log2LinSpace)
{
    FILE* fin = fopen("../data/nt_log_profile.prf", "r");
    Profile<Nucleotide> profile(fin);
    fclose(fin);

    ASSERT_EQ(6, profile.num_cols());
    ASSERT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);

    profile.transform_to_linspace();

    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);
}

TEST(ProfileTest, Lin2LogSpace)
{
    FILE* fin = fopen("../data/nt_profile.prf", "r");
    Profile<Nucleotide> profile(fin);
    fclose(fin);

    ASSERT_EQ(6, profile.num_cols());
    ASSERT_EQ(4, profile.alphabet_size());
    EXPECT_FLOAT_EQ(1.0f, profile[0][0]);

    profile.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, profile[0][0]);
}

}  // cs
