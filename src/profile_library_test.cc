#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "alignment.h"
#include "amino_acid.h"
#include "blosum_matrix.h"
#include "profile_library.h"
#include "log.h"
#include "matrix_pseudocounts.h"
#include "nucleotide.h"
#include "profile-inl.h"
#include "shared_ptr.h"

namespace cs
{

class ProfileLibraryTest : public testing::Test
{
  protected:

    virtual void SetUp()
    {
        std::string data;
        data.append("Profile\n");
        data.append("num_cols\t\t5\n");
        data.append("alphabet_size\t4\n");
        data.append("logspace\t0\n");
        data.append(" \tA\tC\tG\tT\n");
        data.append("1\t0\t*\t*\t*\n");
        data.append("2\t*\t0\t*\t*\n");
        data.append("3\t*\t*\t0\t*\n");
        data.append("4\t*\t*\t*\t0\n");
        data.append("5\t0\t*\t*\t*\n");

        std::istringstream ss1(data);
        p1_.read(ss1);
        std::istringstream ss2(data);
        p2_.read(ss2);
        std::istringstream ss3(data);
        p3_.read(ss3);
    }

    Profile<Nucleotide> p1_;
    Profile<Nucleotide> p2_;
    Profile<Nucleotide> p3_;
};

TEST_F(ProfileLibraryTest, SimpleConstruction)
{
    ProfileLibrary<Nucleotide> lib(3, 5);

    EXPECT_FALSE(lib.full());
    EXPECT_EQ(3, lib.num_profiles());
    EXPECT_EQ(5, lib.num_cols());

    int index_p1 = lib.add_profile(p1_);
    EXPECT_EQ(0, index_p1);
    int index_p2 = lib.add_profile(p2_);
    EXPECT_EQ(1, index_p2);
    int index_p3 = lib.add_profile(p3_);
    EXPECT_EQ(2, index_p3);

    EXPECT_FLOAT_EQ(1.0f, lib[0][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[0][1][0]);
    EXPECT_FLOAT_EQ(1.0f, lib[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[1][1][0]);
    EXPECT_FLOAT_EQ(1.0f, lib[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[2][1][0]);

    lib[2][0][0] = 0.0f;
    lib[2][0][1] = 1.0f;
    EXPECT_FLOAT_EQ(0.0f, lib[2][0][0]);
    EXPECT_FLOAT_EQ(1.0f, lib[2][0][1]);
}

TEST_F(ProfileLibraryTest, ConstructionFromSerializedLibrary)
{
    ProfileLibrary<Nucleotide> lib1(3, 5);
    lib1.add_profile(p1_);
    lib1.add_profile(p2_);
    lib1.add_profile(p3_);

    std::ostringstream out;
    lib1.write(out);
    std::istringstream in(out.str());

    ProfileLibrary<Nucleotide> lib2(in);

    EXPECT_EQ(lib1.num_profiles(), lib2.num_profiles());
    EXPECT_FLOAT_EQ(lib1[0][0][0], lib2[0][0][0]);
    EXPECT_FLOAT_EQ(lib1[0][1][0], lib2[0][1][0]);
    EXPECT_FLOAT_EQ(lib1[1][0][0], lib2[1][0][0]);
    EXPECT_FLOAT_EQ(lib1[1][1][0], lib2[1][1][0]);
    EXPECT_FLOAT_EQ(lib1[2][0][0], lib2[2][0][0]);
    EXPECT_FLOAT_EQ(lib1[2][1][0], lib2[2][1][0]);
}

TEST_F(ProfileLibraryTest, LinLogTransformation)
{
    ProfileLibrary<Nucleotide> lib(3, 5);
    lib.add_profile(p1_);
    lib.add_profile(p2_);
    lib.add_profile(p3_);

    lib.transform_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, lib[0][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[2][0][0]);

    lib.transform_to_linspace();

    EXPECT_FLOAT_EQ(1.0f, lib[0][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[0][1][0]);
    EXPECT_FLOAT_EQ(1.0f, lib[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[1][1][0]);
    EXPECT_FLOAT_EQ(1.0f, lib[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, lib[2][1][0]);
}

TEST(ProfileLibraryTestInitialization, RandomSampleInitializer)
{
    std::ifstream fin("../data/1Q7L.fas");
    Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
    fin.close();
    ali.assign_match_columns_by_gap_rule();

    typedef shared_ptr< CountProfile<AminoAcid> > profile_ptr;
    CountProfile<AminoAcid> p_full(ali, true);
    profile_ptr p_window(new CountProfile<AminoAcid>(p_full, 0, 13));
    std::vector<profile_ptr> profiles(100, p_window);

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> pc(&m);
    SamplingProfileInitializer<AminoAcid> profile_init(profiles, &pc, 0.2f);
    ProfileLibrary<AminoAcid> lib(10, 13, profile_init);

    EXPECT_EQ(10, lib.num_profiles());
    EXPECT_EQ(13, lib.num_cols());
    EXPECT_EQ(13, lib[0].num_cols());
}

}  // cs
