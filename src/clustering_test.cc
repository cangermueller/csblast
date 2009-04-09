#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "alignment-inl.h"
#include "amino_acid.h"
#include "clustering-inl.h"
#include "blosum_matrix.h"
#include "count_profile-inl.h"
#include "log.h"
#include "matrix_pseudocounts-inl.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs
{

const float kFloatDelta = 0.01;

class ClusteringTest : public testing::Test
{
  protected:
    static const int WINDOW_LENGTH = 5;
    static const int NUM_STATES    = 10;

    ClusteringTest()
            : lib_(NUM_STATES, WINDOW_LENGTH)
    { }

    virtual void SetUp()
    {
        // Build generic zinc finger HMM
        Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
        Profile<AminoAcid> profile(seq.length());

        BlosumMatrix m;
        MatrixPseudocounts<AminoAcid> mpc(&m);
        mpc.add_to_sequence(seq, ConstantAdmixture(1.0f), &profile);

        // Initialize profile library
        for (int i = 0; i < NUM_STATES; ++i) {
            Profile<AminoAcid> p(profile, i, WINDOW_LENGTH);
            lib_.add_profile(p);
        }
        lib_.transform_to_logspace();

        // Read zinc finger alignments
        std::ifstream ali_in("../data/zinc_finger_alignments.fas");
        std::vector< shared_ptr< Alignment<AminoAcid> > > alis;
        Alignment<AminoAcid>::readall(ali_in, Alignment<AminoAcid>::FASTA, alis);
        ali_in.close();

        // Convert alignments to counts and add pseudocounts
        for (std::vector< shared_ptr< Alignment<AminoAcid> > >::iterator ai = alis.begin(); ai != alis.end(); ++ai) {
            CountProfile<AminoAcid> p_full(**ai, true);
            for (int i = 0; i < p_full.num_cols() - WINDOW_LENGTH + 1; ++i) {
                shared_ptr< CountProfile<AminoAcid> > p_ptr(new CountProfile<AminoAcid>(p_full, i, WINDOW_LENGTH));
                mpc.add_to_profile(ConstantAdmixture(0.01f), p_ptr.get());
                p_ptr->convert_to_counts();
                counts_.push_back(p_ptr);
            }
        }
    }

    ProfileLibrary<AminoAcid> lib_;
    std::vector< shared_ptr< CountProfile<AminoAcid> > > counts_;
};

TEST_F(ClusteringTest, ZincFingerAlisClustering)
{
    ClusteringParams p;
    p.num_blocks   = 1;
    p.epsilon_null = 1.0f;
    p.beta         = 0.0f;
    p.min_scans    = 5;
    Clustering<AminoAcid, CountProfile> em_clust(p, counts_, lib_, std::cout);
    em_clust.run();

    lib_.transform_to_linspace();
    EXPECT_NEAR(0.9948, lib_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
    EXPECT_NEAR(0.4975, lib_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

TEST_F(ClusteringTest, ZincFingerAlisOnlineClustering)
{
    ClusteringParams p;
    p.num_blocks   = 2;
    p.epsilon_null = 0.05f;
    p.beta         = 0.5f;
    p.min_scans    = 5;
    Clustering<AminoAcid, CountProfile> em_clust(p, counts_, lib_, std::cout);
    em_clust.run();

    lib_.transform_to_linspace();
    EXPECT_NEAR(0.9948, lib_[0][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
    EXPECT_NEAR(0.4975, lib_[5][0][AminoAcid::instance().ctoi('C')], kFloatDelta);
}

}  // namespace cs
