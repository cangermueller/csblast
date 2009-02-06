#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "alignment.h"
#include "amino_acid.h"
#include "blosum_matrix.h"
#include "hmm.h"
#include "log.h"
#include "matrix_pseudocounts.h"
#include "nucleotide.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs
{

class HMMTest : public testing::Test
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

TEST_F(HMMTest, SimpleConstruction)
{
    HMM<Nucleotide> hmm(3);

    EXPECT_EQ(0, hmm.num_states());
    EXPECT_EQ(0, hmm.num_transitions());
    EXPECT_EQ(3, hmm.size());

    int index_p1 = hmm.add_profile(p1_);
    EXPECT_EQ(1, index_p1);
    int index_p2 = hmm.add_profile(p2_);
    EXPECT_EQ(2, index_p2);
    int index_p3 = hmm.add_profile(p3_);
    EXPECT_EQ(3, index_p3);

    EXPECT_FLOAT_EQ(1.0f, hmm[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[1][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[2][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[3][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[3][1][0]);

    hmm[3][0][0] = 0.0f;
    hmm[3][0][1] = 1.0f;
    EXPECT_FLOAT_EQ(0.0f, hmm[3][0][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[3][0][1]);

    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(0, 1));
    hmm.set_transition_probability(0, 1, 0.5f);
    EXPECT_EQ(1, hmm[0].num_out_transitions());
    EXPECT_EQ(0, hmm[0].num_in_transitions());
    EXPECT_EQ(1, hmm[1].num_in_transitions());
    EXPECT_EQ(0, hmm[1].num_out_transitions());
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 1));
    EXPECT_FLOAT_EQ(0.5f, hmm[0].transition_probability_to(1));
    EXPECT_FLOAT_EQ(0.5f, hmm[1].transition_probability_from(0));
    EXPECT_FLOAT_EQ(0.5f, (**hmm[0].out_transitions_begin()).probability);
    EXPECT_FLOAT_EQ(0.5f, (**hmm[1].in_transitions_begin()).probability);

    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(0, 2));
    hmm.set_transition_probability(0, 2, 0.5f);
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 2));
    EXPECT_EQ(2, hmm[0].num_out_transitions());
    EXPECT_EQ(1, hmm[2].num_in_transitions());

    hmm.set_transition_probability(1, 3, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(1, 3));
    EXPECT_EQ(1, hmm[1].num_out_transitions());
    EXPECT_EQ(1, hmm[3].num_in_transitions());

    hmm.set_transition_probability(2, 3, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(2, 3));
    EXPECT_EQ(1, hmm[2].num_out_transitions());
    EXPECT_EQ(2, hmm[3].num_in_transitions());

    hmm.set_transition_probability(3, 0, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(3, 0));
    EXPECT_EQ(1, hmm[3].num_out_transitions());
    EXPECT_EQ(1, hmm[0].num_in_transitions());

    EXPECT_EQ(5, hmm.num_transitions());
    EXPECT_EQ(3, hmm.num_states());
}

TEST_F(HMMTest, ConstructionFromSerializedHMM)
{
    HMM<Nucleotide> hmm1(3);
    hmm1.add_profile(p1_);
    hmm1.add_profile(p2_);
    hmm1.add_profile(p3_);
    hmm1.set_transition_probability(0, 1, 0.5f);
    hmm1.set_transition_probability(0, 2, 0.5f);
    hmm1.set_transition_probability(1, 3, 1.0f);
    hmm1.set_transition_probability(2, 3, 1.0f);
    hmm1.set_transition_probability(3, 0, 1.0f);

    std::ostringstream out;
    hmm1.write(out);
    std::istringstream in(out.str());

    HMM<Nucleotide> hmm2(in);

    EXPECT_EQ(hmm1.num_transitions(), hmm2.num_transitions());
    EXPECT_EQ(hmm1.num_states(), hmm2.num_states());
    EXPECT_FLOAT_EQ(hmm1[1][0][0], hmm2[1][0][0]);
    EXPECT_FLOAT_EQ(hmm1[1][1][0], hmm2[1][1][0]);
    EXPECT_FLOAT_EQ(hmm1[2][0][0], hmm2[2][0][0]);
    EXPECT_FLOAT_EQ(hmm1[2][1][0], hmm2[2][1][0]);
    EXPECT_FLOAT_EQ(hmm1[3][0][0], hmm2[3][0][0]);
    EXPECT_FLOAT_EQ(hmm1[3][1][0], hmm2[3][1][0]);
    EXPECT_EQ(hmm1[2].num_out_transitions(), hmm2[2].num_out_transitions());
    EXPECT_EQ(hmm1[3].num_in_transitions(), hmm2[3].num_in_transitions());
    EXPECT_FLOAT_EQ(hmm1.transition_probability(3, 0), hmm2.transition_probability(3, 0));
}

TEST_F(HMMTest, NormalizeTransitions)
{
    HMM<Nucleotide> hmm(3);
    hmm.add_profile(p1_);
    hmm.add_profile(p2_);
    hmm.add_profile(p3_);
    hmm.set_transition_probability(0, 1, 0.3f);
    hmm.set_transition_probability(0, 2, 0.3f);
    hmm.set_transition_probability(1, 3, 0.9f);
    hmm.set_transition_probability(2, 3, 0.7f);

    normalize_transitions(hmm);

    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 1));
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 2));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(1, 3));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(2, 3));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(3, 0));
}

TEST_F(HMMTest, LinLogTransformation)
{
    HMM<Nucleotide> hmm(3);
    hmm.add_profile(p1_);
    hmm.add_profile(p2_);
    hmm.add_profile(p3_);
    hmm.set_transition_probability(0, 1, 0.5f);
    hmm.set_transition_probability(0, 2, 0.5f);
    hmm.set_transition_probability(1, 3, 1.0f);
    hmm.set_transition_probability(2, 3, 1.0f);
    hmm.set_transition_probability(3, 0, 1.0f);

    hmm.transform_transitions_to_logspace();
    hmm.transform_states_to_logspace();

    EXPECT_FLOAT_EQ(0.0f, hmm[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[3][0][0]);

    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(1, 3));
    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(2, 3));
    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(3, 0));

    hmm.transform_transitions_to_linspace();
    hmm.transform_states_to_linspace();

    EXPECT_FLOAT_EQ(1.0f, hmm[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[1][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[2][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[3][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[3][1][0]);

    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 1));
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 2));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(1, 3));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(2, 3));
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(3, 0));
}

TEST(HMMTestInitialization, RandomSampleInitializer)
{
    std::ifstream fin("../data/1Q7L.fas");
    Alignment<AminoAcid> ali(fin, Alignment<AminoAcid>::FASTA);
    fin.close();
    ali.assign_match_columns_by_gap_rule();

    const int HMM_SIZE = 10;
    const int NUM_COLS = 5;
    const float PC_ADMIXTURE = 0.6f;
    const float SAMPLE_RATE = 0.2f;

    // setup training profiles
    typedef shared_ptr< CountsProfile<AminoAcid> > profile_ptr;
    profile_ptr p(new CountsProfile<AminoAcid>(ali, true));
    std::vector<profile_ptr> profiles(10, p);

    // setup substitution matrix pseudocounts
    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(m, shared_ptr<Admixture>(new ConstantAdmixture(PC_ADMIXTURE)));

    RandomSampleStateInitializer<AminoAcid> st_init(profiles, NUM_COLS, SAMPLE_RATE, mpc);
    RandomTransitionInitializer<AminoAcid> tr_init;
    HMM<AminoAcid> hmm(HMM_SIZE, st_init, tr_init);
    sparsify(hmm, 1.0f / HMM_SIZE);

    EXPECT_EQ(HMM_SIZE, hmm.num_states());
}

}  // cs