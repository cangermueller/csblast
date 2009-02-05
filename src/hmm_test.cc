#include <gtest/gtest.h>

#include <iostream>
#include <sstream>
#include <string>

#include "hmm.h"
#include "nucleotide.h"
#include "profile.h"

namespace cs
{

class HMMTest : public testing::Test {
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

    int index_p1 = hmm.add_state(p1_);
    EXPECT_EQ(1, index_p1);
    int index_p2 = hmm.add_state(p2_);
    EXPECT_EQ(2, index_p2);
    int index_p3 = hmm.add_state(p3_);
    EXPECT_EQ(3, index_p3);

    EXPECT_FLOAT_EQ(1.0f, hmm[1][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[1][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[2][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[2][1][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[3][0][0]);
    EXPECT_FLOAT_EQ(0.0f, hmm[3][1][0]);

    hmm[3][0][0] = 0.0f;
    hmm[3][1][0] = 1.0f;
    EXPECT_FLOAT_EQ(0.0f, hmm[3][0][0]);
    EXPECT_FLOAT_EQ(1.0f, hmm[3][1][0]);

    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(0, 1));
    hmm.set_transition(0, 1, 0.5f);
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 1));
    EXPECT_EQ(1, hmm[0].num_out_transitions());
    EXPECT_EQ(0, hmm[0].num_in_transitions());
    EXPECT_EQ(1, hmm[1].num_in_transitions());
    EXPECT_EQ(0, hmm[1].num_out_transitions());

    EXPECT_FLOAT_EQ(0.0f, hmm.transition_probability(0, 2));
    hmm.set_transition(0, 2, 0.5f);
    EXPECT_FLOAT_EQ(0.5f, hmm.transition_probability(0, 2));
    EXPECT_EQ(2, hmm[0].num_out_transitions());
    EXPECT_EQ(1, hmm[2].num_in_transitions());

    hmm.set_transition(1, 3, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(1, 3));
    EXPECT_EQ(1, hmm[1].num_out_transitions());
    EXPECT_EQ(1, hmm[3].num_in_transitions());

    hmm.set_transition(2, 3, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(2, 3));
    EXPECT_EQ(1, hmm[2].num_out_transitions());
    EXPECT_EQ(2, hmm[3].num_in_transitions());

    hmm.set_transition(3, 0, 1.0f);
    EXPECT_FLOAT_EQ(1.0f, hmm.transition_probability(3, 0));
    EXPECT_EQ(1, hmm[3].num_out_transitions());
    EXPECT_EQ(1, hmm[0].num_in_transitions());

    EXPECT_EQ(5, hmm.num_transitions());
    EXPECT_EQ(3, hmm.num_states());
}

}  // cs

