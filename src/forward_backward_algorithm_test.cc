#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "forward_backward_algorithm.h"
#include "hmm.h"
#include "log.h"
#include "matrix_pseudocounts.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

const float DELTA = 0.0001;

class ForwardBackwardAlgorithmTest : public testing::Test
{
  protected:
    ForwardBackwardAlgorithmTest()
            : hmm_(13, 1)
    { }

    virtual void SetUp()
    {
        Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
        Profile<AminoAcid> profile(seq.length());

        BlosumMatrix m;
        MatrixPseudocounts<AminoAcid> mpc(&m);
        mpc.add_to_sequence(seq, profile, ConstantAdmixture(0.1f));

        for (int i = 0; i < seq.length(); ++i) {
            Profile<AminoAcid> p(profile, i, 1);
            hmm_.add_profile(p);
        }
        hmm_.init_transitions(HomogeneousTransitionInitializer<AminoAcid>());
        hmm_.transform_states_to_logspace();
    }

    HMM<AminoAcid> hmm_;
};

TEST_F(ForwardBackwardAlgorithmTest, ZincFingerMotif)
{
    Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
    ForwardBackwardAlgorithm<AminoAcid, Sequence> fb;
    shared_ptr<ForwardBackwardMatrices> m = fb.run(hmm_, seq);

    EXPECT_NEAR(0.9566, m->f[0][0] * m->b[0][0], DELTA);
    EXPECT_NEAR(0.4920, m->f[1][1] * m->b[1][1], DELTA);
    EXPECT_NEAR(0.9326, m->f[2][2] * m->b[2][2], DELTA);
}

TEST_F(ForwardBackwardAlgorithmTest, 1Q7L)
{
    std::ifstream fin("../data/1Q7L.fas");
    Alignment<AminoAcid> alignment(fin, Alignment<AminoAcid>::FASTA);
    fin.close();
    alignment.assign_match_columns_by_gap_rule();
    CountsProfile<AminoAcid> profile(alignment, true);

    BlosumMatrix m;
    MatrixPseudocounts<AminoAcid> mpc(&m);
    mpc.add_to_profile(profile, ConstantAdmixture(0.01f));

    HMM<AminoAcid> hmm(profile.length(), 1);
    for (int i = 0; i < profile.length(); ++i) {
            Profile<AminoAcid> p(profile, i, 1);
            hmm.add_profile(p);
    }
    hmm.init_transitions(HomogeneousTransitionInitializer<AminoAcid>());
    hmm.transform_states_to_logspace();

    profile.convert_to_counts();
    ForwardBackwardAlgorithm<AminoAcid, CountsProfile> fb;
    LOG(INFO) << hmm.num_states();
    LOG(INFO) << "start";
    shared_ptr<ForwardBackwardMatrices> mat = fb.run(hmm, profile);
    LOG(INFO) << "end";
}

}  // cs
