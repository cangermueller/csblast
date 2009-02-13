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
            : hmm_(13)
    { }

    virtual void SetUp()
    {
        Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
        Profile<AminoAcid> profile(seq.length());

        BlosumMatrix m;
        ConstantAdmixture pca(0.1f);
        MatrixPseudocounts<AminoAcid> mpc(&m, &pca);
        mpc.add_to_sequence(seq, profile);

        for (int i = 0; i < seq.length(); ++i) {
            Profile<AminoAcid> p(profile, i, 1);
            hmm_.add_profile(p);
        }
        hmm_.init_transitions(ConstantTransitionInitializer<AminoAcid>());
        hmm_.transform_states_to_logspace();
    }

    HMM<AminoAcid> hmm_;
};

TEST_F(ForwardBackwardAlgorithmTest, ZincFingerMotif)
{
    Sequence<AminoAcid> seq("zinc finger motif", "GQKPFQCRICMRN\n");
    ForwardBackwardAlgorithm<AminoAcid, Sequence> fb;
    shared_ptr<ForwardBackwardMatrices> m = fb.run(hmm_, seq);

    EXPECT_NEAR(0.9566, m->f[1][1] * m->b[1][1] / m->likelihood, DELTA);
    EXPECT_NEAR(0.4920, m->f[2][2] * m->b[2][2] / m->likelihood, DELTA);
    EXPECT_NEAR(0.9326, m->f[3][3] * m->b[3][3] / m->likelihood, DELTA);
}

}  // cs
