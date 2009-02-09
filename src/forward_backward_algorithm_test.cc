#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid.h"
#include "blosum_matrix.h"
#include "forward_backward_algorithm.h"
#include "hmm.h"
#include "log.h"
#include "matrix_pseudocounts.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"

namespace cs
{

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
        ConstantAdmixture pca(0.3f);
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
    ForwardBackwardAlgorithm<AminoAcid, Sequence> fb = ForwardBackwardParams().ignore_profile_context(true);
    shared_ptr<ForwardBackwardMatrices> fb_mat = fb.run(hmm_, seq);

    EXPECT_EQ(195, hmm_.num_transitions());
    EXPECT_EQ(13, hmm_.num_states());
}

}  // cs
