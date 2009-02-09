#include <gtest/gtest.h>

#include <fstream>
#include <iomanip>
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
    Sequence<AminoAcid> seq("zinc finger motif", "CRIC\n");
    ForwardBackwardAlgorithm<AminoAcid, Sequence> fb;
    shared_ptr<ForwardBackwardMatrices> fb_mat = fb.run(hmm_, seq);

    for (int i = 1; i <= seq.length(); ++i) {
        for (int k = 1; k <= hmm_.num_states(); ++k) {
            //            double p = ( fb_mat->f[i][k] * fb_mat->b[i][k] ) / fb_mat->p_forward;
            double p = fb_mat->f[i][k];
            std::cout << std::setw(8) << std::right << std::fixed << std::setprecision(5) << p;
        }
        std::cout << std::endl;
    }

    EXPECT_EQ(195, hmm_.num_transitions());
    EXPECT_EQ(13, hmm_.num_states());
}

}  // cs
