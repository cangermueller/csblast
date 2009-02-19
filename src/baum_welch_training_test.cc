#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "amino_acid.h"
#include "baum_welch_training.h"
#include "blosum_matrix.h"
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

class BaumWelchTrainingTest : public testing::Test
{
  protected:
    BaumWelchTrainingTest()
            : hmm_(23)
    { }

    virtual void SetUp()
    {
        Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
        Profile<AminoAcid> profile(seq.length());

        BlosumMatrix m;
        ConstantAdmixture pca(1.0f);
        MatrixPseudocounts<AminoAcid> mpc(&m, &pca);
        mpc.add_to_sequence(seq, profile);

        for (int i = 0; i < seq.length(); ++i) {
            Profile<AminoAcid> p(profile, i, 1);
            hmm_.add_profile(p);
        }
        hmm_.init_transitions(HomogeneousTransitionInitializer<AminoAcid>());
        hmm_.transform_states_to_logspace();

        std::ifstream fin("../data/zinc_finger_proteins.fas");
        seqs_ = Sequence<AminoAcid>::readall(fin);
        fin.close();
    }

    HMM<AminoAcid> hmm_;
    std::vector< shared_ptr<Sequence<AminoAcid> > > seqs_;
};

TEST_F(BaumWelchTrainingTest, ZincFingersTraining)
{
    BaumWelchTraining<AminoAcid, Sequence> bwt;
    bwt.run(hmm_, seqs_);
}

}  // cs
