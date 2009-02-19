#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "alignment.h"
#include "amino_acid.h"
#include "baum_welch_training.h"
#include "blosum_matrix.h"
#include "counts_profile.h"
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
        // Build generic zinc finger HMM
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

        // Read zinc finger sequences
        std::ifstream seq_in("../data/zinc_finger_proteins.fas");
        seqs_ = Sequence<AminoAcid>::readall(seq_in);
        seq_in.close();

        // Read zinc finger alignments and construct count profiles
        std::ifstream ali_in("../data/zinc_finger_alignments.fas");
        std::vector< shared_ptr< Alignment<AminoAcid> > > alis = Alignment<AminoAcid>::readall(ali_in, Alignment<AminoAcid>::FASTA);
        ali_in.close();

        // Convert alignments to counts and add pseudocounts
        DivergenceDependentAdmixture pca2(1.0f, 5.0f);
        mpc.set_admixture(&pca2);
        for (std::vector< shared_ptr< Alignment<AminoAcid> > >::iterator ai = alis.begin(); ai != alis.end(); ++ai) {
            shared_ptr< CountsProfile<AminoAcid> > p_ptr(new CountsProfile<AminoAcid>(**ai, true));
            mpc.add_to_profile(*p_ptr);
            p_ptr->convert_to_counts();
            counts_.push_back(p_ptr);
        }
    }

    HMM<AminoAcid> hmm_;
    std::vector< shared_ptr< Sequence<AminoAcid> > > seqs_;
    std::vector< shared_ptr< CountsProfile<AminoAcid> > > counts_;
};

TEST_F(BaumWelchTrainingTest, ZincFingerSeqsTraining)
{
    BaumWelchTraining<AminoAcid, Sequence> bwt;
    bwt.run(hmm_, seqs_);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisTraining)
{
    BaumWelchTraining<AminoAcid, CountsProfile> bwt;
    bwt.run(hmm_, counts_);
}

}  // cs
