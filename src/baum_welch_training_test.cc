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
#include "utils.h"

namespace cs
{

const float DELTA = 0.01;

class BaumWelchTrainingTest : public testing::Test
{
  protected:
    BaumWelchTrainingTest()
            : hmm_(23, 1)
    { }

    virtual void SetUp()
    {
        // Build generic zinc finger HMM
        Sequence<AminoAcid> seq("zinc finger motif", "CPVESCDRRFSRSDELTRHIRIH\n");
        Profile<AminoAcid> profile(seq.length());

        BlosumMatrix m;
        MatrixPseudocounts<AminoAcid> mpc(&m);
        mpc.add_to_sequence(seq, profile, ConstantAdmixture(1.0f));

        for (int i = 0; i < seq.length(); ++i) {
            Profile<AminoAcid> p(profile, i, 1);
            hmm_.add_profile(p);
        }
        hmm_.init_transitions(HomogeneousTransitionInitializer<AminoAcid>());
        hmm_.transform_states_to_logspace();

        // Read zinc finger sequences
        std::ifstream seq_in("../data/zinc_finger_proteins.fas");
        Sequence<AminoAcid>::readall(seq_in, seqs_);
        seq_in.close();

        // Read zinc finger alignments and construct count profiles
        std::ifstream ali_in("../data/zinc_finger_alignments.fas");
        std::vector< shared_ptr< Alignment<AminoAcid> > > alis;
        Alignment<AminoAcid>::readall(ali_in, Alignment<AminoAcid>::FASTA, alis);
        ali_in.close();

        // Convert alignments to counts and add pseudocounts
        for (std::vector< shared_ptr< Alignment<AminoAcid> > >::iterator ai = alis.begin(); ai != alis.end(); ++ai) {
            shared_ptr< CountsProfile<AminoAcid> > p_ptr(new CountsProfile<AminoAcid>(**ai, true));
            mpc.add_to_profile(*p_ptr, ConstantAdmixture(0.01f));
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
    BaumWelchParams p;
    p.num_blocks   = 1;
    p.epsilon_null = 1.0;
    p.beta         = 0.0;
    TrainingProgressInfo prg_info(std::cout);
    BaumWelchTraining<AminoAcid, Sequence> bwt(p, seqs_, hmm_, &prg_info);
    bwt.run();

    hmm_.transform_states_to_linspace();
    EXPECT_NEAR(0.8977, hmm_[0][0][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(0.8977, hmm_[5][0][AminoAcid::instance().ctoi('C')], DELTA);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisTraining)
{
    BaumWelchParams p;
    p.transition_pseudocounts  = 0.8;
    p.log_likelihood_change    = 0.02;
    p.max_connectivity         = 5;
    p.num_blocks               = 1;
    p.epsilon_null             = 1.0;
    p.beta                     = 0.0;
    TrainingProgressInfo prg_info(std::cout);
    BaumWelchTraining<AminoAcid, CountsProfile> bwt(p, counts_, hmm_, &prg_info);
    bwt.run();

    hmm_.transform_states_to_linspace();
    EXPECT_NEAR(0.9948, hmm_[0][0][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(0.9948, hmm_[5][0][AminoAcid::instance().ctoi('C')], DELTA);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisOnlineTraining)
{
    BaumWelchParams p;
    p.transition_pseudocounts  = 0.8;
    p.log_likelihood_change    = 0.02;
    p.max_connectivity         = 5;
    p.num_blocks               = 2;
    p.epsilon_null             = 0.05;
    p.beta                     = 0.5;
    TrainingProgressInfo prg_info(std::cout);
    BaumWelchTraining<AminoAcid, CountsProfile> bwt(p, counts_, hmm_, &prg_info);
    bwt.run();

    hmm_.transform_states_to_linspace();
    // TODO: write assertions!!!
}

}  // cs
