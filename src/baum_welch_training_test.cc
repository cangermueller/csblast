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
        seqs_ = Sequence<AminoAcid>::readall(seq_in);
        seq_in.close();

        // Read zinc finger alignments and construct count profiles
        std::ifstream ali_in("../data/zinc_finger_alignments.fas");
        std::vector< shared_ptr< Alignment<AminoAcid> > > alis = Alignment<AminoAcid>::readall(ali_in, Alignment<AminoAcid>::FASTA);
        ali_in.close();

        // Convert alignments to counts and add pseudocounts
        for (std::vector< shared_ptr< Alignment<AminoAcid> > >::iterator ai = alis.begin(); ai != alis.end(); ++ai) {
            shared_ptr< CountsProfile<AminoAcid> > p_ptr(new CountsProfile<AminoAcid>(**ai, true));
            mpc.add_to_profile(*p_ptr, DivergenceDependentAdmixture(0.5f, 5.0f));
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
    TrainingProgressInfo<AminoAcid> prg_info(hmm_, std::cout);
    BaumWelchTraining<AminoAcid, Sequence> bwt;
    bwt.run(hmm_, seqs_, &prg_info);

    hmm_.transform_states_to_linspace();
    EXPECT_NEAR(1.0, hmm_[0][0][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(1.0, hmm_[5][0][AminoAcid::instance().ctoi('C')], DELTA);
}

TEST_F(BaumWelchTrainingTest, ZincFingerAlisTraining)
{
    TrainingProgressInfo<AminoAcid> prg_info(hmm_, std::cout);
    BaumWelchTraining<AminoAcid, CountsProfile> bwt;
    bwt.run(hmm_, counts_, &prg_info);

    hmm_.transform_states_to_linspace();
    EXPECT_NEAR(0.7424, hmm_[0][0][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(0.7424, hmm_[3][0][AminoAcid::instance().ctoi('C')], DELTA);
}

}  // cs
