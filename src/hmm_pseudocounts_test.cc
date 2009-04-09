#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "amino_acid.h"
#include "count_profile-inl.h"
#include "emitter.h"
#include "hmm-inl.h"
#include "hmm_pseudocounts.h"
#include "library_pseudocounts.h"
#include "log.h"

namespace cs
{

const float DELTA = 0.01f;

TEST(HMMPseudocountsTest, AddToSequence)
{
    Sequence<AminoAcid> seq("triple zinc finger",
                            "KPSRMRKYPNRPSKTPPHERPYACPVESCDRRFSRSDELTRHIRIHTGQKPFQCRICMRNFSRSDHLTTH");
    Profile<AminoAcid> profile(seq.length());

    std::ifstream in("../data/scop20_K100.hmm");
    HMM<AminoAcid> hmm(in);
    ASSERT_EQ(100, hmm.num_states());

    EmissionParams params;
    HMMPseudocounts<AminoAcid> pc(&hmm, params);
    pc.add_to_sequence(seq, DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

    EXPECT_NEAR(0.8121f, profile[23][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(0.8121f, profile[28][AminoAcid::instance().ctoi('C')], DELTA);
}

TEST(HMMPseudocountsTest, AddProfileSequence)
{
    std::ifstream ali_in("../data/zinc_finger.fas");
    Alignment<AminoAcid> ali(ali_in, Alignment<AminoAcid>::FASTA);
    ali_in.close();
    CountProfile<AminoAcid> profile(ali, false);

    std::ifstream in("../data/scop20_K100.hmm");
    HMM<AminoAcid> hmm(in);
    ASSERT_EQ(100, hmm.num_states());

    EmissionParams params;
    HMMPseudocounts<AminoAcid> pc(&hmm, params);
    pc.add_to_profile(DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

    EXPECT_NEAR(0.7756f, profile[53][AminoAcid::instance().ctoi('C')], DELTA);
    EXPECT_NEAR(0.7720f, profile[56][AminoAcid::instance().ctoi('C')], DELTA);
}

}  // namespace cs


