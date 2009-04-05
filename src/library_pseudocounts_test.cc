#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "amino_acid.h"
#include "emitter.h"
#include "log.h"
#include "profile_library.h"
#include "library_pseudocounts.h"

namespace cs
{

const float DELTA = 0.01f;

TEST(LibraryPseudocountsTest, AddToSequence)
{
    const Sequence<AminoAcid> seq("header", "ARNDCQEGHILKMFPSTWYV");
    Profile<AminoAcid> profile(seq.length());

    ASSERT_EQ(AminoAcid::instance().size(), seq.length());
    ASSERT_EQ(AminoAcid::instance().ctoi('R'), seq[1]);
    ASSERT_EQ(seq.length(), profile.num_cols());

    std::ifstream in("../data/scop20_1.73_opt_N100000_W13.lib");
    ProfileLibrary<AminoAcid> lib(in);
    ASSERT_EQ(50, lib.num_profiles());

    EmitterParams params;
    Emitter<AminoAcid> emitter(lib.num_cols(), params);

    LibraryPseudocounts<AminoAcid> pc(&lib, &emitter);
    pc.add_to_sequence(seq, DivergenceDependentAdmixture(1.0f, 10.0f), &profile);

    EXPECT_NEAR(0.0736f, profile[0][AminoAcid::instance().ctoi('V')], DELTA);
}

}  // cs
