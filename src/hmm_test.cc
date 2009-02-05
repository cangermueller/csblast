#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "amino_acid.h"
#include "hmm.h"
#include "log.h"
#include "profile.h"

namespace cs
{

const float DELTA = 0.01f;

TEST(HMMTest, SimpleConstruction)
{
    HMM<AminoAcid> hmm(3);

    EXPECT_EQ(0, hmm.num_states());
    EXPECT_EQ(0, hmm.num_transitions());
    EXPECT_EQ(3, hmm.size());
}

}  // cs

