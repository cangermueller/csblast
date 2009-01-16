/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "substitution_matrix.h"

namespace cs
{

SubstitutionMatrix::SubstitutionMatrix(int size)
  : size_(size),
    p_(size, size, 0.0f),
    s_(size, size, 0.0f),
    r_(size, size, 0.0f),
    f_(size, 0.0f)
{}

void SubstitutionMatrix::init_from_target_frequencies()
{
    // // Check transition probability matrix, renormalize p and calculate f[a]
    // float sumab = 0.0f;
    // for (int a = 0; a <20; a++)
    //     for (b=0; b<20; ++b) sumab+=P[a][b];
    // for (a=0; a<20; a++)
    //     for (b=0; b<20; ++b) P[a][b]/=sumab;
    // for (a=0; a<20; a++)
    // for (pb[a]=0.0f, b=0; b<20; ++b) pb[a]+=P[a][b];

    // //Compute similarity matrix for amino acid pairs (for calculating consensus sequence)
    // for (a=0; a<20; ++a)
    //     for (b=0; b<20; ++b)
    //   Sim[a][b] = P[a][b]*P[a][b]/P[a][a]/P[b][b];

    // //Precompute matrix R for amino acid pseudocounts:
    // for (a=0; a<20; ++a)
    //     for (b=0; b<20; ++b)
    //         R[a][b] = P[a][b]/pb[b]; //R[a][b]=P(a|b)

    // //Precompute matrix R for amino acid pseudocounts:
    // for (a=0; a<20; ++a)
    //     for (b=0; b<20; ++b)
    //         S[a][b] = log2(R[a][b]/pb[a]); // S[a][b] = log2(P(a,b)/P(a)/P(b))

}

}  // cs
