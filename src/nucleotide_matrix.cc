/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "nucleotide_matrix.h"

#include <iostream>

#include "nucleotide_alphabet.h"
#include "matrix.h"
#include "substitution_matrix.h"
#include "util.h"

namespace cs
{

NucleotideMatrix::NucleotideMatrix(float match, float mismatch)
        : SubstitutionMatrix(NucleotideAlphabet::instance()),
          match_(match),
          mismatch_(mismatch)
{
    init();
}

void NucleotideMatrix::init()
{
    // Set background vector
    f_.assign(size_, 1.0f / size_);

    // Fill substitution matrix
    for (int a = 0; a < size_; ++a)
        for (int b = 0; b < size_; ++b)
            s_[a][b] = a==b ? match_ : mismatch_;

    init_from_substitution_matrix_and_background_frequencies();
}

}  // cs
