/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "alignment_profile.h"

#include <cmath>

#include <iostream>

#include "alignment.h"
#include "my_exception.h"
#include "profile.h"
#include "sequence_alphabet.h"
#include "smart_ptr.h"

namespace cs
{

AlignmentProfile::AlignmentProfile(std::istream& in, const SequenceAlphabet* alphabet)
        : Profile(0, alphabet)  // dont't call input stream constructor
{ init(in); }

AlignmentProfile::AlignmentProfile(const Sequence& sequence)
        : Profile(sequence.length(), sequence.alphabet())
{
    for(int i = 0; i < ncols(); ++i) {
        for(int j = 0; j < ndim(); ++j) (*this)(i,j) = 0.0f;
        (*this)(i, sequence(i));
    }
}

AlignmentProfile::AlignmentProfile(const Alignment& alignment, bool position_dependent_weights)
        : Profile(alignment.ncols(), alignment.alphabet())
{
    const int ncols = alignment.ncols();
    const int nseqs = alignment.nseqs();
    const int any   = alphabet()->any();

    if (position_dependent_weights) {
        std::pair< Matrix<float>, std::vector<float> > wi_neff = cs::position_dependent_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), wi_neff.second.begin(), wi_neff.second.end());
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment(k,i) < any)
                    (*this)(i, alignment(k,i)) += wi_neff.first(i,k);
    } else {
        std::pair<std::vector<float>, float> wg_neff = cs::global_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), ncols, wg_neff.second);
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment(k,i) < any)
                    (*this)(i, alignment(k,i)) += wg_neff.first[k];
    }
    normalize(*this);
}

AlignmentProfile::AlignmentProfile(const AlignmentProfile& other,
                                   int index,
                                   int length)
        : Profile(other, index, length)
{
    neff_.insert(neff_.begin(), other.neff_.begin() + index, other.neff_.begin() + index + length);
}

std::vector< SmartPtr<AlignmentProfile> > AlignmentProfile::read(std::istream& in,
                                                                 const SequenceAlphabet* alphabet)
{
    std::vector< SmartPtr<AlignmentProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        SmartPtr<AlignmentProfile> p(new AlignmentProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void AlignmentProfile::convert_to_counts()
{
    if (!has_counts_) {
        for (int i = 0; i < ncols(); ++i)
            for (int a = 0; a < ndim(); ++a)
                (*this)(i,a) *= neff_[i];
    }
}

void AlignmentProfile::convert_to_frequencies()
{
    if (has_counts_) normalize(*this);
}

void AlignmentProfile::init(std::istream& in)
{
    Profile::init(in);  // Initialize base class part
    // TODO
}

void AlignmentProfile::print(std::ostream& out) const
{
    Profile::print(out);  // Print the base class part
    // TODO
}

}//cs


