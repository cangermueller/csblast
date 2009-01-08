/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_profile.h"

namespace cs
{

SequenceProfile::SequenceProfile(int ncols,
                                 const SequenceAlphabet* alphabet)
        : RowMajorMatrix<float>(ncols, alphabet->size()), alphabet_(alphabet)
{}

SequenceProfile::SequenceProfile(const Sequence& sequence,
                                 const SequenceAlphabet* alphabet)
        : RowMajorMatrix<float>(sequence.length(), alphabet->size()), alphabet_(alphabet)
{
    const int rows = ncols();
    const int cols = nalph();
    for(int i=0; i<rows; ++i) {
        for(int j=0; j<cols; ++j) (*this)(i,j) = 0.0f;
        (*this)(i, sequence(i));
    }
}

SequenceProfile::~SequenceProfile()
{}

// TODO: implement input operator >> for class SequenceProfile
std::istream& operator>> (std::istream& i, SequenceProfile& profile)
{ return i; }

// TODO: implement output operator << for class SequenceProfile
std::ostream& operator<< (std::ostream& o, const SequenceProfile& profile)
{ return o; }

void reset(SequenceProfile& profile)
{
    const int ncols = profile.ncols();
    const int nalph = profile.nalph();
    for(int i=0; i<ncols; ++i)
        for(int j=0; i<nalph; ++j)
            profile(i,j) = 0.0f;
}

}//cs
