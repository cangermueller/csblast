/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_profile.h"

namespace cs
{

SequenceProfile::SequenceProfile(int ncols, int ndim)
        : Profile(ncols, ndim),
          alphabet_(0)
{}

SequenceProfile::SequenceProfile(int ncols,
                                 const SequenceAlphabet* alphabet)
        : Profile(ncols, alphabet->size()),
          alphabet_(alphabet)
{}

SequenceProfile::SequenceProfile(const Sequence& sequence,
                                 const SequenceAlphabet* alphabet)
        : Profile(sequence.length(), alphabet->size()),
          alphabet_(alphabet)
{
    const int cols = ncols();
    const int dim  = ndim();
    for(int i=0; i<cols; ++i) {
        for(int j=0; j<dim; ++j) (*this)(i,j) = 0.0f;
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

}//cs
