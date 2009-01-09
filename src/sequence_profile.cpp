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
        : Profile(ncols, alphabet->size()-1),
          alphabet_(alphabet)
{}

SequenceProfile::SequenceProfile(const Sequence& sequence)
        : Profile(sequence.length(), sequence.alphabet().size()-1),
          alphabet_(&sequence.alphabet())
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

std::istream& operator>> (std::istream& in, SequenceProfile& profile)
{
    std::string data((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());

    return in;
}

// TODO: implement output operator << for class SequenceProfile
std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile)
{ return out; }

}//cs
