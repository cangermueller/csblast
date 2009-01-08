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

SequenceProfile::~SequenceProfile()
{}

// TODO: implement input operator >> for class SequenceProfile
std::istream& operator>> (std::istream& i, SequenceProfile& profile)
{ return i; }

// TODO: implement output operator << for class SequenceProfile
std::ostream& operator<< (std::ostream& o, const SequenceProfile& profile)
{ return o; }

}//cs
