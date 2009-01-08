/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_profile.h"

namespace cs
{

SequenceProfile::SequenceProfile(size_t ncols, const SequenceAlphabet* const alphabet)
    : RowMajorMatrix<float>(ncols, alphabet->size()), alphabet_(alphabet)
{ }


SequenceProfile::~SequenceProfile()
{ }


const SequenceAlphabet& SequenceProfile::alphabet()
{ return *alphabet_; }

}//cs
