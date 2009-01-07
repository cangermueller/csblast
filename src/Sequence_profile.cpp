/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "Sequence_profile.h"


Sequence_profile::Sequence_profile(size_t ncols, const Sequence_alphabet* const alphabet)
    : Row_major_matrix<float>(ncols, alphabet->size()), alphabet_(alphabet)
{ }


Sequence_profile::~Sequence_profile()
{ }


const Sequence_alphabet& Sequence_profile::alphabet()
{ return *alphabet_; }

