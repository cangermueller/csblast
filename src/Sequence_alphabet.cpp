/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "Sequence_alphabet.h"


void Sequence_alphabet::init()
{
    itoc_ = itoc();
    const size_t char_size = static_cast<int>(pow(2, 8*sizeof(char)));
    for(size_t i=0; i<char_size; ++i) ctoi_.push_back(-1);
    for(size_t i=0; i<itoc_.size(); ++i) ctoi_[itoc_[i]]=i;
}

