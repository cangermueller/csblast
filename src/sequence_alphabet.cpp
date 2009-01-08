/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alphabet.h"

namespace cs
{

void SequenceAlphabet::init()
{
    init_itoc();
    const int char_size = static_cast<int>(pow(2, 8*sizeof(char)));
    for(int i=0; i<char_size; ++i) ctoi_.push_back(-1);
    const int itoc_size = itoc_.size();
    for(int i=0; i<itoc_size; ++i) ctoi_[itoc_[i]]=i;
}

}//cs


