/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alphabet.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

#include <vector>

namespace cs
{

void SequenceAlphabet::init()
{
    std::vector<char>().swap(itoc_);  // clear and minimize capazity

    const char* itoc = get_itoc();  // derived classes decide how to fil itoc array
    itoc_.insert(itoc_.begin(), itoc, itoc + strlen(itoc));

    const int ctoi_size = static_cast<int>( pow(2, 8*sizeof(char)) );
    const int itoc_size = itoc_.size();
    ctoi_ = std::vector<int>(ctoi_size, kInvalidChar);
    for (int i = 0; i < itoc_size; ++i) ctoi_[static_cast<int>(itoc_[i])]=i;
}

}//cs


