/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alphabet.h"

#include <cmath>
#include <cstddef>

#include <vector>

namespace cs
{

SequenceAlphabet::SequenceAlphabet()
{}

SequenceAlphabet::~SequenceAlphabet()
{}

void SequenceAlphabet::init()
{
    std::vector<char>().swap(itoc_); // clear and minimize capazity
    init_itoc(); //subclasses define their alphabet and its integer representation

    const int ctoi_size = static_cast<int>(pow(2, 8*sizeof(char)));
    ctoi_ = std::vector<int>(ctoi_size, kInvalidChar);
    const int itoc_size = itoc_.size();
    for(int i = 0; i < itoc_size; ++i) ctoi_[static_cast<int>(itoc_[i])]=i;
}

}//cs


