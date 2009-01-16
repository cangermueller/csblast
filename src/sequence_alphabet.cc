/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alphabet.h"

#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

#include <algorithm>
#include <vector>

#include "my_exception.h"

namespace cs
{

SequenceAlphabet::SequenceAlphabet(int size)
        : size_(size),
          ctoi_(static_cast<int>(pow(2, 8*sizeof(char))), kInvalidChar),
          itoc_(size, '\0')
{}

void SequenceAlphabet::init()
{
    const char* itoc = get_itoc();  // derived classes decide how to fil itoc array
    if (static_cast<int>(strlen(itoc)) != size_)
        throw MyException("Provided itoc array has length %i but should have length %i!", strlen(itoc), size_);
    copy(itoc, itoc + strlen(itoc), itoc_.begin());
    for (int i = 0; i < size_; ++i) ctoi_[static_cast<int>(itoc_[i])]=i;
}

}//cs


