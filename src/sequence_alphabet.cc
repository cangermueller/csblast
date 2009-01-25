/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alphabet.h"

#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>

#include <algorithm>
#include <string>
#include <vector>

#include "exception.h"

namespace cs
{

SequenceAlphabet::SequenceAlphabet(int size, char any)
        : size_(size),
          any_(any),
          ctoi_(static_cast<int>(pow(2, 8*sizeof(char))), kInvalidChar),
          itoc_(size + 3, '\0')
{}

void SequenceAlphabet::init()
{
    const char* itoc = get_itoc();  // derived classes decide how to fill itoc array
    if (static_cast<int>(strlen(itoc)) != size_)
        throw Exception("Provided itoc array has length %i but should have length %i!", strlen(itoc), size_);

    // Setup itoc vector
    copy(itoc, itoc + strlen(itoc), itoc_.begin());
    itoc_[size_]   = any_;   // ANY
    itoc_[size_ + 1] = '-';  // GAP
    itoc_[size_ + 2] = '-';  // ENDGAP

    // Setup ctoi vector
    for (int i = 0; i < size_; ++i) ctoi_[toupper(itoc_[i])] = i;
    ctoi_[toupper(any_)] = size_;      // ANY
    ctoi_['-']           = size_ + 1;  // MATCH GAP
    ctoi_['.']           = size_ + 1;  // INSERT GAP
}

void SequenceAlphabet::print(std::ostream& out, const std::string& delim) const
{
    out << *begin();
    for (const_iterator iter = begin() + 1; iter != end(); ++iter)
        out << delim << *iter;
}

}//cs


