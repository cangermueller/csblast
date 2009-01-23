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
#include <string>
#include <vector>

#include "exception.h"

namespace cs
{

SequenceAlphabet::SequenceAlphabet(int size, char any, char gap)
        : size_(size),
          any_(any),
          gap_(gap),
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
    itoc_[size_]   = any_;    // ANY
    itoc_[size_ + 1] = gap_;  // GAP
    itoc_[size_ + 2] = gap_;  // ENDGAP

    // Setup ctoi vector
    for (int i = 0; i < size_; ++i) ctoi_[static_cast<int>(itoc_[i])] = i;
    ctoi_[static_cast<int>(any_)] = size_;      // ANY
    ctoi_[static_cast<int>(gap_)] = size_ + 1;  // GAP
}

void SequenceAlphabet::print(std::ostream& out, const std::string& delim) const
{
    out << *begin();
    for (const_iterator iter = begin() + 1; iter != end(); ++iter)
        out << delim << *iter;
}

}//cs


