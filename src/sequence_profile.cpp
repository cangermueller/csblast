/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_profile.h"

#include <cctype>
#include <cmath>

#include <iostream>

#include "my_exception.h"
#include "profile.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "smart_ptr.h"
#include "util.h"

namespace cs
{

SequenceProfile::SequenceProfile(int ncols,
                                 const SequenceAlphabet* alphabet)
        : ncols_(ncols),
          ndim_(alphabet->size()-1),
          data_(ncols * alphabet->size()-1),
          alphabet_(alphabet)
{ reset(*this); }

SequenceProfile::SequenceProfile(std::istream& in, const SequenceAlphabet* alphabet)
            : ncols_(0),
              ndim_(alphabet->size()-1),
              alphabet_(alphabet)
{ in >> *this; }

SequenceProfile::SequenceProfile(const Sequence& sequence)
        : ncols_(sequence.length()),
          ndim_(sequence.alphabet()->size()-1),
          data_(sequence.length() * sequence.alphabet()->size()-1),
          alphabet_(sequence.alphabet())
{
    for(int i = 0; i < ncols_; ++i) {
        for(int j = 0; j < ndim_; ++j) (*this)(i,j) = 0.0f;
        (*this)(i, sequence(i));
    }
}

SequenceProfile::SequenceProfile(const SequenceProfile& other,
                                 int index,
                                 int length)
        : ncols_(length),
          ndim_(other.ndim_),
          data_(length * other.ndim_),
          alphabet_(other.alphabet_)
{
    if (index + length > other.ncols_)
        throw MyException("Arguments index=%i and length=%i for construction of sub-profile are out of bounds!", index, length);
    for(int i = 0; i < ncols_; ++i)
        for(int j = 0; j < ndim_; ++j)
            (*this)(i,j) = other(i+index,j);
}

std::vector< SmartPtr<SequenceProfile> > SequenceProfile::read(std::istream& in,
                                                               const SequenceAlphabet* alphabet)
{
    std::vector< SmartPtr<SequenceProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        SmartPtr<SequenceProfile> p(new SequenceProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void SequenceProfile::init(std::istream& in)
{
    const int kBufferSize = 1048576;  //1MB
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    //read column data records line by line
    std::vector< std::vector<int> > data;
    std::vector<int> column;
    int i_prev = 0;  //index of previous column
    while (in.getline(buffer, kBufferSize)) {
        int len = strlen(buffer);
        if (!strscn(buffer) || len > 0 && buffer[0] == '#') continue;
        if (len > 1 && buffer[0] == '/' && buffer[1] == '/') break;

        const char* ptr = buffer;
        std::vector<int> column;
        column.reserve(alphabet_->size()-1);
        int i = strtoi(ptr);
        if(i != i_prev+1)
            throw MyException("Bad format: column record %i is followed by column record %i!", i_prev, i);
        while ((ptr = strscn(ptr)))  //read column values one by one
            column.push_back(strtoi_asterix(ptr));
        if(static_cast<int>(column.size()) != alphabet_->size()-1)
            throw MyException("Bad format: column record contains %i values but should have %i!",
                              static_cast<int>(column.size()), alphabet_->size()-1);
        data.push_back(column);
        std::vector<int>().swap(column);  //clear for next use
        i_prev = i;
    }

    //transform data to lin space anf fill profile matrix
    resize(data.size(), alphabet_->size()-1);
    for (int i = 0; i < ncols_; ++i)
        for (int j = 0; j < ndim_; ++j)
            (*this)(i,j) = pow(2.0, static_cast<float>(-data[i][j]) / kScaleFactor);
}

void SequenceProfile::resize(int ncols, int ndim)
{
    if (ncols == 0 || ndim == 0)
        throw MyException("Bad dimensions for profile resizing: ncols=%i ndim=%i", ncols, ndim);
    ncols_ = ncols;
    ndim_  = ndim;
    data_.resize(ncols * ndim);
}

std::istream& operator>> (std::istream& in, SequenceProfile& profile)
{
    profile.init(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const SequenceProfile& profile)
{
    const int ncols = profile.ncols();
    const int ndim  = profile.ndim();

    // print header with character alphabet
    out << "#";
    for (int j = 0; j < ndim; ++j)
        out << "\t" << profile.alphabet()->itoc(j);
    out << std::endl;
    // print profile values in log representation
    for (int i = 0; i < ncols; ++i) {
        out << i+1;
        for (int j = 0; j < ndim; ++j) {
            double logval = log2(profile(i,j));
            if (-logval == std::numeric_limits<double>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * profile.kScaleFactor);
        }
        out << std::endl;
    }
    out << "//\n";

    return out;
}

void reset(SequenceProfile& profile, float value)
{
    const int ncols = profile.ncols();
    const int ndim = profile.ndim();
    for(int i=0; i<ncols; ++i)
        for(int j=0; i<ndim; ++j)
            profile(i,j) = value;
}

}//cs


