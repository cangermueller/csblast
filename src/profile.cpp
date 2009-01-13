/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "profile.h"

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

Profile::Profile(int ncols,
                 const SequenceAlphabet* alphabet)
        : ncols_(ncols),
          ndim_(alphabet->size()-1),
          data_(ncols * alphabet->size()-1),
          alphabet_(alphabet)
{ reset(*this); }

Profile::Profile(std::istream& in, const SequenceAlphabet* alphabet)
            : ncols_(0),
              ndim_(alphabet->size()-1),
              data_(),
              alphabet_(alphabet)
{ init(in); }

Profile::Profile(const Profile& other,
                 int index,
                 int length)
        : ncols_(length),
          ndim_(other.ndim_),
          data_(length * other.ndim_),
          alphabet_(other.alphabet_)
{
    if (index + length > other.ncols_)
        throw MyException("Arguments index=%i and length=%i for construction of sub-profile are out of bounds!", index, length);
    for (int i = 0; i < ncols_; ++i)
        for (int j = 0; j < ndim_; ++j)
            (*this)(i,j) = other(i+index,j);
}

std::vector< SmartPtr<Profile> > Profile::read(std::istream& in,
                                               const SequenceAlphabet* alphabet)
{
    std::vector< SmartPtr<Profile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        SmartPtr<Profile> p(new Profile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void Profile::init(std::istream& in)
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

void Profile::print(std::ostream& out) const
{
    // print header with character alphabet
    out << "#";
    for (int j = 0; j < ndim_; ++j)
        out << "\t" << alphabet_->itoc(j);
    out << std::endl;
    // print profile values in log representation
    for (int i = 0; i < ncols_; ++i) {
        out << i+1;
        for (int j = 0; j < ndim_; ++j) {
            double logval = log2((*this)(i,j));
            if (-logval == std::numeric_limits<double>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * kScaleFactor);
        }
        out << std::endl;
    }
    out << "//\n";
}

void Profile::resize(int ncols, int ndim)
{
    if (ncols == 0 || ndim == 0)
        throw MyException("Bad dimensions for profile resizing: ncols=%i ndim=%i", ncols, ndim);
    ncols_ = ncols;
    ndim_  = ndim;
    data_.resize(ncols * ndim);
}

std::istream& operator>> (std::istream& in, Profile& profile)
{
    profile.init(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const Profile& profile)
{
    profile.print(out);
    return out;
}

void reset(Profile& profile, float value)
{
    const int ncols = profile.ncols();
    const int ndim = profile.ndim();
    for(int i=0; i<ncols; ++i)
        for(int j=0; i<ndim; ++j)
            profile(i,j) = value;
}

void normalize(Profile& profile, float value)
{
    const int ncols = profile.ncols();
    const int ndim  = profile.ndim();

    for (int i = 0; i < ncols; ++i) {
        float sum = 0.0f;
        for (int a = 0; a < ndim; ++a) sum += profile(i,a);
        if (sum != 0.0f) {
            float fac = value / sum;
            for (int a = 0; a < ndim; ++a) profile(i,a) *= fac;
        } else {
            throw MyException("Unable to normalize profile to one. Sum of column %i is zero!", i);
        }
    }
}

}//cs


