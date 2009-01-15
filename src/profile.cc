/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "profile.h"

#include <cmath>
#include <cstring>

#include <iostream>

#include "my_exception.h"
#include "profile.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "smart_ptr.h"
#include "util.h"

namespace cs
{

Profile::Profile(const SequenceAlphabet* alphabet)
        : ncols_(0),
          ndim_(0),
          alphabet_(alphabet)
{}

Profile::Profile(int ncols,
                 const SequenceAlphabet* alphabet)
        : ncols_(ncols),
          ndim_(alphabet->size()-1),
          data_(ncols * alphabet->size()-1, 0.0f),
          alphabet_(alphabet)
{}

Profile::Profile(std::istream& in, const SequenceAlphabet* alphabet)
            : ncols_(0),
              ndim_(alphabet->size()-1),
              alphabet_(alphabet)
{ unserialize(in); }

Profile::Profile(const Profile& other,
                 int index,
                 int length)
        : ncols_(length),
          ndim_(other.ndim_),
          data_(length * other.ndim_, 0.0f),
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

void Profile::unserialize(std::istream& in)
{
    const int kBufferSize = 1000;
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    // Check if stream actually contains a serialized profile
    while (in.getline(buffer, kBufferSize) && !strscn(buffer)) continue;
    if (strcmp(buffer, "Profile") != 0)
        throw MyException("Bad format: serialized profile does not start with 'Profile' class identifier!");

    // Read ncols
    if (in.getline(buffer, kBufferSize) && strncmp(buffer, "ncols", 5) == 0)
        ncols_ = atoi(buffer+5);
    else
        throw MyException("Bad format: serialized profile does not contain 'ncols' record!");

    // Read ndim
    if (in.getline(buffer, kBufferSize) && strncmp(buffer, "ndim", 4) == 0)
        ndim_ = atoi(buffer+4);
    else
        throw MyException("Bad format: serialized profile does not contain 'ndim' record!");
    if (ndim_ != alphabet_->size() - 1)
        throw MyException("Bad format: ndim=%i does not fit with provided alphabet!", ndim_);

    // Read column data records line by line
    resize(ncols_, ndim_);
    in.getline(buffer, kBufferSize);  // Skip description line
    const char* ptr;  // for string traversal
    int i = 0;        // column index
    while (in.getline(buffer, kBufferSize)) {
        if (!strscn(buffer)) continue;
        if (strlen(buffer) > 1 && buffer[0] == '/' && buffer[1] == '/') break;

        ptr = buffer;
        i = strtoi(ptr)-1;
        if (!ptr)
            throw MyException("Bad format: malformed line after column record %i!", i-1);

        for (int a = 0; a < ndim_; ++a) {
            int log_p = strtoi_asterix(ptr);
            (*this)(i,a) = pow(2.0, static_cast<float>(-log_p) / kScaleFactor);
        }
    }
    if (i != ncols_-1)
        throw MyException("Bad format: profile has %i column records but should have %i!", i+1, ncols_);
}

void Profile::serialize(std::ostream& out) const
{
    out << "Profile" << std::endl;
    out << "ncols\t" << ncols_ << std::endl;
    out << "ndim\t" << ndim_ << std::endl;

    // print profile values in log representation
    for (int j = 0; j < ndim_; ++j)
        out << "\t" << alphabet_->itoc(j);
    out << std::endl;
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
    out << "//" << std::endl;
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
    profile.unserialize(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const Profile& profile)
{
    profile.serialize(out);
    return out;
}

void reset(Profile& profile, float value)
{
    const int ncols = profile.ncols();
    const int ndim = profile.ndim();
    for(int i = 0; i < ncols; ++i)
        for(int j = 0; i < ndim; ++j)
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


