/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "alignment_profile.h"

#include <cmath>

#include <iostream>

#include "alignment.h"
#include "my_exception.h"
#include "profile.h"
#include "sequence_alphabet.h"
#include "smart_ptr.h"
#include "util.h"

namespace cs
{

AlignmentProfile::AlignmentProfile(std::istream& in, const SequenceAlphabet* alphabet)
        : Profile(alphabet),
          has_counts_(false)
{ unserialize(in); }

AlignmentProfile::AlignmentProfile(const Sequence& sequence)
        : Profile(sequence.length(), sequence.alphabet()),
          has_counts_(false)
{
    for(int i = 0; i < ncols(); ++i) {
        for(int j = 0; j < ndim(); ++j) (*this)(i,j) = 0.0f;
        (*this)(i, sequence(i));
    }
}

AlignmentProfile::AlignmentProfile(const Alignment& alignment, bool position_dependent_weights)
        : Profile(alignment.ncols(), alignment.alphabet()),
          has_counts_(false)
{
    const int ncols = alignment.ncols();
    const int nseqs = alignment.nseqs();
    const int any   = alphabet()->any();

    if (position_dependent_weights) {
        std::pair< Matrix<float>, std::vector<float> > wi_neff = cs::position_dependent_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), wi_neff.second.begin(), wi_neff.second.end());
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment(k,i) < any)
                    (*this)(i, alignment(k,i)) += wi_neff.first(i,k);
    } else {
        std::pair<std::vector<float>, float> wg_neff = cs::global_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), ncols, wg_neff.second);
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment(k,i) < any)
                    (*this)(i, alignment(k,i)) += wg_neff.first[k];
    }

    normalize(*this);
}

AlignmentProfile::AlignmentProfile(const AlignmentProfile& other,
                                   int index,
                                   int length)
        : Profile(other, index, length)
{
    neff_.insert(neff_.begin(), other.neff_.begin() + index, other.neff_.begin() + index + length);
    has_counts_ = other.has_counts_;
}

std::vector< SmartPtr<AlignmentProfile> > AlignmentProfile::read(std::istream& in,
                                                                 const SequenceAlphabet* alphabet)
{
    std::vector< SmartPtr<AlignmentProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        SmartPtr<AlignmentProfile> p(new AlignmentProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void AlignmentProfile::convert_to_counts()
{
    if (!has_counts_) {
        for (int i = 0; i < ncols(); ++i)
            for (int a = 0; a < ndim(); ++a)
                (*this)(i,a) *= neff_[i];
        has_counts_ = true;
    }
}

void AlignmentProfile::convert_to_frequencies()
{
    if (has_counts_) {
        normalize(*this);
        has_counts_ = false;
    }
}

void AlignmentProfile::unserialize(std::istream& in)
{
    const int kBufferSize = 1000;
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    // Check if stream actually contains a serialized alignment profile
    while (in.getline(buffer, kBufferSize) && !strscn(buffer)) continue;
    if (strcmp(buffer, "AlignmentProfile") != 0)
        throw MyException("Bad format: serialized profile does not start with 'AlignmentProfile' class identifier!");

    Profile::unserialize(in);  // Initialize base class part

    // Read has_counts
    if (in.getline(buffer, kBufferSize) && strncmp(buffer, "has_counts", 10) == 0)
        has_counts_ = atoi(buffer+10) == 1;

    // read column data records line by line
    neff_.resize(ncols());
    in.getline(buffer, kBufferSize);  // Skip description line
    const char* ptr;  // for string traversal
    int i = 0;        // column index
    while (in.getline(buffer, kBufferSize)) {
        if (!strscn(buffer)) continue;
        if (strlen(buffer) > 1 && buffer[0] == '/' && buffer[1] == '/' || buffer[0] == '#') break;

        ptr = buffer;
        i = strtoi(ptr)-1;
        if (!ptr)
            throw MyException("Bad format: malformed line after neff record %i!", i-1);

        int logval = strtoi_asterix(ptr);
        neff_[i] = pow(2.0, static_cast<float>(-logval) / kScaleFactor);
    }
    if (i != ncols()-1)
        throw MyException("Bad format: neff vector has %i column records but should have %i!", i+1, ncols());
}

void AlignmentProfile::serialize(std::ostream& out) const
{
    out << "AlignmentProfile" << std::endl;
    Profile::serialize(out);  // Print the profile part

    // print has_counts
    out << "has_counts\t" << (has_counts_ ? 1 : 0) << std::endl;
    // print neff vector
    out << "\tneff\n";
    // print number of effective sequences in log representation
    for (int i = 0; i < ncols(); ++i) {
        double logval = log2(neff_[i]);
        if (-logval == std::numeric_limits<double>::infinity())
            out << i+1 << "\t*" << std::endl;
        else
            out << i+1 << "\t" << -iround(logval * kScaleFactor) << std::endl;
    }
    out << "//" << std::endl;;
}

}//cs


