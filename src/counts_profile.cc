/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "counts_profile.h"

#include <cmath>

#include <iomanip>
#include <iostream>

#include "alignment.h"
#include "exception.h"
#include "profile.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

CountsProfile::CountsProfile(std::istream& in, const SequenceAlphabet* alphabet)
        : Profile(alphabet),
          has_counts_(false)
{ read(in); }

CountsProfile::CountsProfile(const Sequence& sequence)
        : Profile(sequence.length(), sequence.alphabet()),
          has_counts_(false)
{
    for(int i = 0; i < ncols(); ++i)
        data_[i][sequence[i]] = 1.0f;
}

CountsProfile::CountsProfile(const Alignment& alignment, bool position_specific_weights)
        : Profile(alignment.nmatch(), alignment.alphabet()),
          has_counts_(false)
{
    const int ncols = alignment.nmatch();
    const int nseqs = alignment.nseqs();
    const int any   = alphabet()->any();

    if (position_specific_weights) {
        std::pair< matrix<float>, std::vector<float> > wi_neff = position_specific_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), wi_neff.second.begin(), wi_neff.second.end());
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += wi_neff.first[i][k];
    } else {
        std::pair<std::vector<float>, float> wg_neff = global_weights_and_diversity(alignment);
        neff_.insert(neff_.begin(), ncols, wg_neff.second);
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += wg_neff.first[k];
    }

    normalize(*this);
}

CountsProfile::CountsProfile(const CountsProfile& other,
                                   int index,
                                   int length)
        : Profile(other, index, length)
{
    neff_.insert(neff_.begin(), other.neff_.begin() + index, other.neff_.begin() + index + length);
    has_counts_ = other.has_counts_;
}

std::vector< shared_ptr<CountsProfile> > CountsProfile::readall(std::istream& in,
                                                                 const SequenceAlphabet* alphabet)
{
    std::vector< shared_ptr<CountsProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<CountsProfile> p(new CountsProfile(in, alphabet));
        profiles.push_back(p);
    }

    return profiles;
}

void CountsProfile::convert_to_counts()
{
    if (!has_counts_) {
        const bool islog = logspace();
        if (islog) transform_to_linspace();

        for (int i = 0; i < ncols(); ++i)
            for (int a = 0; a < nalph(); ++a)
                data_[i][a] *= neff_[i];
        has_counts_ = true;

        if (islog) transform_to_logspace();
    }
}

void CountsProfile::convert_to_frequencies()
{
    if (has_counts_) {
        normalize(*this);
        has_counts_ = false;
    }
}

void CountsProfile::read_header(std::istream& in)
{
    Profile::read_header(in);

    const int kBufferSize = 1000;
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    // Read has_counts
    if (in.getline(buffer, kBufferSize) && strncmp(buffer, "has_counts", 10) == 0)
        has_counts_ = atoi(buffer + 10) == 1;
}

void CountsProfile::read_body(std::istream& in)
{
    const int kBufferSize = 1000;
    std::vector<char> char_arr(kBufferSize, '\0');
    char* buffer = &char_arr[0];

    // Read column data records line by line
    neff_.resize(ncols());
    in.getline(buffer, kBufferSize);  // Skip description line
    const char* ptr;  // for string traversal
    int i = 0;        // column index
    while (in.getline(buffer, kBufferSize)) {
        if (!strscn(buffer)) continue;
        if (strlen(buffer) > 1 && buffer[0] == '/' && buffer[1] == '/') break;

        ptr = buffer;
        i = strtoi(ptr) - 1;
        if (!ptr)
            throw Exception("Bad format: malformed line after column record %i!", i - 1);
        // Read profile frequencies
        for (int a = 0; a < nalph(); ++a) {
            int log_p = strtoi_asterix(ptr);
            if (logspace_)
                data_[i][a] = static_cast<float>(-log_p) / kScaleFactor;
            else
                data_[i][a] = pow(2.0, static_cast<float>(-log_p) / kScaleFactor);
        }
        // Read neff
        neff_[i] = static_cast<float>(strtoi_asterix(ptr)) / kScaleFactor;
    }
    if (i != ncols() - 1)
        throw Exception("Bad format: alignment profile has %i column records but should have %i!", i+1, ncols());
}

void CountsProfile::write_header(std::ostream& out) const
{
    Profile::write_header(out);
    out << "has_counts\t" << has_counts_ << std::endl;
}

void CountsProfile::write_body(std::ostream& out) const
{
    // print profile values in log representation
    out << "\t";
    alphabet_->print(out, "\t");
    out << "\tneff" << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < nalph(); ++a) {
            float log_p = logspace_ ? data_[i][a] : log2(data_[i][a]);
            if (-log_p == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(log_p * kScaleFactor);
        }
        out << "\t" << iround(neff_[i] * kScaleFactor) << std::endl;
    }
    out << "//" << std::endl;
}

void CountsProfile::print(std::ostream& out) const
{
    const int kWidth = 6;
    std::ios_base::fmtflags flags = out.flags(); // save old flags

    out << std::string(2*kWidth-2, ' ');
    alphabet_->print(out, std::string(kWidth-1, ' '));
    out << std::string(kWidth-4, ' ') + "neff" << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << std::left << std::setw(kWidth-1) << i+1;
        for (int a = 0; a < nalph(); ++a)
            out << std::right << std::setw(kWidth) << std::fixed << std::setprecision(2)
                << (logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
        // print neff
        out << std::right << std::setw(kWidth-1) << std::setprecision(1) << neff_[i] << std::endl;
    }

    out.flags(flags);
}

}  // cs


