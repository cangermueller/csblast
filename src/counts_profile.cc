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
#include "matrix.h"
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
{
    read(in);
}

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
        matrix<float> w;  // position-specific sequence weights
        neff_ = position_specific_weights_and_diversity(alignment, w);
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += w[i][k];
    } else {
        std::vector<float> wg;  // global sequence weights
        neff_.assign(ncols, global_weights_and_diversity(alignment, wg));
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += wg[k];
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
    neff_.resize(ncols());

    // Read has_counts
   std::string tmp;
   if (getline(in, tmp) && tmp.find("has_counts") != std::string::npos)
        has_counts_ = atoi(tmp.c_str() + 10) == 1;
}

void CountsProfile::read_body(std::istream& in)
{
    std::string tmp;
    std::vector<std::string> tokens;
    int i = 0;
    getline(in, tmp);  // skip alphabet description line
    while (getline(in, tmp)) {
        if (tmp.empty()) continue;
        if (tmp.length() > 1 && tmp[0] == '/' && tmp[1] == '/') break;

        split(tmp, '\t', tokens);
        i = atoi(tokens[0].c_str()) - 1;
        for (int a = 0; a < nalph(); ++a) {
            float log_p = tokens[a+1][0] == '*' ? std::numeric_limits<int>::max() : atoi(tokens[a+1].c_str());
            data_[i][a] = (logspace_ ? -log_p / kScaleFactor : pow(2.0, -log_p / kScaleFactor)) ;
        }
        neff_[i] = atof(tokens[nalph()+1].c_str()) / kScaleFactor;
        tokens.clear();
    }
    if (i != ncols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, ncols());
}

void CountsProfile::write_header(std::ostream& out) const
{
    Profile::write_header(out);
    out << "has_counts\t" << has_counts_ << std::endl;
}

void CountsProfile::write_body(std::ostream& out) const
{
    out << "\t" << *alphabet_ << "\tneff" << std::endl;
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
    std::ios_base::fmtflags flags = out.flags();  // save old flags

    out << "\t" << *alphabet_ << "\tneff" << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < nalph(); ++a)
            out << '\t' << std::fixed << std::setprecision(4)
                << (logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
        // print neff
        out << '\t' << std::setprecision(2) << neff_[i] << std::endl;
    }

    out.flags(flags);
}

}  // cs


