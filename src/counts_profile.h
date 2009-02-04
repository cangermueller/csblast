#ifndef CS_COUNTS_PROFILE_H
#define CS_COUNTS_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for profiles derived from alignments.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "alignment.h"
#include "exception.h"
#include "matrix.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

template<class AlphabetType>
class CountsProfile : public Profile<AlphabetType>
{
  public:
    // Constructs profile from serialized profile read from input stream.
    CountsProfile(std::istream& in);
    // Constructs a profile of the given sequence.
    explicit CountsProfile(const Sequence<AlphabetType>& sequence);
    // Constructs a profile of the given alignment with specified sequence weighting method.
    explicit CountsProfile(const Alignment<AlphabetType>& alignment, bool position_specific_weights = true);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    CountsProfile(const CountsProfile& other, int index, int length);

    virtual ~CountsProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<CountsProfile> > readall(std::istream& in);
    // Returns the number of effective sequences in alignment column i
    float neff(int i) const { return neff_[i]; }
    // Converts the profile to counts of alphabet letters.
    void convert_to_counts();
    // Converts the profile back to relative frequencies of alphabet letters.
    void convert_to_frequencies();
    // Returns true if the profile contains counts.
    bool has_counts() const { return has_counts_; }

  protected:
    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized array data members to stream.
    virtual void write_body(std::ostream& out) const;

  private:
    // Disallow copy and assign
    CountsProfile(const CountsProfile&);
    void operator=(const CountsProfile&);

    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;
    // Return serialization class identity.
    virtual const std::string class_identity() const { static std::string id("CountsProfile"); return id;}

    // Flag indicating if the profile contains counts or (relative) frequencies.
    bool has_counts_;
    // Number of effective sequences in each alignment column.
    std::vector<float> neff_;
};  // CountsProfile



template<class AlphabetType>
CountsProfile<AlphabetType>::CountsProfile(std::istream& in)
        : has_counts_(false)
{
    this->read(in);
}

template<class AlphabetType>
CountsProfile<AlphabetType>::CountsProfile(const Sequence<AlphabetType>& sequence)
        : Profile<AlphabetType>(sequence.length()),
          has_counts_(false)
{
    for(int i = 0; i < this->ncols(); ++i)
        this->data_[i][sequence[i]] = 1.0f;
}

template<class AlphabetType>
CountsProfile<AlphabetType>::CountsProfile(const Alignment<AlphabetType>& alignment, bool position_specific_weights)
        : Profile<AlphabetType>(alignment.nmatch()),
          has_counts_(false)
{
    const int ncols = alignment.nmatch();
    const int nseqs = alignment.nseqs();
    const int any   = AlphabetType::instance().any();

    if (position_specific_weights) {
        matrix<float> w;  // position-specific sequence weights
        neff_ = position_specific_weights_and_diversity(alignment, w);
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    this->data_[i][alignment[i][k]] += w[i][k];
    } else {
        std::vector<float> wg;  // global sequence weights
        neff_.assign(ncols, global_weights_and_diversity(alignment, wg));
        for (int i = 0; i < ncols; ++i)
            for (int k = 0; k < nseqs; ++k)
                if (alignment[i][k] < any)
                    this->data_[i][alignment[i][k]] += wg[k];
    }

    normalize(*this);
}

template<class AlphabetType>
CountsProfile<AlphabetType>::CountsProfile(const CountsProfile& other,
                                           int index,
                                           int length)
        : Profile<AlphabetType>(other, index, length)
{
    neff_.insert(neff_.begin(), other.neff_.begin() + index, other.neff_.begin() + index + length);
    has_counts_ = other.has_counts_;
}

template<class AlphabetType>
std::vector< shared_ptr< CountsProfile<AlphabetType> > > CountsProfile<AlphabetType>::readall(std::istream& in)
{
    std::vector< shared_ptr<CountsProfile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<CountsProfile> p(new CountsProfile(in));
        profiles.push_back(p);
    }

    return profiles;
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::convert_to_counts()
{
    if (!has_counts_) {
        const bool islog = this->logspace();
        if (islog) this->transform_to_linspace();

        for (int i = 0; i < this->ncols(); ++i)
            for (int a = 0; a < this->nalph(); ++a)
                this->data_[i][a] *= neff_[i];
        has_counts_ = true;

        if (islog) this->transform_to_logspace();
    }
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::convert_to_frequencies()
{
    if (has_counts_) {
        normalize(*this);
        has_counts_ = false;
    }
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::read_header(std::istream& in)
{
    Profile<AlphabetType>::read_header(in);
    neff_.resize(this->ncols());

    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("has_counts") != std::string::npos)
        has_counts_ = atoi(tmp.c_str() + 10) == 1;
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::read_body(std::istream& in)
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
        for (int a = 0; a < this->nalph(); ++a) {
            float log_p = tokens[a+1][0] == '*' ? std::numeric_limits<int>::max() : atoi(tokens[a+1].c_str());
            this->data_[i][a] = (this->logspace_ ? -log_p / this->SCALE_FACTOR : pow(2.0, -log_p / this->SCALE_FACTOR)) ;
        }
        neff_[i] = atof(tokens[this->nalph()+1].c_str()) / this->SCALE_FACTOR;
        tokens.clear();
    }
    if (i != this->ncols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, this->ncols());
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::write_header(std::ostream& out) const
{
    Profile<AlphabetType>::write_header(out);
    out << "has_counts\t" << has_counts_ << std::endl;
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::write_body(std::ostream& out) const
{
    out << "\t" << AlphabetType::instance() << "\tNeff" << std::endl;
    for (int i = 0; i < this->ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < this->nalph(); ++a) {
            float log_p = this->logspace_ ? this->data_[i][a] : log2(this->data_[i][a]);
            if (-log_p == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(log_p * this->SCALE_FACTOR);
        }
        out << "\t" << iround(neff_[i] * this->SCALE_FACTOR) << std::endl;
    }
    out << "//" << std::endl;
}

template<class AlphabetType>
void CountsProfile<AlphabetType>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save old flags

    out << "\t" << AlphabetType::instance() << "\tNeff" << std::endl;
    for (int i = 0; i < this->ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < this->nalph(); ++a)
            out << '\t' << std::fixed << std::setprecision(4)
                << (this->logspace_ ? pow(2.0, this->data_[i][a]) : this->data_[i][a]);
        // print neff
        out << '\t' << std::setprecision(2) << neff_[i] << std::endl;
    }

    out.flags(flags);
}

}  // cs

#endif
