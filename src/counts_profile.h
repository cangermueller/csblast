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
#include <string>
#include <vector>

#include "alignment.h"
#include "exception.h"
#include "matrix.h"
#include "profile.h"
#include "sequence.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class CountsProfile : public Profile<Alphabet_T>
{
  public:
    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::num_cols;
    using Profile<Alphabet_T>::alphabet_size;
    using Profile<Alphabet_T>::read;
    using Profile<Alphabet_T>::logspace;
    using Profile<Alphabet_T>::transform_to_logspace;
    using Profile<Alphabet_T>::transform_to_linspace;

    // Constructs profile from serialized profile read from input stream.
    explicit CountsProfile(std::istream& in);
    // Constructs a profile of the given sequence.
    explicit CountsProfile(const Sequence<Alphabet_T>& sequence);
    // Constructs a profile of the given alignment with specified sequence weighting method.
    explicit CountsProfile(const Alignment<Alphabet_T>& alignment, bool position_specific_weights = true);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    CountsProfile(const CountsProfile& other, int index, int length);

    virtual ~CountsProfile() { }

    // Reads all available profiles from the input stream and returns them in a vector.
    static void readall(std::istream& in, std::vector< shared_ptr<CountsProfile> >& v);
    // Returns the number of effective sequences in alignment column i
    float neff(int i) const { return neff_[i]; }
    // Converts the profile to counts of alphabet letters.
    void convert_to_counts();
    // Converts the profile back to relative frequencies of alphabet letters.
    void convert_to_frequencies();
    // Returns true if the profile contains counts.
    bool has_counts() const { return has_counts_; }

  protected:
    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::data_;
    using Profile<Alphabet_T>::SCALE_FACTOR;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized array data members to stream.
    virtual void write_body(std::ostream& out) const;
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

  private:
    // Return serialization class identity.
    virtual const std::string class_identity() const { static std::string id("CountsProfile"); return id; }

    // Number of effective sequences in each alignment column.
    std::vector<float> neff_;
    // Flag indicating if the profile contains counts or (relative) frequencies.
    bool has_counts_;
};  // CountsProfile



template<class Alphabet_T>
inline CountsProfile<Alphabet_T>::CountsProfile(std::istream& in)
        : neff_(),
          has_counts_(false)
{
    read(in);
}

template<class Alphabet_T>
inline CountsProfile<Alphabet_T>::CountsProfile(const Sequence<Alphabet_T>& sequence)
        : Profile<Alphabet_T>(sequence.length()),
          neff_(sequence.length(), 1.0f),
          has_counts_(false)
{
    for(int i = 0; i < num_cols(); ++i)
        data_[i][sequence[i]] = 1.0f;
}

template<class Alphabet_T>
inline CountsProfile<Alphabet_T>::CountsProfile(const Alignment<Alphabet_T>& alignment, bool position_specific_weights)
        : Profile<Alphabet_T>(alignment.num_match_cols()),
          neff_(alignment.num_match_cols()),
          has_counts_(false)
{
    const int num_cols = alignment.num_match_cols();
    const int num_seqs = alignment.num_seqs();
    const int any   = Alphabet_T::instance().any();

    if (position_specific_weights) {
        matrix<float> w;  // position-specific sequence weights
        neff_ = position_specific_weights_and_diversity(alignment, w);
        for (int i = 0; i < num_cols; ++i)
            for (int k = 0; k < num_seqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += w[i][k];
    } else {
        std::vector<float> wg;  // global sequence weights
        neff_.assign(num_cols, global_weights_and_diversity(alignment, wg));
        for (int i = 0; i < num_cols; ++i)
            for (int k = 0; k < num_seqs; ++k)
                if (alignment[i][k] < any)
                    data_[i][alignment[i][k]] += wg[k];
    }

    normalize(*this);
}

template<class Alphabet_T>
inline CountsProfile<Alphabet_T>::CountsProfile(const CountsProfile& other,
                                                int index,
                                                int length)
        : Profile<Alphabet_T>(other, index, length)
{
    neff_.insert(neff_.begin(), other.neff_.begin() + index, other.neff_.begin() + index + length);
    has_counts_ = other.has_counts_;
}

template<class Alphabet_T>
inline void CountsProfile<Alphabet_T>::readall(std::istream& in, std::vector< shared_ptr<CountsProfile> >& v)
{
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<CountsProfile> p(new CountsProfile(in));
        v.push_back(p);
    }
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::convert_to_counts()
{
    if (!has_counts_) {
        const bool islog = logspace();
        if (islog) transform_to_linspace();

        for (int i = 0; i < num_cols(); ++i)
            for (int a = 0; a < alphabet_size(); ++a)
                data_[i][a] *= neff_[i];
        has_counts_ = true;

        if (islog) transform_to_logspace();
    }
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::convert_to_frequencies()
{
    if (has_counts_) {
        normalize(*this);
        has_counts_ = false;
    }
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::read_header(std::istream& in)
{
    Profile<Alphabet_T>::read_header(in);
    neff_.resize(num_cols());

    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("has_counts") != std::string::npos)
        has_counts_ = atoi(tmp.c_str() + 10) == 1;
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::read_body(std::istream& in)
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
        for (int a = 0; a < alphabet_size(); ++a) {
            float log_p = tokens[a+1][0] == '*' ? std::numeric_limits<int>::max() : atoi(tokens[a+1].c_str());
            data_[i][a] = (logspace() ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR)) ;
        }
        neff_[i] = atof(tokens[alphabet_size()+1].c_str()) / SCALE_FACTOR;
        tokens.clear();
    }
    if (i != num_cols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, num_cols());
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::write_header(std::ostream& out) const
{
    Profile<Alphabet_T>::write_header(out);
    out << "has_counts\t" << has_counts_ << std::endl;
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::write_body(std::ostream& out) const
{
    out << "\t" << Alphabet_T::instance() << "\tNeff" << std::endl;
    for (int i = 0; i < num_cols(); ++i) {
        out << i+1;
        for (int a = 0; a < alphabet_size(); ++a) {
            float log_p = logspace() ? data_[i][a] : log2(data_[i][a]);
            if (-log_p == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(log_p * SCALE_FACTOR);
        }
        out << "\t" << iround(neff_[i] * SCALE_FACTOR) << std::endl;
    }
    out << "//" << std::endl;
}

template<class Alphabet_T>
void CountsProfile<Alphabet_T>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save old flags

    out << "\t" << Alphabet_T::instance() << "\tNeff" << std::endl;
    for (int i = 0; i < num_cols(); ++i) {
        out << i+1;
        for (int a = 0; a < alphabet_size(); ++a)
            out << '\t' << std::fixed << std::setprecision(4)
                << (logspace() ? pow(2.0, data_[i][a]) : data_[i][a]);
        // print neff
        out << '\t' << std::setprecision(2) << neff_[i] << std::endl;
    }

    out.flags(flags);
}

}  // cs

#endif
