#ifndef CS_PROFILE_H
#define CS_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies over a sequence alphabet.

#include <cmath>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "exception.h"
#include "log.h"
#include "matrix.h"
#include "profile.h"
#include "shared_ptr.h"
#include "util.h"

namespace cs
{

// Forward declarations
template<class T>
class Sequence;

template<class Alphabet_T>
class Profile
{
  public:
    typedef matrix<float>::row_type col_type;
    typedef matrix<float>::const_row_type const_col_type;
    typedef matrix<float>::iterator iterator;
    typedef matrix<float>::const_iterator const_iterator;

    // Constructs a dummy profile.
    Profile();
    // Constructs a profile with num_cols columns initialized to zero.
    explicit Profile(int num_cols);
    // Constructs profile from serialized profile read from input stream.
    explicit Profile(std::istream& in);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    Profile(const Profile& other, int index, int length);

    virtual ~Profile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<Profile> > readall(std::istream& in);

    // Access methods to get the (i,j) element
    col_type operator[](int i) { return data_[i]; }
    const_col_type operator[](int i) const { return data_[i]; }
    // Returns #columns in the profile
    int num_cols() const { return data_.num_rows(); }
    // Returns #columns in the profile
    int length() const { return data_.num_rows(); }
    // Returns #entries per column
    int alphabet_size() const { return data_.num_cols(); }
    // Returns the total number of elements in the profile.
    int size() const { return data_.size(); }
    // Transforms profile to logspace
    virtual void transform_to_logspace();
    // Transforms profile to linspace
    virtual void transform_to_linspace();
    // Returns true if the profile is in logspace
    bool logspace() const { return logspace_; }
    // Returns an iterator to the first element in profile column i.
    col_type col_begin(int i) { return data_.row_begin(i); }
    // Returns an iterator just past the end of profile column i.
    col_type col_end(int i) { return data_.row_end(i); }
    // Returns a const iterator to the first element in profile column i.
    const_col_type col_begin(int i) const { return data_.row_begin(i); }
    // Returns a const iterator just past the end of profile column i.
    const_col_type col_end(int i) const { return data_.row_end(i); }
    // Returns an iterator to the first element in the profile matrix.
    iterator begin() { return data_.begin(); }
    // Returns an iterator just past the end of the profile matrix.
    iterator end() { return data_.end(); }
    // Returns a const iterator to the first element in the profile matrix.
    const_iterator begin() const { return data_.begin(); }
    // Returns a const iterator just past the end of the profile matrix.
    const_iterator end() const { return data_.end(); }
    // Initializes the profile object with a serialized profile read from stream.
    void read(std::istream& in);
    // Writes the profile in serialization format to output stream.
    void write(std::ostream& out) const;

    // Prints profile in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const Profile& p)
    {
        p.print(out);
        return out;
    }

  protected:
    // Scaling factor for serialization of profile log values
    static const int SCALE_FACTOR = 1000;

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
    // Resize the profile matrix to given dimensions. Attention: old data is lost!
    void resize(int num_cols, int alphabet_size);

     // Profile matrix in row major format
    matrix<float> data_;
    // Flag indicating if profile is in log- or linspace
    bool logspace_;

  private:
    // Returns serialization class identity.
    virtual const std::string class_identity() const { static std::string id("Profile"); return id;}
};  // Profile



// Resets all entries in given profile to the given value or zero if none is given.
template<class Alphabet_T>
void reset(Profile<Alphabet_T>& profile, float value = 0.0f);

// Normalize profile columns to value or to one if none provided.
template<class Alphabet_T>
bool normalize(Profile<Alphabet_T>& profile, float value = 1.0f);



template<class Alphabet_T>
Profile<Alphabet_T>::Profile()
        : data_(),
          logspace_(false)
{}

template<class Alphabet_T>
Profile<Alphabet_T>::Profile(int num_cols)
        : data_(num_cols, Alphabet_T::instance().size(), 0.0f),
          logspace_(false)
{}

template<class Alphabet_T>
Profile<Alphabet_T>::Profile(std::istream& in)
        : data_(),
          logspace_(false)
{
    read(in);
}

template<class Alphabet_T>
Profile<Alphabet_T>::Profile(const Profile& other,
                 int index,
                 int length)
        : data_(length, other.alphabet_size(), 0.0f),
          logspace_(other.logspace_)
{
    if (index + length > other.num_cols())
        throw Exception("Arguments index=%i and length=%i of sub-profile are out of bounds!", index, length);
    for (int i = 0; i < num_cols(); ++i)
        for (int a = 0; a < alphabet_size(); ++a)
            data_[i][a] = other[i+index][a];
}

template<class Alphabet_T>
std::vector< shared_ptr< Profile<Alphabet_T> > > Profile<Alphabet_T>::readall(std::istream& in)
{
    std::vector< shared_ptr<Profile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<Profile> p(new Profile(in));
        profiles.push_back(p);
    }

    return profiles;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::transform_to_logspace()
{
    if (!logspace_) {
        for(int i = 0; i < num_cols(); ++i)
            for(int a = 0; a < alphabet_size(); ++a)
                data_[i][a] = log2(data_[i][a]);
        logspace_ = true;
    }
}

template<class Alphabet_T>
void Profile<Alphabet_T>::transform_to_linspace()
{
    if (logspace_) {
        for(int i = 0; i < num_cols(); ++i)
            for(int a = 0; a < alphabet_size(); ++a)
                data_[i][a] = pow(2.0, data_[i][a]);
        logspace_ = false;
    }
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading profile from stream ...";

    // Check if stream actually contains a serialized profile
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find(class_identity()) == std::string::npos)
        throw Exception("Bad format: serialized profile does not start with '%s'!", class_identity().c_str());

    read_header(in);
    read_body(in);

    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read_header(std::istream& in)
{
    std::string tmp;

    // Read num_cols
    int num_cols = 0;
    if (getline(in, tmp) && tmp.find("num_cols") != std::string::npos)
        num_cols = atoi(tmp.c_str() + 8);
    else
        throw Exception("Bad format: serialized profile does not contain 'num_cols' record!");

    // Read alphabet_size
    int alphabet_size = 0;
    if (getline(in, tmp) && tmp.find("alphabet_size") != std::string::npos)
        alphabet_size = atoi(tmp.c_str() + 13);
    else
        throw Exception("Bad format: serialized profile does not contain 'alphabet_size' record!");
    if (alphabet_size != Alphabet_T::instance().size())
        throw Exception("Bad format: alphabet_size=%i does not fit with alphabet size %i!", alphabet_size, Alphabet_T::instance().size());

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    resize(num_cols, alphabet_size);
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read_body(std::istream& in)
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
            data_[i][a] = logspace_ ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR);
        }
        tokens.clear();
    }
    if (i != num_cols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, num_cols());
}

template<class Alphabet_T>
void Profile<Alphabet_T>::write(std::ostream& out) const
{
    out << class_identity() << std::endl;
    write_header(out);
    write_body(out);
}

template<class Alphabet_T>
void Profile<Alphabet_T>::write_header(std::ostream& out) const
{
    // print dimensions
    out << "num_cols\t" << num_cols() << std::endl;
    out << "alphabet_size\t" << alphabet_size() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::write_body(std::ostream& out) const
{
    out << "\t" << Alphabet_T::instance() << std::endl;
    for (int i = 0; i < num_cols(); ++i) {
        out << i+1;
        for (int a = 0; a < alphabet_size(); ++a) {
            float logval = logspace_ ? data_[i][a] : log2(data_[i][a]);
            if (-logval == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * SCALE_FACTOR);
        }
        out << std::endl;
    }
    out << "//" << std::endl;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save flags

    out << "\t" << Alphabet_T::instance() << std::endl;
    for (int i = 0; i < num_cols(); ++i) {
        out << i+1;
        for (int a = 0; a < alphabet_size(); ++a)
            out << '\t' << std::fixed << std::setprecision(4)
                << (logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
        // print neff
        out << std::endl;
    }

    out.flags(flags);
}

template<class Alphabet_T>
void Profile<Alphabet_T>::resize(int num_cols, int alphabet_size)
{
    if (num_cols == 0 || alphabet_size == 0)
        throw Exception("Bad dimensions for profile resizing: num_cols=%i alphabet_size=%i", num_cols, alphabet_size);
    data_.resize(num_cols, alphabet_size);
}

template<class Alphabet_T>
void reset(Profile<Alphabet_T>& profile, float value)
{
    const int num_cols = profile.num_cols();
    const int alphabet_size = profile.alphabet_size();
    for(int i = 0; i < num_cols; ++i)
        for(int a = 0; a < alphabet_size; ++a)
            profile[i][a] = value;
}

template<class Alphabet_T>
bool normalize(Profile<Alphabet_T>& profile, float value)
{
    const bool logspace = profile.logspace();
    if (logspace) profile.transform_to_linspace();

    const int num_cols = profile.num_cols();
    const int alphabet_size  = profile.alphabet_size();
    bool rv = true;
    for (int i = 0; i < num_cols; ++i) {
        float sum = 0.0f;
        for (int a = 0; a < alphabet_size; ++a) sum += profile[i][a];
        if (sum != 0.0f) {
            float fac = value / sum;
            for (int a = 0; a < alphabet_size; ++a) profile[i][a] *= fac;
        } else {
            rv = false;  // couldn't normalize at least one column
        }
    }

    if (logspace) profile.transform_to_logspace();
    return rv;
}

}  // cs

#endif
