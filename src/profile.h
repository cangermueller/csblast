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

template<class AlphabetType>
class Profile
{
  public:
    typedef matrix<float>::row_type col_type;
    typedef matrix<float>::const_row_type const_col_type;
    typedef matrix<float>::iterator iterator;
    typedef matrix<float>::const_iterator const_iterator;

    // Constructs a dummy profile.
    Profile();
    // Constructs a profile with ncols columns initialized to zero.
    Profile(int ncols);
    // Constructs profile from serialized profile read from input stream.
    Profile(std::istream& in);
    // Creates a profile from subprofile in other, starting at column index and length columns long.
    Profile(const Profile& other, int index, int length);

    virtual ~Profile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static std::vector< shared_ptr<Profile> > readall(std::istream& in);

    // Access methods to get the (i,j) element
    col_type operator[](int i) { return data_[i]; }
    const_col_type operator[](int i) const { return data_[i]; }
    // Returns #columns in the profile
    int ncols() const { return data_.nrows(); }
    // Returns #entries per column
    int nalph() const { return data_.ncols(); }
    // Returns the total number of elements in the profile.
    int size() const { return data_.size(); }
    // Transforms profile to logspace
    virtual void transform_to_logspace();
    // Transforms profile to linspace
    virtual void transform_to_linspace();
    // Returns true if the profile is in logspace
    bool logspace() const { return logspace_; }
    // Sets logspace flag.
    void set_logspace(bool flag) { logspace_ = flag; }
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
    static const int kScaleFactor = 1000;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized array data members to stream.
    virtual void write_body(std::ostream& out) const;
    // Resize the profile matrix to given dimensions. Attention: old data is lost!
    void resize(int ncols, int nalph);
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

     // Profile matrix in row major format
    matrix<float> data_;
    // Flag indicating if profile is in log- or linspace
    bool logspace_;

  private:
    // Disallow copy and assign
    Profile(const Profile&);
    void operator=(const Profile&);

    // Returns serialization class identity.
    virtual const std::string class_identity() const { static std::string id("Profile"); return id;}
};  // Profile



// Resets all entries in given profile to the given value or zero if none is given.
template<class AlphabetType>
void reset(Profile<AlphabetType>& profile, float value);

// Normalize profile columns to value or to one if none provided.
template<class AlphabetType>
void normalize(Profile<AlphabetType>& profile, float value);



template<class AlphabetType>
Profile<AlphabetType>::Profile()
        : logspace_(false)
{}

template<class AlphabetType>
Profile<AlphabetType>::Profile(int ncols)
        : data_(ncols, AlphabetType::instance().size(), 0.0f),
          logspace_(false)
{}

template<class AlphabetType>
Profile<AlphabetType>::Profile(std::istream& in)
{ read(in); }

template<class AlphabetType>
Profile<AlphabetType>::Profile(const Profile& other,
                 int index,
                 int length)
        : data_(length, other.nalph(), 0.0f)
{
    if (index + length > other.ncols())
        throw Exception("Arguments index=%i and length=%i of sub-profile are out of bounds!", index, length);
    for (int i = 0; i < ncols(); ++i)
        for (int a = 0; a < nalph(); ++a)
            data_[i][a] = other[i+index][a];
}

template<class AlphabetType>
std::vector< shared_ptr< Profile<AlphabetType> > > Profile<AlphabetType>::readall(std::istream& in)
{
    std::vector< shared_ptr<Profile> > profiles;
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<Profile> p(new Profile(in));
        profiles.push_back(p);
    }

    return profiles;
}

template<class AlphabetType>
void Profile<AlphabetType>::transform_to_logspace()
{
    for(int i = 0; i < ncols(); ++i)
        for(int a = 0; a < nalph(); ++a)
            data_[i][a] = data_[i][a] == 0.0f ? -std::numeric_limits<float>::infinity() : log2(data_[i][a]);
    logspace_ = true;
}

template<class AlphabetType>
void Profile<AlphabetType>::transform_to_linspace()
{
    for(int i = 0; i < ncols(); ++i)
        for(int a = 0; a < nalph(); ++a)
            data_[i][a] = pow(2.0, data_[i][a]);
    logspace_ = false;
}

template<class AlphabetType>
void Profile<AlphabetType>::read(std::istream& in)
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

template<class AlphabetType>
void Profile<AlphabetType>::read_header(std::istream& in)
{
    std::string tmp;

    // Read ncols
    int ncols = 0;
    if (getline(in, tmp) && tmp.find("ncols") != std::string::npos)
        ncols = atoi(tmp.c_str() + 5);
    else
        throw Exception("Bad format: serialized profile does not contain 'ncols' record!");

    // Read nalph
    int nalph = 0;
    if (getline(in, tmp) && tmp.find("nalph") != std::string::npos)
        nalph = atoi(tmp.c_str() + 5);
    else
        throw Exception("Bad format: serialized profile does not contain 'nalph' record!");
    if (nalph != AlphabetType::instance().size())
        throw Exception("Bad format: nalph=%i does not fit with alphabet size %i!", nalph, AlphabetType::instance().size());

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    resize(ncols, nalph);
}

template<class AlphabetType>
void Profile<AlphabetType>::read_body(std::istream& in)
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
        tokens.clear();
    }
    if (i != ncols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!", i+1, ncols());
}

template<class AlphabetType>
void Profile<AlphabetType>::write(std::ostream& out) const
{
    out << class_identity() << std::endl;
    write_header(out);
    write_body(out);
}

template<class AlphabetType>
void Profile<AlphabetType>::write_header(std::ostream& out) const
{
    // print dimensions
    out << "ncols\t\t" << ncols() << std::endl;
    out << "nalph\t\t" << nalph() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;
}

template<class AlphabetType>
void Profile<AlphabetType>::write_body(std::ostream& out) const
{
    out << "\t" << AlphabetType::instance() << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < nalph(); ++a) {
            float logval = logspace_ ? data_[i][a] : log2(data_[i][a]);
            if (-logval == std::numeric_limits<float>::infinity())
                out << "\t*";
            else
                out << "\t" << -iround(logval * kScaleFactor);
        }
        out << std::endl;
    }
}

template<class AlphabetType>
void Profile<AlphabetType>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save old flags

    out << "\t" << AlphabetType::instance().size() << std::endl;
    for (int i = 0; i < ncols(); ++i) {
        out << i+1;
        for (int a = 0; a < nalph(); ++a)
            out << '\t' << std::fixed << std::setprecision(4)
                << (logspace_ ? pow(2.0, data_[i][a]) : data_[i][a]);
        // print neff
        out << std::endl;
    }

    out.flags(flags);
}

template<class AlphabetType>
void Profile<AlphabetType>::resize(int ncols, int nalph)
{
    if (ncols == 0 || nalph == 0)
        throw Exception("Bad dimensions for profile resizing: ncols=%i nalph=%i", ncols, nalph);
    data_.resize(ncols, nalph);
}

template<class AlphabetType>
void reset(Profile<AlphabetType>& profile, float value)
{
    const int ncols = profile.ncols();
    const int nalph = profile.nalph();
    for(int i = 0; i < ncols; ++i)
        for(int a = 0; a < nalph; ++a)
            profile[i][a] = value;
}

template<class AlphabetType>
void normalize(Profile<AlphabetType>& profile, float value)
{
    const bool logspace = profile.logspace();
    if (logspace) profile.transform_to_linspace();

    const int ncols = profile.ncols();
    const int nalph  = profile.nalph();
    for (int i = 0; i < ncols; ++i) {
        float sum = 0.0f;
        for (int a = 0; a < nalph; ++a) sum += profile[i][a];
        if (sum != 0.0f) {
            float fac = value / sum;
            for (int a = 0; a < nalph; ++a) profile[i][a] *= fac;
        } else {
            throw Exception("Unable to normalize profile to one. Sum of column %i is zero!", i);
        }
    }

    if (logspace) profile.transform_to_logspace();
}

}  // cs

#endif
