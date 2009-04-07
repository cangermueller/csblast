#ifndef CS_PROFILE_H
#define CS_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A profile class representing columns of frequencies over a sequence alphabet.

#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
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
#include "utils.h"

namespace cs
{

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
    // Constructs profile from serialized profile read from input stream.
    explicit Profile(FILE* fin);
    // Creates a profile from subprofile starting at column index and length columns long.
    Profile(const Profile& other, int index, int length);

    virtual ~Profile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static void readall(std::istream& in, std::vector< shared_ptr<Profile> >& v);
    // Reads all available profiles from the input stream and returns them in a vector.
    static void readall(FILE* in, std::vector< shared_ptr<Profile> >* v);

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
    // Initializes the profile object with a serialized profile read from stream.
    void read(FILE* fin);
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
    static const int LOG_SCALE = 1000;
    // Buffer size for reading
    static const int BUFFER_SIZE = 1024;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(FILE* fin);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(FILE* fin);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized scalar data members to stream.
    virtual void write_header(FILE* fout) const;
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
    // Class identifier
    static const char* CLASS_ID;

    // Returns serialization class identity.
    virtual const std::string class_identity() const
    { static std::string id("Profile"); return id; }
    virtual const char* class_id() const { return CLASS_ID; }

};  // Profile



// Resets all entries in given profile to the given value or zero if none is given.
template<class Alphabet_T>
void reset(Profile<Alphabet_T>* p, float value = 0.0f);

// Normalize profile columns to value or to one if none provided.
template<class Alphabet_T>
bool normalize(Profile<Alphabet_T>* p, float value = 1.0f);



template<class Alphabet_T>
const char* Profile<Alphabet_T>::CLASS_ID = "Profile";

template<class Alphabet_T>
inline Profile<Alphabet_T>::Profile()
        : data_(),
          logspace_(false)
{ }

template<class Alphabet_T>
inline Profile<Alphabet_T>::Profile(int num_cols)
        : data_(num_cols, Alphabet_T::instance().size(), 0.0f),
          logspace_(false)
{ }

template<class Alphabet_T>
inline Profile<Alphabet_T>::Profile(std::istream& in)
        : data_(),
          logspace_(false)
{
    read(in);
}

template<class Alphabet_T>
inline Profile<Alphabet_T>::Profile(FILE* fin)
        : data_(),
          logspace_(false)
{
    read(fin);
}

template<class Alphabet_T>
inline Profile<Alphabet_T>::Profile(const Profile& other,
                                    int index,
                                    int length)
        : data_(length, other.alphabet_size(), 0.0f),
          logspace_(other.logspace_)
{
    if (index + length > other.num_cols())
        throw Exception("Arguments index=%i and length=%i of sub-profile are out of bounds!",
                        index, length);
    for (int i = 0; i < num_cols(); ++i)
        for (int a = 0; a < alphabet_size(); ++a)
            data_[i][a] = other[i+index][a];
}

template<class Alphabet_T>
inline void Profile<Alphabet_T>::readall(std::istream& in,
                                         std::vector< shared_ptr<Profile> >& v)
{
    while (in.peek() && in.good()) { //peek first to make sure that we don't read beyond '//'
        shared_ptr<Profile> p(new Profile(in));
        v.push_back(p);
    }
}

template<class Alphabet_T>
inline void Profile<Alphabet_T>::readall(FILE* fin, std::vector< shared_ptr<Profile> >* v)
{
    while (!feof(fin)) {
        shared_ptr<Profile> p(new Profile(fin));
        v->push_back(p);

        int c = getc(fin);
        if (c == EOF) break;
        ungetc(c, fin);
    }
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
        throw Exception("Bad format: profile does not start with '%s'!", class_id());

    read_header(in);
    read_body(in);

    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read(FILE* fin)
{
    LOG(DEBUG1) << "Reading profile from stream ...";

    // Check if stream actually contains a serialized profile
    char buffer[BUFFER_SIZE];
    while (fgetline(buffer, BUFFER_SIZE, fin))
        if (strscn(buffer)) break;
    if (!strstr(buffer, class_id()))
        throw Exception("Bad format: profile does not start with '%s'!", class_id());

    read_header(fin);
    read_body(fin);

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
        throw Exception("Bad format: alphabet_size=%i does not fit with alphabet size %i!",
                        alphabet_size, Alphabet_T::instance().size());

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    resize(num_cols, alphabet_size);
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read_header(FILE* fin)
{
    char buffer[BUFFER_SIZE];
    const char* ptr = buffer;

    // Read num_cols
    int num_cols = 0;
    if (fgetline(buffer, BUFFER_SIZE, fin) && strstr(buffer, "num_cols")) {
        ptr = buffer;
        num_cols = strtoi(ptr);
    } else {
        throw Exception("Bad format: profile does not contain 'num_cols' record!");
    }
    // Read alphabet_size
    int alphabet_size = 0;
    if (fgetline(buffer, BUFFER_SIZE, fin) && strstr(buffer, "alphabet_size")) {
        ptr = buffer;
        alphabet_size = strtoi(ptr);
    } else {
        throw Exception("Bad format: profile does not contain 'alphabet_size' record!");
    }
    if (alphabet_size != Alphabet_T::instance().size())
        throw Exception("Bad format: profile alphabet_size should be %i but is %i!",
                        Alphabet_T::instance().size(), alphabet_size);
    // Read logspace
    if (fgetline(buffer, BUFFER_SIZE, fin) && strstr(buffer, "logspace")) {
        ptr = buffer;
        logspace_ = strtoi(ptr) == 1;
    } else {
        throw Exception("Bad format: profile does not contain 'logspace' record!");
    }

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
            float log_p =
                tokens[a+1][0] == '*' ? INT_MAX : atoi(tokens[a+1].c_str());
            data_[i][a] = logspace_ ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR);
        }
        tokens.clear();
    }
    if (i != num_cols() - 1)
        throw Exception("Bad format: profile has %i column records but should have %i!",
                        i+1, num_cols());
}

template<class Alphabet_T>
void Profile<Alphabet_T>::read_body(FILE* fin)
{
    const int alph_size = alphabet_size();
    char buffer[BUFFER_SIZE];
    const char* ptr = buffer;
    int i = 0;

    fgetline(buffer, BUFFER_SIZE, fin);  // skip alphabet description line
    while (fgetline(buffer, BUFFER_SIZE, fin) && buffer[0] != '/' && buffer[1] != '/') {
        ptr = buffer;
        i = strtoi(ptr) - 1;
        for (int a = 0; a < alph_size; ++a) {
            if (logspace_)
                data_[i][a] = static_cast<float>(-strtoi_ast(ptr)) / LOG_SCALE;
            else
                data_[i][a] = pow(2.0, static_cast<float>(-strtoi_ast(ptr)) / LOG_SCALE);
        }
    }
    if (i != num_cols() - 1)
        throw Exception("Bad format: profile has %i columns but should have %i!",
                        i+1, num_cols());
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
    out << "num_cols\t" << num_cols() << std::endl;
    out << "alphabet_size\t" << alphabet_size() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;
}

template<class Alphabet_T>
void Profile<Alphabet_T>::write_header(FILE* fout) const
{
    fprintf(fout, "num_cols\t%i\n", num_cols());
    fprintf(fout, "alphabet_size\t%i\n", alphabet_size());
    fprintf(fout, "logspace\t%i\n", logspace() ? 1 : 0);
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
inline void reset(Profile<Alphabet_T>* p)
{
    Profile<Alphabet_T>& profile = *p;
    const int num_cols = profile.num_cols();
    const int alphabet_size = profile.alphabet_size();
    for(int i = 0; i < num_cols; ++i)
        for(int a = 0; a < alphabet_size; ++a)
            profile[i][a] = 0.0f;
}

template<class Alphabet_T>
bool normalize(Profile<Alphabet_T>* p, float value)
{
    Profile<Alphabet_T>& profile = *p;
    const bool logspace = profile.logspace();
    if (logspace) profile.transform_to_linspace();

    const int num_cols       = profile.num_cols();
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
