#ifndef CS_CONTEXT_PROFILE_H
#define CS_CONTEXT_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for profiles derived from alignments.

#include <iostream>
#include <vector>

#include "alignment.h"
#include "exception.h"
#include "profile.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

template<class Alphabet_T>
class ContextProfile : public Profile<Alphabet_T>
{
  public:
    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::num_cols;
    using Profile<Alphabet_T>::alphabet_size;
    using Profile<Alphabet_T>::logspace;
    using Profile<Alphabet_T>::read;

     // Constructs a dummy context profile.
    ContextProfile();
    // Constructs a context profile with num_cols columns initialized to zero.
    ContextProfile(int index, int num_cols);
    // Constructs profile from serialized profile read from input stream.
    explicit ContextProfile(std::istream& in);
    // Constructs a context profile from simple profile and checks if length is valid.
    ContextProfile(int index, const Profile<Alphabet_T>& profile);

    virtual ~ContextProfile() {}

    // Reads all available profiles from the input stream and returns them in a vector.
    static void readall(std::istream& in, std::vector< shared_ptr<ContextProfile> >& v);

    // Returns the index of this context profile.
    int index() const { return index_; }
    // Sets the index of this context profile.
    void set_index(int i) { index_ = i; }
    // Returns the prior probability of this context profile.
    float prior() const { return prior_; }
    // Sets the prior probability of this context profile.
    void set_prior(float p) { prior_ = p; }
    // Returns index of central profile column.
    int center() const { return (num_cols() - 1) / 2; }

  protected:
    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::data_;
    using Profile<Alphabet_T>::logspace_;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

    // Index of context-profile.
    int index_;
    // Prior probability of context-profile.
    float prior_;

  private:
    // Checks if profile has odd number of columns.
    void check();

    // Return serialization class identity.
    virtual const std::string class_identity() { static std::string id("ContextProfile"); return id; }
};  // ContextProfile



template<class Alphabet_T>
inline ContextProfile<Alphabet_T>::ContextProfile()
        : Profile<Alphabet_T>(),
          index_(0),
          prior_(0.0f)
{ }

template<class Alphabet_T>
inline ContextProfile<Alphabet_T>::ContextProfile(int index, int num_cols)
        : Profile<Alphabet_T>(num_cols),
          index_(index),
          prior_(0.0f)
{
    check();
}

template<class Alphabet_T>
inline ContextProfile<Alphabet_T>::ContextProfile(int index, const Profile<Alphabet_T>& profile)
        : Profile<Alphabet_T>(profile),
          index_(index),
          prior_(0.0f)
{
    check();
}

template<class Alphabet_T>
inline ContextProfile<Alphabet_T>::ContextProfile(std::istream& in)
        : Profile<Alphabet_T>(),
          index_(0),
          prior_(0.0f)
{
    read(in);
    check();
}

template<class Alphabet_T>
inline void ContextProfile<Alphabet_T>::readall(std::istream& in, std::vector< shared_ptr<ContextProfile> >& v)
{
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<ContextProfile> p(new ContextProfile(in));
        v.push_back(p);
    }
}

template<class Alphabet_T>
void ContextProfile<Alphabet_T>::check()
{
    if (num_cols() % 2 != 1)
        throw Exception("Context profiles must have odd number of columns, but num_cols=%i!", num_cols());
}

template<class Alphabet_T>
void ContextProfile<Alphabet_T>::read_header(std::istream& in)
{
    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("index") != std::string::npos)
        index_ = atoi(tmp.c_str() + 5);
    // Read prior
    if (getline(in, tmp) && tmp.find("prior") != std::string::npos)
        prior_ = atof(tmp.c_str() + 5);

    Profile<Alphabet_T>::read_header(in);
}

template<class Alphabet_T>
void ContextProfile<Alphabet_T>::write_header(std::ostream& out) const
{
    out << "index\t\t" << index() << std::endl;
    out << "prior\t\t" << strprintf("%-10.8g", prior()) << std::endl;
    Profile<Alphabet_T>::write_header(out);
}

template<class Alphabet_T>
void ContextProfile<Alphabet_T>::print(std::ostream& out) const
{
    out << "Context-Profile " << index() << ":" << std::endl;
    out << "prior probability = " << strprintf("%-10.8g", prior()) << std::endl;
    Profile<Alphabet_T>::print(out);
}



template<class Alphabet_T>
void reset(ContextProfile<Alphabet_T>& profile)
{
    const int num_cols = profile.num_cols();
    const int alphabet_size = profile.alphabet_size();

    profile.set_prior(0.0f);
    for(int i = 0; i < num_cols; ++i)
        for(int a = 0; a < alphabet_size; ++a)
            profile[i][a] = 0.0f;
}

}  // cs

#endif
