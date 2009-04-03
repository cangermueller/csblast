#ifndef CS_PROFILE_LIBRARY_H
#define CS_PROFILE_LIBRARY_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for context profiles to represent the most common
// sequence motifs in a training database of proteins/DNA sequences.

#include <cstdlib>
#include <ctime>

#include <iostream>

#include "exception.h"
#include "profile.h"
#include "context_profile.h"
#include "counts_profile.h"
#include "log.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class ProfileLibrary;


template<class Alphabet_T>
class ProfileInitializer
{
  public:
    ProfileInitializer() {}
    virtual ~ProfileInitializer() {};
    virtual void init(ProfileLibrary<Alphabet_T>& lib) const = 0;
};


template<class Alphabet_T>
class ProfileLibrary
{
  public:
    // Public typedefs
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > >::iterator profile_iterator;
    typedef typename std::vector< shared_ptr< ContextProfile<Alphabet_T> > >::const_iterator const_profile_iterator;

    // Constructs an empty profile library of given dimenions.
    ProfileLibrary(int num_profiles, int num_cols);
    // Constructs a profile library from serialized data read from input stream.
    ProfileLibrary(std::istream& in);
    // Constructs profile library with a specific init-strategy encapsulated by an initializer.
    ProfileLibrary(int num_profiles,
                   int num_cols,
                   const ProfileInitializer<Alphabet_T>& profile_init);

    virtual ~ProfileLibrary() {}

    bool full() const { return static_cast<int>(profiles_.size()) == num_profiles_; }
    // Returns the number of profiles in the profile library
    int num_profiles() const { return num_profiles_; }
    // Returns the number of profiles in the profile library
    int size() const { return num_profiles_; }
    // Returns the number of columns in each context profile.
    int num_cols() const { return num_cols_; }
    // Returns the size of the alphabet of the profile library.
    int alphabet_size() const { return Alphabet_T::instance().size();  }
    // Returns the number of clustering iterations.
    int iterations() const { return iterations_; }
    // Accessor methods for state i, where i is from interval [0,num_profiles].
    ContextProfile<Alphabet_T>& operator[](int i) { return *profiles_[i]; }
    const ContextProfile<Alphabet_T>& operator[](int i) const { return *profiles_[i]; }
    // Clears the library.
    void clear();
    // Adds the given profile to the library and returns its profile index. Note that number of profile columns must be odd!
    int add_profile(const Profile<Alphabet_T>& profile);
    // Returns an iterator to a list of pointers to context profiles.
    profile_iterator begin() { return profiles_.begin(); }
    // Returns an iterator pointing past the end of a list of pointers of states.
    profile_iterator end() { return profiles_.end(); }
    // Returns a const iterator to a list of pointers of states.
    const_profile_iterator begin() const { return profiles_.begin(); }
    // Returns a const iterator pointing past the end of a list of pointers of states.
    const_profile_iterator end() const { return profiles_.end(); }
    // Writes the profile library in serialization format to output stream.
    void write(std::ostream& out) const;
    // Returns true if state profiles are in logspace.
    bool logspace() const { return logspace_; }
    // Transforms profiles to logspace.
    void transform_to_logspace();
    // Transforms profiles to linspace.
    void transform_to_linspace();
    // Increments the EM-clustering iteration counter.
    ProfileLibrary& operator++() { ++iterations_; return *this; }

    // Prints the library in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const ProfileLibrary& lib)
    {
        lib.print(out);
        return out;
    }

  private:
    // Prints the library in human-readable format to output stream.
    void print(std::ostream& out) const;
    // Initializes the library from serialized data read from stream.
    void read(std::istream& in);

    // Number of profiles in the fully assembled library.
    int num_profiles_;
    // Number of columns in each context profile.
    int num_cols_;
    // Number of EM-clustering iterations performed on this profile library.
    int iterations_;
    // Context profiles ordered by index.
    std::vector< shared_ptr< ContextProfile<Alphabet_T> > > profiles_;
    // Flag indicating if profile probabilities are in log- or linspace
    bool logspace_;
};  // ProfileLibrary



template<class Alphabet_T>
ProfileLibrary<Alphabet_T>::ProfileLibrary(int num_profiles, int num_cols)
        : num_profiles_(num_profiles),
          num_cols_(num_cols),
          iterations_(0),
          profiles_(),
          logspace_(false)
{
    profiles_.reserve(num_profiles);
}

template<class Alphabet_T>
ProfileLibrary<Alphabet_T>::ProfileLibrary(std::istream& in)
        : num_profiles_(0),
          num_cols_(0),
          iterations_(0),
          profiles_(),
          logspace_(false)
{
    read(in);
}

template<class Alphabet_T>
ProfileLibrary<Alphabet_T>::ProfileLibrary(int num_profiles,
                                           int num_cols,
                                           const ProfileInitializer<Alphabet_T>& profile_init)
        : num_profiles_(num_profiles),
          num_cols_(num_cols),
          iterations_(0),
          profiles_(),
          logspace_(false)
{
    profiles_.reserve(num_profiles);
    profile_init.init(*this);
}

template<class Alphabet_T>
inline void ProfileLibrary<Alphabet_T>::clear()
{
    profiles_.clear();
    profiles_.reserve(num_profiles());
}

template<class Alphabet_T>
inline int ProfileLibrary<Alphabet_T>::add_profile(const Profile<Alphabet_T>& profile)
{
    if (full())
        throw Exception("Unable to add profile: profile library contains already %i profiles!", num_profiles());
    if (profile.num_cols() != num_cols())
        throw Exception("Unable to add profile: provided profile has %i columns but should have %i!",
                        profile.num_cols(), num_cols());

    shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(profiles_.size(), profile));
    profile_ptr->set_prior(1.0f / num_profiles());  // start with unbiased prior probabilities

    profiles_.push_back(profile_ptr);
    return profiles_.size() - 1;
}

template<class Alphabet_T>
inline void ProfileLibrary<Alphabet_T>::transform_to_logspace()
{
    if (!logspace()) {
        for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
            (*pi)->transform_to_logspace();
        logspace_ = true;
    }
}

template<class Alphabet_T>
inline void ProfileLibrary<Alphabet_T>::transform_to_linspace()
{
    if (logspace()) {
        for (profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
            (*pi)->transform_to_linspace();
        logspace_ = false;
    }
}

template<class Alphabet_T>
void ProfileLibrary<Alphabet_T>::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading profile library from stream ...";

    // Check if stream actually contains a serialized profile library
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find("ProfileLibrary") == std::string::npos)
        throw Exception("Bad format: serialized profile library does not start with 'ProfileLibrary' keyword!");
    // Read number of profiles
    if (getline(in, tmp) && tmp.find("num_profiles") != std::string::npos)
        num_profiles_ = atoi(tmp.c_str() + 12);
    else
        throw Exception("Bad format: serialized profile library does not contain 'num_profiles' record!");
    // Read number of columns
    if (getline(in, tmp) && tmp.find("num_cols") != std::string::npos)
        num_cols_= atoi(tmp.c_str() + 8);
    else
        throw Exception("Bad format: serialized profile library does not contain 'num_cols' record!");
    // Read number of iterations
    if (getline(in, tmp) && tmp.find("iterations") != std::string::npos)
        iterations_= atoi(tmp.c_str() + 10);
    else
        throw Exception("Bad format: serialized profile library does not contain 'num_cols' record!");
    // Read profiles_logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;
    // Read profile records
    profiles_.reserve(num_profiles());
    while (!full() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr< ContextProfile<Alphabet_T> > profile_ptr(new ContextProfile<Alphabet_T>(in));
        profiles_.push_back(profile_ptr);
    }
    if (!full())
        throw Exception("Error while reading profile library: number of profile records is %i but should be %i!",
                        profiles_.size(), num_profiles());

    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void ProfileLibrary<Alphabet_T>::write(std::ostream& out) const
{
    // Write header
    out << "ProfileLibrary"   << std::endl;
    out << "num_profiles\t\t" << num_profiles() << std::endl;
    out << "num_cols\t\t"     << num_cols() << std::endl;
    out << "iterations\t\t"   << iterations() << std::endl;
    out << "logspace\t\t"     << (logspace() ? 1 : 0) << std::endl;

    // Serialize profiles
    for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
        (*pi)->write(out);
}

template<class Alphabet_T>
void ProfileLibrary<Alphabet_T>::print(std::ostream& out) const
{
    out << "ProfileLibrary" << std::endl;
    out << "Total number of profiles: " << num_profiles() << std::endl;
    out << "Context profile columns:  " << num_cols() << std::endl;
    out << "Clustering iterations:    " << iterations() << std::endl;

    for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end(); ++pi)
        out << **pi;
}



template<class Alphabet_T>
class SamplingProfileInitializer : public ProfileInitializer<Alphabet_T>
{
  public:
    typedef typename std::vector< shared_ptr< CountsProfile<Alphabet_T> > > profile_vector;
    typedef typename profile_vector::const_iterator profile_iterator;

    SamplingProfileInitializer(profile_vector profiles,
                               const Pseudocounts<Alphabet_T>* pc = NULL,
                               float pc_admixture = 1.0f)
            : profiles_(profiles),
              pc_(pc),
              pc_admixture_(pc_admixture)
    {
        random_shuffle(profiles_.begin(), profiles_.end());
    }

    virtual ~SamplingProfileInitializer() { };

    virtual void init(ProfileLibrary<Alphabet_T>& lib) const
    {
        LOG(DEBUG) << "Initializing profile library with " << lib.num_profiles() << " profile windows randomly sampled from "
                  << profiles_.size() << " training profiles ...";

        // Iterate over randomly shuffled profiles and add them to the library until the library is full.
        for (profile_iterator pi = profiles_.begin(); pi != profiles_.end() && !lib.full(); ++pi) {
            if (lib.num_cols() != (*pi)->num_cols())
                throw Exception("Unable to initialize profiles: library profiles have length %i but found training profile with length %i!",
                                lib.num_cols(), (*pi)->num_cols());

            CountsProfile<Alphabet_T> p(**pi);
            p.convert_to_frequencies(); // make sure that profile contains frequencies not counts
            if (pc_) pc_->add_to_profile(ConstantAdmixture(pc_admixture_), &p);
            lib.add_profile(p);
        }
        if (!lib.full())
            throw Exception("Could not fully initialize all %i library profiles. Maybe too few training profiles provided?",
                            lib.num_profiles());

        LOG(DEBUG) << "Profile library after profile initialization:";
        LOG(DEBUG) << lib;
    }

  private:
    // Pool of profile windows to sample from.
    profile_vector profiles_;
    // Pseudocount factory for library profiles.
    const Pseudocounts<Alphabet_T>* pc_;
    // Constant pseudocount admixture for library profiles.
    float pc_admixture_;
};  // SamplingProfileInitializer

}  // cs

#endif
