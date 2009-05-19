// Copyright 2009, Andreas Biegert

#ifndef SRC_PROFILE_LIBRARY_H_
#define SRC_PROFILE_LIBRARY_H_

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <vector>

#include "globals.h"
#include "profile.h"
#include "context_profile.h"
#include "count_profile.h"
#include "pseudocounts.h"
#include "shared_ptr.h"

namespace cs {

// Forward declarations
template<class Alphabet>
class ProfileLibrary;

template<class Alphabet>
class ProfileInitializer {
 public:
  ProfileInitializer() {}
  virtual ~ProfileInitializer() {}
  virtual void Init(ProfileLibrary<Alphabet>& lib) const = 0;
};


// A container class for context profiles to represent the most common
// sequence motifs in a training database of proteins/DNA sequences.
template<class Alphabet>
class ProfileLibrary {
 public:
  // Public typedefs
  typedef std::vector< shared_ptr< ContextProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::iterator profile_iterator;
  typedef typename profile_vector::const_iterator const_profile_iterator;

  // Constructs an empty profile library of given dimenions.
  ProfileLibrary(int num_profiles, int num_cols);
  // Constructs a profile library from serialized data read from input stream.
  explicit ProfileLibrary(FILE* fin);
  // Constructs profile library with a specific init-strategy encapsulated by an
  // initializer.
  ProfileLibrary(int num_profiles,
                 int num_cols,
                 const ProfileInitializer<Alphabet>& profile_init);

  virtual ~ProfileLibrary() {}

  bool full() const {
    return static_cast<int>(profiles_.size()) == num_profiles_;
  }
  // Returns the number of profiles in the profile library
  int num_profiles() const { return num_profiles_; }
  // Returns the number of profiles in the profile library
  int size() const { return num_profiles_; }
  // Returns the number of columns in each context profile.
  int num_cols() const { return num_cols_; }
  // Returns index of central profile column.
  int center() const { return (num_cols() - 1) / 2; }
  // Returns the size of the alphabet of the profile library.
  int alphabet_size() const { return Alphabet::instance().size(); }
  // Returns the number of clustering iterations.
  int iterations() const { return iterations_; }
  // Accessor methods for state i, where i is from interval [0,num_profiles].
  ContextProfile<Alphabet>& operator[](int i) { return *profiles_[i]; }
  const ContextProfile<Alphabet>& operator[](int i) const {
    return *profiles_[i];
  }
  // Clears the library.
  void clear();
  // Adds the given profile to the library and returns its profile index.
  int AddProfile(const Profile<Alphabet>& profile);
  // Returns an iterator to a list of pointers to context profiles.
  profile_iterator begin() { return profiles_.begin(); }
  // Returns an iterator pointing past the end of a list of pointers of states.
  profile_iterator end() { return profiles_.end(); }
  // Returns a const iterator to a list of pointers of states.
  const_profile_iterator begin() const { return profiles_.begin(); }
  // Returns a const iterator pointing past the end of a list of pointers of
  // states.
  const_profile_iterator end() const { return profiles_.end(); }
  // Writes the profile library in serialization format to output stream.
  void Write(FILE* fout) const;
  // Returns true if state profiles are in logspace.
  bool logspace() const { return logspace_; }
  // Transforms profiles to logspace.
  void transform_to_logspace();
  // Transforms profiles to linspace.
  void transform_to_linspace();
  // Increments the EM-clustering iteration counter.
  void increment_iterations() { ++iterations_; }

  // Prints the library in human-readable format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const ProfileLibrary& lib) {
    lib.print(out);
    return out;
  }

 private:
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Prints the library in human-readable format to output stream.
  void print(std::ostream& out) const;
  // Initializes the library from serialized data read from stream.
  void Read(FILE* fin);

  // Number of profiles in the fully assembled library.
  int num_profiles_;
  // Number of columns in each context profile.
  int num_cols_;
  // Number of EM-clustering iterations performed on this profile library.
  int iterations_;
  // Context profiles ordered by index.
  std::vector< shared_ptr< ContextProfile<Alphabet> > > profiles_;
  // Flag indicating if profile probabilities are in log- or linspace
  bool logspace_;
};  // ProfileLibrary


template<class Alphabet>
class SamplingProfileInitializer : public ProfileInitializer<Alphabet> {
 public:
  typedef typename
  std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;

  SamplingProfileInitializer(profile_vector profiles,
                             const Pseudocounts<Alphabet>* pc = NULL,
                             float pc_admixture = 1.0f)
      : profiles_(profiles),
        pc_(pc),
        pc_admixture_(pc_admixture) {
    random_shuffle(profiles_.begin(), profiles_.end());
  }
  virtual ~SamplingProfileInitializer() {}
  virtual void Init(ProfileLibrary<Alphabet>& lib) const;

 private:
  // Pool of profile windows to sample from.
  profile_vector profiles_;
  // Pseudocount factory for library profiles.
  const Pseudocounts<Alphabet>* pc_;
  // Constant pseudocount admixture for library profiles.
  float pc_admixture_;
};  // SamplingProfileInitializer

}  // namespace cs

#endif  // SRC_PROFILE_LIBRARY_H_
