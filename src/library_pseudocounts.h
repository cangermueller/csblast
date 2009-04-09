// Copyright 2009, Andreas Biegert

#ifndef SRC_LIBRARY_PSEUDOCOUNTS_H_
#define SRC_LIBRARY_PSEUDOCOUNTS_H_

#include <valarray>

#include "count_profile-inl.h"
#include "emitter.h"
#include "matrix.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from a library
// of context profiles.
template<class Alphabet>
class LibraryPseudocounts : public Pseudocounts<Alphabet> {
 public:
  LibraryPseudocounts(const ProfileLibrary<Alphabet>* lib,
                      const EmissionParams& params);
  ~LibraryPseudocounts() { }

  // Adds context-specific pseudocounts to sequence and stores resulting
  // frequencies in given profile.
  virtual void add_to_sequence(const Sequence<Alphabet>& seq,
                               const Admixture& pca,
                               Profile<Alphabet>* profile) const;
  // Adds context-specific pseudocounts to alignment derived profile.
  virtual void add_to_profile(const Admixture& pca,
                              CountProfile<Alphabet>* p) const;

 private:
  // Profile library with context profiles.
  const ProfileLibrary<Alphabet>& lib_;
  // Needed to compute emission probabilities of context profiles.
  const Emitter<Alphabet> emitter_;

  DISALLOW_COPY_AND_ASSIGN(LibraryPseudocounts);
};  // LibraryPseudocounts

}  // namespace cs

#endif  // SRC_LIBRARY_PSEUDOCOUNTS_H_
