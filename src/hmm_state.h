// Copyright 2009, Andreas Biegert

#ifndef CS_HMM_STATE_H_
#define CS_HMM_STATE_H_

#include <google/sparsetable>

#include "context_profile-inl.h"

namespace cs {

// Forward declarations
class Transition;

// Simple container class for an HMM state deriving from a context-profile with
// additional arrays containing pointers to in- and out-transitions. The transitions
// themselves are stored in the HMM object to which the state belongs to.
template<class Abc>
struct HmmState : public ContextProfile<Abc> {
  typedef google::sparsetable<Transition*> TransitionTable;

  // Default construction
  HmmState() : ContextProfile<Abc>() {}

  // Constructs a context profile with 'len' columns
  explicit HmmState(size_t len) ContextProfile<Abc>(len) {}

  // Construction from serialized profile read from input stream.
  explicit HmmState(FILE* fin) : ContextProfile<Abc>(fin) {}

  // Constructs a context profile from normalized values in given profile 'p'.
  explicit HmmState(const Profile<Abc>& p) : ContextProfile<Abc>(p) {}

  // Construct a context profile by copying subprofile starting at index 'start'
  // for 'len' columns and normalizing afterwards.
  HmmState(const Profile<Abc>& p, size_t start, size_t len)
      : ContextProfile<Abc>(p, start, len) {}

  TransitionTable in_tr;   // sparse table with pointers to in-transitions
  TransitionTable out_tr;  // sparse table with pointers to out-transitions
};

}  // namespace cs

#endif  // CS_HMM_STATE_H_
