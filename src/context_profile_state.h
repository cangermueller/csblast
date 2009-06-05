// Copyright 2009, Andreas Biegert

#ifndef SRC_CONTEXT_PROFILE_STATE_H_
#define SRC_CONTEXT_PROFILE_STATE_H_

#include <cstdio>

#include <functional>
#include <google/sparsetable>

#include "context_profile.h"
#include "profile.h"
#include "transition.h"

using google::sparsetable;

namespace cs {

// Forward declarations
template< class Alphabet, template<class> class State >
class FactorGraph;

// A class representing a context state in a context HMM.
template<class Alphabet>
class ContextProfileState : public ContextProfile<Alphabet> {
 public:
  typedef sparsetable<AnchoredTransition> TransitionTable;
  typedef typename TransitionTable::const_nonempty_iterator ConstTransitionIter;

  // Needed to access names in templatized Profile base class
  using ContextProfile<Alphabet>::num_cols;
  using ContextProfile<Alphabet>::alphabet_size;
  using ContextProfile<Alphabet>::logspace;
  using ContextProfile<Alphabet>::prior;
  using ContextProfile<Alphabet>::set_prior;
  using ContextProfile<Alphabet>::index;
  using ContextProfile<Alphabet>::set_index;
  using ContextProfile<Alphabet>::Read;

  // Constructs HMM state from serialized state read from input stream.
  explicit ContextProfileState(FILE* fin);
  // Constructs HMM state with given profile and all transitions initialized to
  // zero.
  ContextProfileState(int index, int num_states, const Profile<Alphabet>& profile);

  virtual ~ContextProfileState() {}

  // Returns number of in-transitions.
  int num_in_transitions() const { return in_transitions_.num_nonempty(); }
  // Returns number of out-transitions.
  int num_out_transitions() const { return out_transitions_.num_nonempty(); }
  // Clears all in- and out-transitions.
  void ClearTransitions();
  // Returns a const iterator to start of list with non-null in-transition
  // pointers.
  ConstTransitionIter in_transitions_begin() const
  { return in_transitions_.nonempty_begin(); }
  // Returns a const iterator past the end of list with non-null
  // in-transition pointers.
  ConstTransitionIter in_transitions_end() const
  { return in_transitions_.nonempty_end(); }
  // Returns a const iterator to start of list with non-null out-transition
  // pointers.
  ConstTransitionIter out_transitions_begin() const
  { return out_transitions_.nonempty_begin(); }
  // Returns a const iterator past the end of list with non-null out-transition
  // pointers.
  ConstTransitionIter out_transitions_end() const
  { return out_transitions_.nonempty_end(); }

 protected:
  // HMM needs access to transition tables.
  friend class FactorGraph<Alphabet, ::cs::ContextProfileState>;

  // Needed to access names in templatized Profile base class
  using ContextProfile<Alphabet>::kBufferSize;
  using ContextProfile<Alphabet>::ReadHeader;
  using ContextProfile<Alphabet>::ReadBody;
  using ContextProfile<Alphabet>::WriteHeader;
  using ContextProfile<Alphabet>::WriteBody;
  using ContextProfile<Alphabet>::index_;

  // Reads and initializes serialized scalar data members from stream.
  virtual void ReadHeader(FILE* in);
  // Writes serialized scalar data members to stream.
  virtual void WriteHeader(FILE* fout) const;

 private:
  // Class identifier
  static const char* kClassID;

  // Returns serialization class identity.
  virtual const char* class_id() const { return kClassID; }

  // List of in-transitions.
  sparsetable<AnchoredTransition> in_transitions_;
  // List of out-transitions.
  sparsetable<AnchoredTransition> out_transitions_;
};  // class ContextProfileState

}  // namespace cs

#endif  // SRC_CONTEXT_PROFILE_STATE_H_
