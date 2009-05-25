// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_STATE_H_
#define SRC_HMM_STATE_H_

#include <cstdio>

#include <functional>
#include <google/sparsetable>

#include "context_profile.h"
#include "profile.h"
#include "transition.h"

using google::sparsetable;

namespace cs {

// Forward declarations
template<class Alphabet>
class HMM;

// A class representing a context state in a context HMM.
template<class Alphabet>
class HMMState : public ContextProfile<Alphabet> {
 public:
  typedef typename
  sparsetable<AnchoredTransition>::const_nonempty_iterator ConstTransitionIter;

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
  explicit HMMState(FILE* fin);
  // Constructs HMM state with given profile and all transitions initialized to
  // zero.
  HMMState(int index, const Profile<Alphabet>& profile, int num_states);
  // Constructs HMM state with given context profile and all transitions
  // initialized to zero.
  HMMState(int index, const ContextProfile<Alphabet>& profile, int num_states);

  virtual ~HMMState() {}

  // Returns number of in-transitions.
  int num_in_transitions() const
  { return in_transitions_.num_nonempty(); }
  // Returns number of out-transitions.
  int num_out_transitions() const
  { return out_transitions_.num_nonempty(); }
  // Returns the transition probability from this state to state k.
  float to(int k) const;
  // Returns the transition probability from state k to this state.
  float from(int k) const;
  // Clears all in- and out-transitions.
  void clear_transitions();
  // Resizes the transition tables to new HMM size.
  void resize(int num_states);

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
  friend class HMM<Alphabet>;

  // Needed to access names in templatized Profile base class
  using ContextProfile<Alphabet>::kBufferSize;
  using ContextProfile<Alphabet>::read_header;
  using ContextProfile<Alphabet>::read_body;
  using ContextProfile<Alphabet>::write_header;
  using ContextProfile<Alphabet>::write_body;
  using ContextProfile<Alphabet>::index_;

  // Reads and initializes serialized scalar data members from stream.
  virtual void read_header(FILE* in);
  // Writes serialized scalar data members to stream.
  virtual void write_header(FILE* fout) const;

 private:
  // Class identifier
  static const char* kClassID;

  // Returns serialization class identity.
  virtual const char* class_id() const { return kClassID; }

  // Size of HMM to which the state belongs.
  int num_states_;
  // List of in-transitions.
  sparsetable<AnchoredTransition> in_transitions_;
  // List of out-transitions.
  sparsetable<AnchoredTransition> out_transitions_;
};  // HMMState

// Comparison functor that compares the index of two states.
template<class Alphabet>
struct HMMStateIndexCompare : public std::binary_function< HMMState<Alphabet>,
                                                        HMMState<Alphabet>,
                                                        bool > {
  bool operator()(const HMMState<Alphabet>& lhs, const HMMState<Alphabet>& rhs) const {
    return lhs.index() < rhs.index();
  }
};

// Comparison functor that compares the number of in-transitions of two states.
template<class Alphabet>
struct NumInTransitionsCompare : public std::binary_function< HMMState<Alphabet>,
                                                              HMMState<Alphabet>,
                                                              bool > {
  bool operator()(const HMMState<Alphabet>& lhs, const HMMState<Alphabet>& rhs) const {
    return lhs.num_in_transitions() < rhs.num_in_transitions();
  }
};

// Comparison functor that compares the number of out-transitions of two states.
template<class Alphabet>
struct NumOutTransitionsCompare : public std::binary_function< HMMState<Alphabet>,
                                                               HMMState<Alphabet>,
                                                               bool > {
  bool operator()(const HMMState<Alphabet>& lhs, const HMMState<Alphabet>& rhs) const {
    return lhs.num_out_transitions() < rhs.num_out_transitions();
  }
};

}  // namespace cs

#endif  // SRC_HMM_STATE_H_
