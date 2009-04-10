// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_H_
#define SRC_HMM_H_

#include <cstdlib>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "profile.h"
#include "context_profile.h"
#include "count_profile.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"
#include "state.h"
#include "transition.h"

namespace cs {

// Forward declarations
template<class Alphabet>
class HMM;

template<class Alphabet>
class TransitionAdaptor;

template<class Alphabet>
class StateInitializer {
 public:
  StateInitializer() { }
  virtual ~StateInitializer() { }
  virtual void init(HMM<Alphabet>& hmm) const = 0;
};

template<class Alphabet>
class TransitionInitializer {
 public:
  TransitionInitializer() { }
  virtual ~TransitionInitializer() { }
  virtual void init(HMM<Alphabet>& hmm) const = 0;
};


// A hidden Markov model that stores context information in the form of
// context states and state transition probabilities.
template<class Alphabet>
class HMM {
 public:
  // Public typedefs
  typedef std::vector< shared_ptr< State<Alphabet> > > state_vector;
  typedef sparse_matrix<Transition> transition_matrix;
  typedef typename state_vector::iterator state_iterator;
  typedef typename state_vector::const_iterator const_state_iterator;
  typedef typename transition_matrix::nonempty_iterator transition_iterator;
  typedef typename transition_matrix::const_nonempty_iterator const_transition_iterator;

  // Constructs an empty HMM of given size without any states or transitions.
  HMM(int num_states, int num_cols);
  // Constructs context HMM from serialized HMM read from input stream.
  explicit HMM(FILE* fin);
  // Constructs context HMM with the help of a state- and a
  // transition-initializer.
  HMM(int num_states,
      int num_cols,
      const StateInitializer<Alphabet>& st_init,
      const TransitionInitializer<Alphabet>& tr_init);

  virtual ~HMM() {}

  // Initializes HMM states with the provided initializer.
  void init_states(const StateInitializer<Alphabet>& st_init);
  // Initializes HMM transitions with the provided initializer.
  void init_transitions(const TransitionInitializer<Alphabet>& tr_init);
  // Returns true if all states have been fully assembled.
  bool full() const { return static_cast<int>(states_.size()) == num_states_; }
  // Returns the number of states in the HMM
  int num_states() const { return num_states_; }
  // Returns the number of states in the HMM
  int size() const { return num_states_; }
  // Returns the number of columns in each context state.
  int num_cols() const { return num_cols_; }
  // Returns index of central profile column in states.
  int center() const { return (num_cols() - 1) / 2; }
  // Returns the size of the alphabet of the HMM.
  int alphabet_size() const { return Alphabet::instance().size(); }
  // Returns the number of training iterations.
  int iterations() const { return iterations_; }
  // Returns the number of non-null transitions in the HMM.
  int num_transitions() const { return transitions_.num_nonempty(); }
  // Returns the average state connectivity.
  float connectivity() const {
    return static_cast<float>(num_transitions()) / num_states();
  }
  // Accessor methods for state i, where i is from interval [0,num_states].
  State<Alphabet>& operator[](int i) { return *states_[i]; }
  const State<Alphabet>& operator[](int i) const { return *states_[i]; }
  // Accessor methods for transition probability (k,l)
  TransitionAdaptor<Alphabet> operator() (int k, int l) {
    return TransitionAdaptor<Alphabet>(this, k, l);
  }
  float operator() (int k, int l) const {
    return transitions_.get(k,l).probability;
  }
  // Sets the transition probability from state k to state l.
  void set_transition(int k, int l, float prob);
  // Returns the transition probability from state k to state l.
  float transition_probability(int k, int l) const {
    return transitions_.get(k,l).probability;
  }
  // Removes the transition between state k and state l from the HMM.
  void erase_transition(int k, int l);
  // Returns true if there is a transition between state k and state l.
  bool test_transition(int k, int l) const { return transitions_.test(k,l); }
  // Clears all states and transitions.
  void clear();
  // Clears all transitions but leaves profile of states untouched.
  void clear_transitions();
  // Adds the given profile as state to the HMM and returns its state index.
  int add_profile(const Profile<Alphabet>& profile);
  // Returns an iterator to a list of pointers of states.
  state_iterator states_begin() { return states_.begin(); }
  // Returns an iterator pointing past the end of a list of pointers of states.
  state_iterator states_end() { return states_.end(); }
  // Returns a const iterator to a list of pointers of states.
  const_state_iterator states_begin() const { return states_.begin(); }
  // Returns a const iterator pointing past the end of a list of pointers of
  // states.
  const_state_iterator states_end() const { return states_.end(); }
  // Returns an iterator to a list of transitions.
  transition_iterator transitions_begin() {
    return transitions_.nonempty_begin();
  }
  // Returns an iterator pointing past the end of a list of transitions.
  transition_iterator transitions_end() { return transitions_.nonempty_end(); }
  // Returns a const iterator to a list of transitions.
  const_transition_iterator transitions_begin() const {
    return transitions_.nonempty_begin();
  }
  // Returns a const iterator pointing past the end of a list of transitions.
  const_transition_iterator transitions_end() const {
    return transitions_.nonempty_end();
  }
  // Writes the HMM in serialization format to output stream.
  void write(FILE* fout) const;
  // Returns true if transitions are in logspace.
  bool transitions_logspace() const { return transitions_logspace_; }
  // Returns true if state profiles are in logspace.
  bool states_logspace() const { return states_logspace_; }
  // Transforms transitions to logspace.
  void transform_transitions_to_logspace();
  // Transforms transitions to linspace.
  void transform_transitions_to_linspace();
  // Transforms state profiles to logspace.
  void transform_states_to_logspace();
  // Transforms state profiles to linspace.
  void transform_states_to_linspace();
  // Increments the training iteration counter.
  HMM& operator++() {
    ++iterations_;
    return *this;
  }

  // Prints HMM in human-readable format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const HMM& hmm) {
    hmm.print(out);
    return out;
  }

 private:
  // Scaling factor for serialization of profile log values
  static const int kLogScale = 1000;
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Prints the HMM in human-readable format to output stream.
  void print(std::ostream& out) const;
  // Initializes the HMM from a serialized HMM read from stream.
  void read(FILE* fin);
  // Initializes the HMM.
  void init();

  // Number states in the fully assembled HMM
  int num_states_;
  // Number of columns in each context state.
  int num_cols_;
  // Number of training iterations performed on this HMM.
  int iterations_;
  // HMM states ordered by index.
  std::vector< shared_ptr< State<Alphabet> > > states_;
  // Sparse matrix with state transitions.
  sparse_matrix<Transition> transitions_;
  // Flag indicating if HMM transitions are log- or linspace
  bool transitions_logspace_;
  // Flag indicating if HMM profile probabilities are in log- or linspace
  bool states_logspace_;
};  // HMM

template<class Alphabet>
class TransitionAdaptor {
 public:
  TransitionAdaptor(HMM<Alphabet>* hmm, int k, int l)
      : hmm_(hmm), k_(k), l_(l)
  { }

  TransitionAdaptor& operator= (float val) {
    hmm_->set_transition(k_, l_, val);
    return *this;
  }

  operator float() { return hmm_->transition_probability(k_, l_); }

 private:
  HMM<Alphabet>* hmm_;
  int k_;
  int l_;
};

template<class Alphabet>
class SamplingStateInitializer : public StateInitializer<Alphabet> {
 public:
  typedef typename
  std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;

  SamplingStateInitializer(profile_vector profiles,
                           float sample_rate,
                           const Pseudocounts<Alphabet>* pc = NULL,
                           float pc_admixture = 1.0f)
      : profiles_(profiles),
        sample_rate_(sample_rate),
        pc_(pc),
        pc_admixture_(pc_admixture) {
    random_shuffle(profiles_.begin(), profiles_.end());
  }

  virtual ~SamplingStateInitializer() { };
  virtual void init(HMM<Alphabet>& hmm) const;

 private:
  // Pool of full length sequence profiles to sample from.
  profile_vector profiles_;
  // Fraction of profile windows sampled from each subject.
  float sample_rate_;
  // Pseudocount factory for state profiles.
  const Pseudocounts<Alphabet>* pc_;
  // Constant pseudocount admixture for state profiles.
  float pc_admixture_;
};

template<class Alphabet>
class HomogeneousTransitionInitializer : public TransitionInitializer<Alphabet> {
 public:
  HomogeneousTransitionInitializer() { }
  virtual ~HomogeneousTransitionInitializer() { }

  virtual void init(HMM<Alphabet>& hmm) const {
    float prob = 1.0f / hmm.num_states();
    for (int k = 0; k < hmm.num_states(); ++k) {
      for (int l = 0; l < hmm.num_states(); ++l) {
        hmm(k,l) = prob;
      }
    }
  }
};

template<class Alphabet>
class RandomTransitionInitializer : public TransitionInitializer<Alphabet> {
 public:
  RandomTransitionInitializer() { }
  virtual ~RandomTransitionInitializer() { }

  virtual void init(HMM<Alphabet>& hmm) const {
    srand(static_cast<unsigned int>(clock()));
    for (int k = 0; k < hmm.num_states(); ++k)
      for (int l = 0; l < hmm.num_states(); ++l)
        hmm(k,l) =
          static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
    normalize_transitions(hmm);
  }
};

}  // namespace cs

#endif  // SRC_HMM_H_
