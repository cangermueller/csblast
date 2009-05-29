// Copyright 2009, Andreas Biegert

#ifndef SRC_CRF_H_
#define SRC_CRF_H_

#include <cstdlib>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "co_emission-inl.h"
#include "context_profile.h"
#include "count_profile.h"
#include "profile.h"
#include "profile_library.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"
#include "substitution_matrix.h"
#include "crf_state.h"
#include "transition.h"

namespace cs {

// Forward declarations
template<class Alphabet>
class CRF;

template<class Alphabet>
class CRFTransitionAdaptor;


template<class Alphabet>
class CRFStateInitializer {
 public:
  CRFStateInitializer() {}
  virtual ~CRFStateInitializer() {}
  virtual void Init(CRF<Alphabet>& crf) const = 0;
};  // class CRFStateInitializer

template<class Alphabet>
class CRFTransitionInitializer {
 public:
  CRFTransitionInitializer() {}
  virtual ~CRFTransitionInitializer() {}
  virtual void Init(CRF<Alphabet>& crf) const = 0;
};  // class CRFTransitionInitializer


// A hidden Markov model that stores context information in the form of
// context states and state transition probabilities.
template<class Alphabet>
class CRF {
 public:
  // Public typedefs
  typedef std::vector< shared_ptr< CRFState<Alphabet> > > StateVec;
  typedef sparse_matrix<Transition> TransitionMatrix;
  typedef typename StateVec::iterator StateIter;
  typedef typename StateVec::const_iterator ConstStateIter;
  typedef typename TransitionMatrix::nonempty_iterator TransitionIter;
  typedef typename TransitionMatrix::const_nonempty_iterator ConstTransitionIter;

  // Constructs an empty CRF of given size without any states or transitions.
  CRF(int num_states, int num_cols);
  // Constructs context CRF from serialized CRF read from input stream.
  explicit CRF(FILE* fin);
  // Constructs context CRF with the help of a state- and a transition-
  // initializer.
  CRF(int num_states,
      int num_cols,
      const CRFStateInitializer<Alphabet>& st_init,
      const CRFTransitionInitializer<Alphabet>& tr_init);

  virtual ~CRF() {}

  // Initializes CRF states with the provided initializer.
  void InitStates(const CRFStateInitializer<Alphabet>& st_init);
  // Initializes CRF transitions with the provided initializer.
  void InitTransitions(const CRFTransitionInitializer<Alphabet>& tr_init);
  // Returns true if all states have been fully assembled.
  bool full() const { return static_cast<int>(states_.size()) == num_states_; }
  // Returns the number of states in the CRF
  int num_states() const { return num_states_; }
  // Returns the number of states in the CRF
  int size() const { return num_states_; }
  // Returns the number of columns in each context state.
  int num_cols() const { return num_cols_; }
  // Returns index of central context column in states.
  int center() const { return (num_cols() - 1) / 2; }
  // Returns the size of the alphabet of the CRF.
  int alphabet_size() const { return Alphabet::instance().size(); }
  // Returns the number of optimization iterations.
  int iterations() const { return iterations_; }
  // Returns the number of non-null transitions in the CRF.
  int num_transitions() const { return transitions_.num_nonempty(); }
  // Returns the mean state connectivity.
  float connectivity() const {
    return static_cast<float>(num_transitions()) / num_states();
  }
  // Accessor methods for state i, where i is from interval [0,num_states].
  CRFState<Alphabet>& operator[](int i) { return *states_[i]; }
  const CRFState<Alphabet>& operator[](int i) const { return *states_[i]; }
  CRFState<Alphabet>& st(int i) { return *states_[i]; }
  const CRFState<Alphabet>& st(int i) const { return *states_[i]; }
  // Accessor methods for transition probability (k,l)
  CRFTransitionAdaptor<Alphabet> operator() (int k, int l) {
    return CRFTransitionAdaptor<Alphabet>(this, k, l);
  }
  float operator() (int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  CRFTransitionAdaptor<Alphabet> tr(int k, int l) {
    return CRFTransitionAdaptor<Alphabet>(this, k, l);
  }
  float tr(int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  // Sets the transition from state k to state l to value w.
  void set_transition(int k, int l, float w);
  // Removes the transition between state k and state l from the CRF.
  void erase_transition(int k, int l);
  // Returns true if there is a transition between state k and state l.
  bool test_transition(int k, int l) const { return transitions_.test(k,l); }
  // Clears all states and transitions.
  void Clear();
  // Clears all transitions but leaves states untouched.
  void ClearTransitions();
  // Adds a states to the CRF initialized with given profile and returns the
  // index of the new state.
  int AddState(const Profile<Alphabet>& profile);
  // Returns true if transitions are in logspace.
  bool transitions_logspace() const { return transitions_logspace_; }
  // Increments the training iteration counter.
  void increment_iterations() { ++iterations_; }
  // Transforms transitions to logspace.
  void TransformTransitionsToLogSpace();
  // Transforms transitions to linspace.
  void TransformTransitionsToLinSpace();
  // Writes the CRF in serialization format to output stream.
  void Write(FILE* fout) const;
  // Returns an iterator to a list of pointers of states.
  StateIter states_begin() { return states_.begin(); }
  // Returns an iterator pointing past the end of a list of pointers of states.
  StateIter states_end() { return states_.end(); }
  // Returns a const iterator to a list of pointers of states.
  ConstStateIter states_begin() const { return states_.begin(); }
  // Returns a const iterator pointing past the end of a list of pointers of
  // states.
  ConstStateIter states_end() const { return states_.end(); }
  // Returns an iterator to a list of transitions.
  TransitionIter transitions_begin() {
    return transitions_.nonempty_begin();
  }
  // Returns an iterator pointing past the end of a list of transitions.
  TransitionIter transitions_end() { return transitions_.nonempty_end(); }
  // Returns a const iterator to a list of transitions.
  ConstTransitionIter transitions_begin() const {
    return transitions_.nonempty_begin();
  }
  // Returns a const iterator pointing past the end of a list of transitions.
  ConstTransitionIter transitions_end() const {
    return transitions_.nonempty_end();
  }

  // Prints CRF in human-readable format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const CRF& crf) {
    crf.Print(out);
    return out;
  }

 private:
  // Scaling factor for serialization of profile log values
  static const int kLogScale = 1000;
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Prints the CRF in human-readable format to output stream.
  void Print(std::ostream& out) const;
  // Initializes the CRF from a serialized CRF read from stream.
  void Read(FILE* fin);
  // Initializes the CRF.
  void Init();

  // Number states in the fully assembled CRF
  int num_states_;
  // Number of columns in each context state.
  int num_cols_;
  // Number of training iterations performed on this CRF.
  int iterations_;
  // CRF states ordered by index.
  std::vector< shared_ptr< CRFState<Alphabet> > > states_;
  // Sparse matrix with state transitions.
  sparse_matrix<Transition> transitions_;
  // Flag indicating if HMM transitions are log- or linspace
  bool transitions_logspace_;

  friend class CRFTransitionAdaptor<Alphabet>;
};  // class CRF

template<class Alphabet>
class CRFTransitionAdaptor {
 public:
  CRFTransitionAdaptor(CRF<Alphabet>* crf, int k, int l)
      : crf_(crf), k_(k), l_(l) {}
  CRFTransitionAdaptor& operator= (float val) {
    crf_->set_transition(k_, l_, val);
    return *this;
  }
  operator float() { return crf_->transitions_.get(k_, l_).weight; }
  Transition* operator& () { return &hmm_->transitions_.mutating_get(k_, l_); }

 private:
  CRF<Alphabet>* crf_;
  int k_;
  int l_;
};  // class CRFTransitionAdapter

template<class Alphabet>
class SamplingCRFStateInitializer : public CRFStateInitializer<Alphabet> {
 public:
  typedef typename
  std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;

  SamplingCRFStateInitializer(profile_vector profiles,
                              float sample_rate,
                              const Pseudocounts<Alphabet>* pc = NULL,
                              float pc_admixture = 1.0f)
      : profiles_(profiles),
        sample_rate_(sample_rate),
        pc_(pc),
        pc_admixture_(pc_admixture) {
    random_shuffle(profiles_.begin(), profiles_.end());
  }

  virtual ~SamplingCRFStateInitializer() {};
  virtual void Init(CRF<Alphabet>& crf) const;

 private:
  // Pool of full length sequence profiles to sample from.
  profile_vector profiles_;
  // Fraction of profile windows sampled from each subject.
  float sample_rate_;
  // Pseudocount factory for state profiles.
  const Pseudocounts<Alphabet>* pc_;
  // Constant pseudocount admixture for state profiles.
  float pc_admixture_;
};  // class SamplingCRFStateInitializer

template<class Alphabet>
class RandomCRFTransitionInitializer : public CRFTransitionInitializer<Alphabet> {
 public:
  RandomCRFTransitionInitializer() {}
  virtual ~RandomCRFTransitionInitializer() {}

  virtual void Init(CRF<Alphabet>& crf) const {
    srand(static_cast<unsigned int>(clock()));
    for (int k = 0; k < crf.num_states(); ++k)
      for (int l = 0; l < crf.num_states(); ++l)
        crf(k,l) =
          static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
    NormalizeTransitions(crf);
  }
};  // class RandomCRFTransitionInitializer

template<class Alphabet>
class CoEmissionCRFTransitionInitializer
    : public CRFTransitionInitializer<Alphabet> {
 public:
  CoEmissionCRFTransitionInitializer(const SubstitutionMatrix<Alphabet>* sm,
                                     float score_thresh)
      : co_emission_(sm), score_thresh_(score_thresh) {}
  virtual ~CoEmissionCRFTransitionInitializer() {}

  virtual void Init(CRF<Alphabet>& crf) const;

 private:
  // Function object for calculation of co-emission scores
  CoEmission<Alphabet> co_emission_;
  // Minimal co-emission score for inclusion in transition set
  float score_thresh_;
};  // class CoEmissionCRFTransitionInitializer

}  // namespace cs

#endif  // SRC_CRF_H_
