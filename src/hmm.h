// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_H_
#define SRC_HMM_H_

#include <cstdlib>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "co_emission-inl.h"
#include "context_profile.h"
#include "count_profile.h"
#include "initializer-inl.h"
#include "profile.h"
#include "profile_library.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"
#include "substitution_matrix.h"
#include "hmm_state.h"
#include "transition.h"

namespace cs {

// Forward declarations
template<class Alphabet>
class HMMTransitionAdaptor;


// A Hidden Markov Model that stores context information in the form of
// context states and inter state transition probabilities.
template<class Alphabet>
class HMM {
 public:
  // Public typedefs
  typedef std::vector< shared_ptr< HMMState<Alphabet> > > StateVec;
  typedef sparse_matrix<Transition> TransitionMatrix;
  typedef typename StateVec::iterator StateIter;
  typedef typename StateVec::const_iterator ConstStateIter;
  typedef typename TransitionMatrix::nonempty_iterator TransitionIter;
  typedef typename TransitionMatrix::const_nonempty_iterator ConstTransitionIter;

  // Constructs an empty HMM of given size without any states or transitions.
  HMM(int num_states, int num_cols);
  // Constructs context HMM from serialized HMM read from input stream.
  explicit HMM(FILE* fin);
  // Constructs context HMM with the help of a state- and a transition-
  // initializer. State profiles are initially set to lin-space.
  HMM(int num_states,
      int num_cols,
      const StateInitializer<Alphabet, ::cs::HMM >& st_init,
      const TransitionInitializer<Alphabet, ::cs::HMM >& tr_init);

  virtual ~HMM() {}

  // Initializes HMM states with the provided initializer.
  void InitStates(const StateInitializer<Alphabet, ::cs::HMM>& st_init);
  // Initializes HMM transitions with the provided initializer.
  void InitTransitions(const TransitionInitializer<Alphabet, ::cs::HMM>& tr_init);
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
  // Returns the mean state connectivity.
  float connectivity() const {
    return static_cast<float>(num_transitions()) / num_states();
  }
  // Accessor methods for state i, where i is from interval [0,num_states].
  HMMState<Alphabet>& operator[](int i) { return *states_[i]; }
  const HMMState<Alphabet>& operator[](int i) const { return *states_[i]; }
  HMMState<Alphabet>& st(int i) { return *states_[i]; }
  const HMMState<Alphabet>& st(int i) const { return *states_[i]; }
  // Accessor methods for transition probability (k,l)
  HMMTransitionAdaptor<Alphabet> operator() (int k, int l) {
    return HMMTransitionAdaptor<Alphabet>(this, k, l);
  }
  float operator() (int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  HMMTransitionAdaptor<Alphabet> tr(int k, int l) {
    return HMMTransitionAdaptor<Alphabet>(this, k, l);
  }
  float tr(int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  // Sets the transition probability from state k to state l.
  void set_transition(int k, int l, float prob);
  // Removes the transition between state k and state l from the HMM.
  void erase_transition(int k, int l);
  // Returns true if there is a transition between state k and state l.
  bool test_transition(int k, int l) const { return transitions_.test(k,l); }
  // Clears all states and transitions.
  void Clear();
  // Clears all transitions but leaves profile of states untouched.
  void ClearTransitions();
  // Adds a profile as state to the HMM and returns its state index.
  int AddState(const Profile<Alphabet>& profile);
  // Adds a context profile as state to the HMM and returns its state index.
  // The prior probability of the context profile becomes the prior probability
  // of the state.
  int AddState(const ContextProfile<Alphabet>& profile);
  // Returns true if transitions are in logspace.
  bool transitions_logspace() const { return transitions_logspace_; }
  // Returns true if state profiles are in logspace.
  bool states_logspace() const { return states_logspace_; }
  // Increments the training iteration counter.
  void increment_iterations() { ++iterations_; }
  // Transforms transitions to logspace.
  void TransformTransitionsToLogSpace();
  // Transforms transitions to linspace.
  void TransformTransitionsToLinSpace();
  // Transforms state profiles to logspace.
  void TransformStatesToLogSpace();
  // Transforms state profiles to linspace.
  void TransformStatesToLinSpace();
  // Writes the HMM in serialization format to output stream.
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

  // Prints HMM in human-readable format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const HMM& hmm) {
    hmm.Print(out);
    return out;
  }

 private:
  // Scaling factor for serialization of profile log values
  static const int kLogScale = 1000;
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Prints the HMM in human-readable format to output stream.
  void Print(std::ostream& out) const;
  // Initializes the HMM from a serialized HMM read from stream.
  void Read(FILE* fin);
  // Initializes the HMM.
  void Init();

  // Number states in the fully assembled HMM
  int num_states_;
  // Number of columns in each context state.
  int num_cols_;
  // Number of training iterations performed on this HMM.
  int iterations_;
  // HMM states ordered by index.
  std::vector< shared_ptr< HMMState<Alphabet> > > states_;
  // Sparse matrix with state transitions.
  sparse_matrix<Transition> transitions_;
  // Flag indicating if HMM transitions are log- or linspace
  bool transitions_logspace_;
  // Flag indicating if HMM profile probabilities are in log- or linspace
  bool states_logspace_;

  friend class HMMTransitionAdaptor<Alphabet>;
  friend class LibraryStateInitializer<Alphabet, ::cs::HMM>;
};  // class HMM

template<class Alphabet>
class HMMTransitionAdaptor {
 public:
  HMMTransitionAdaptor(HMM<Alphabet>* hmm, int k, int l)
      : hmm_(hmm), k_(k), l_(l) {}
  HMMTransitionAdaptor& operator= (float val) {
    hmm_->set_transition(k_, l_, val);
    return *this;
  }
  operator float() { return hmm_->transitions_.get(k_, l_).weight; }
  Transition* operator& () { return &hmm_->transitions_.mutating_get(k_, l_); }

 private:
  HMM<Alphabet>* hmm_;
  int k_;
  int l_;
};  // class HMMTransitionAdaptor

// template<class Alphabet>
// class SamplingHMMStateInitializer : public HMMStateInitializer<Alphabet> {
//  public:
//   typedef typename
//   std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
//   typedef typename profile_vector::const_iterator profile_iterator;

//   SamplingHMMStateInitializer(profile_vector profiles,
//                               float sample_rate,
//                               const Pseudocounts<Alphabet>* pc = NULL,
//                               float pc_admixture = 1.0f)
//       : profiles_(profiles),
//         sample_rate_(sample_rate),
//         pc_(pc),
//         pc_admixture_(pc_admixture) {
//     random_shuffle(profiles_.begin(), profiles_.end());
//   }

//   virtual ~SamplingHMMStateInitializer() {};
//   virtual void Init(HMM<Alphabet>& hmm) const;

//  private:
//   // Pool of full length sequence profiles to sample from.
//   profile_vector profiles_;
//   // Fraction of profile windows sampled from each subject.
//   float sample_rate_;
//   // Pseudocount factory for state profiles.
//   const Pseudocounts<Alphabet>* pc_;
//   // Constant pseudocount admixture for state profiles.
//   float pc_admixture_;
// };  // class SamplingHMMStateInitializer

// // Compare function to sort states in descending prior probability.
// template<class Alphabet>
// bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
//                   const shared_ptr< ContextProfile<Alphabet> >& rhs);

// template<class Alphabet>
// class LibraryHMMStateInitializer : public HMMStateInitializer<Alphabet> {
//  public:
//   LibraryHMMStateInitializer(const ProfileLibrary<Alphabet>* lib)
//       : lib_(lib) {}

//   virtual ~LibraryHMMStateInitializer() {};
//   virtual void Init(HMM<Alphabet>& hmm) const;

//  private:
//   // Profile library of context profiles.
//   const ProfileLibrary<Alphabet>* lib_;
// };  // class LibraryHMMStateInitializer

// template<class Alphabet>
// class HomogeneousHMMTransitionInitializer
//     : public HMMTransitionInitializer<Alphabet> {
//  public:
//   HomogeneousHMMTransitionInitializer() {}
//   virtual ~HomogeneousHMMTransitionInitializer() {}

//   virtual void Init(HMM<Alphabet>& hmm) const {
//     float prob = 1.0f / hmm.num_states();
//     for (int k = 0; k < hmm.num_states(); ++k) {
//       for (int l = 0; l < hmm.num_states(); ++l) {
//         hmm(k,l) = prob;
//       }
//     }
//   }
// };  // class HomogeneousHMMTransitionInitializer

// template<class Alphabet>
// class RandomHMMTransitionInitializer : public HMMTransitionInitializer<Alphabet> {
//  public:
//   RandomHMMTransitionInitializer() {}
//   virtual ~RandomHMMTransitionInitializer() {}

//   virtual void Init(HMM<Alphabet>& hmm) const {
//     srand(static_cast<unsigned int>(clock()));
//     for (int k = 0; k < hmm.num_states(); ++k)
//       for (int l = 0; l < hmm.num_states(); ++l)
//         hmm(k,l) =
//           static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
//     NormalizeTransitions(hmm);
//   }
// };  // class RandomHMMTransitionInitializer

// template<class Alphabet>
// class CoEmissionHMMTransitionInitializer
//     : public HMMTransitionInitializer<Alphabet> {
//  public:
//   CoEmissionHMMTransitionInitializer(const SubstitutionMatrix<Alphabet>* sm,
//                                   float score_thresh)
//       : co_emission_(sm), score_thresh_(score_thresh) {}
//   virtual ~CoEmissionHMMTransitionInitializer() {}

//   virtual void Init(HMM<Alphabet>& hmm) const;

//  private:
//   // Function object for calculation of co-emission scores
//   CoEmission<Alphabet> co_emission_;
//   // Minimal co-emission score for inclusion in transition set
//   float score_thresh_;
// }; // class CoEmissionHMMTransitionInitializer

}  // namespace cs

#endif  // SRC_HMM_H_
