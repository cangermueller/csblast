// Copyright 2009, Andreas Biegert

#ifndef SRC_FACTOR_GRAPH_H_
#define SRC_FACTOR_GRAPH_H_

#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <vector>

#include "globals.h"
#include "co_emission.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "pseudocounts.h"
#include "profile.h"
#include "shared_ptr.h"
#include "substitution_matrix.h"
#include "sparse_matrix.h"
#include "transition.h"

namespace cs {

// Forward declarations
template< class Alphabet, template<class> class State >
class TransitionAdaptor;

template< class Alphabet, template<class> class State >
class StateInitializer;

template< class Alphabet, template<class> class State >
class TransitionInitializer;


// Abstract base class of a linear chain factor graph from which classes CRF and
// HMM derive.
template< class Alphabet, template<class> class State >
class FactorGraph {
 public:
  // Public typedefs
  typedef shared_ptr< State<Alphabet> > StatePtr;
  typedef std::vector<StatePtr> StateVec;
  typedef sparse_matrix<Transition> TransitionMatrix;
  typedef typename StateVec::iterator StateIter;
  typedef typename StateVec::const_iterator ConstStateIter;
  typedef typename TransitionMatrix::nonempty_iterator TransitionIter;
  typedef typename TransitionMatrix::const_nonempty_iterator ConstTransitionIter;

  // Constructs a dummy factor graph
  FactorGraph();
  // Constructs an empty factor graph without any states or transitions.
  FactorGraph(int num_states, int num_cols);
  // Constructs a factor graph from serialized graph read from input stream.
  explicit FactorGraph(FILE* fin);

  virtual ~FactorGraph() {}

  // Adds a states to the factor graph initialized with given profile and
  // returns the index of the new state. Caller has to guarantee that the profile
  // is in lin-space
  virtual int AddState(const Profile<Alphabet>& profile) = 0;
  // Initializes states with the provided initializer.
  void InitStates(const StateInitializer<Alphabet, State>& st);
  // Initializes factor graph transitions with the provided initializer.
  void InitTransitions(const TransitionInitializer<Alphabet, State>& tr);
  // Returns true if all states have been fully assembled.
  bool full() const { return static_cast<int>(states_.size()) == num_states_; }
  // Returns the number of states of the full graph.
  int num_states() const { return num_states_; }
  // Returns the number of states of the full graph.
  int size() const { return num_states_; }
  // Returns the number of context columns.
  int num_cols() const { return num_cols_; }
  // Returns index of central context column.
  int center() const { return (num_cols() - 1) / 2; }
  // Returns the size of the alphabet of the factor graph.
  int alphabet_size() const { return Alphabet::instance().size(); }
  // Returns the number of training iterations.
  int iterations() const { return iterations_; }
  // Returns the number of transitions.
  int num_transitions() const { return transitions_.num_nonempty(); }
  // Returns the mean connectivity.
  float connectivity() const {
    return static_cast<float>(num_transitions()) / num_states();
  }
  // Accessor methods for state i, where i is from interval [0,num_states].
  State<Alphabet>& operator[](int i) { return *states_[i]; }
  const State<Alphabet>& operator[](int i) const { return *states_[i]; }
  // Accessor methods for transition probability (k,l)
  TransitionAdaptor<Alphabet, State> operator() (int k, int l) {
    return TransitionAdaptor<Alphabet, State>(this, k, l);
  }
  float operator() (int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  TransitionAdaptor<Alphabet, State> tr(int k, int l) {
    return TransitionAdaptor<Alphabet, State>(this, k, l);
  }
  float tr(int k, int l) const {
    return transitions_.get(k,l).weight;
  }
  // Sets the transition from state k to state l to value w.
  void set_transition(int k, int l, float w);
  // Removes the transition between state k and state l from the factor graph.
  void erase_transition(int k, int l);
  // Returns true if there is a transition between state k and state l.
  bool test_transition(int k, int l) const { return transitions_.test(k,l); }
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
  // Clears all states and transitions.
  void Clear();
  // Clears all transitions but leaves states untouched.
  void ClearTransitions();
  // Returns true if transitions are in logspace.
  bool transitions_logspace() const { return transitions_logspace_; }
  // Increments the training iteration counter.
  void increment_iterations() { ++iterations_; }
  // Transforms transitions to logspace.
  void TransformTransitionsToLogSpace();
  // Transforms transitions to linspace.
  void TransformTransitionsToLinSpace();
  // Writes the factor graph in serialization format to output stream.
  void Write(FILE* fout) const;

  // Prints factor graph in human-readable format for debugging.
  friend std::ostream& operator<< (std::ostream& out, const FactorGraph& g) {
    g.Print(out);
    return out;
  }

 protected:
  // Scaling factor for serialization of profile log values
  static const int kLogScale = 1000;
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Initializes the factor graph.
  void Init();
  // Initializes the factor graph from a serialized factor graph read from stream.
  void Read(FILE* fin);
  // Writes header section for serialization.
  virtual void ReadHeader(FILE* fout);
  // Writes states in serialization format.
  virtual void ReadStates(FILE* fout);
  // Writes transitions in serialization format.
  virtual void ReadTransitions(FILE* fout);
  // Writes header section for serialization.
  virtual void WriteHeader(FILE* fout) const;
  // Writes states in serialization format.
  virtual void WriteStates(FILE* fout) const;
  // Writes transitions in serialization format.
  virtual void WriteTransitions(FILE* fout) const;
  // Prints the factor graph in human-readable format to output stream.
  virtual void Print(std::ostream& out) const;

  // Number states in the fully assembled factor graph
  int num_states_;
  // Number of columns in each state.
  int num_cols_;
  // Number of training iterations performed on this factor graph.
  int iterations_;
  // Factor graph states ordered by index.
  StateVec states_;
  // Sparse matrix with state transitions.
  sparse_matrix<Transition> transitions_;
  // Flag indicating if HMM transitions are log- or linspace
  bool transitions_logspace_;

 private:
  // Returns serialization class identity.
  virtual const char* class_id() const = 0;

  friend class TransitionAdaptor<Alphabet, State>;
};  // class FactorGraph


// Opaque object that acts like a transition.
template< class Alphabet, template<class> class State >
class TransitionAdaptor {
 public:
  TransitionAdaptor(FactorGraph<Alphabet, State>* g, int k, int l)
      : g_(g), k_(k), l_(l) {}
  TransitionAdaptor& operator= (float val) {
    g_->set_transition(k_, l_, val);
    return *this;
  }
  operator float() { return g_->transitions_.get(k_, l_).weight; }
  Transition* operator& () { return &g_->transitions_.mutating_get(k_, l_); }

 private:
  FactorGraph<Alphabet, State>* g_;
  int k_;
  int l_;
};  // class TransitionAdapter


// Normalizes transition probabilities to one.
template< class Alphabet, template<class> class State >
void NormalizeTransitions(FactorGraph<Alphabet, State>* graph, float f = 1.0f);

// Compare function to sort profiles in order of descending prior probability.
template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs);


// Abstract base class for state initialization strategies.
template< class Alphabet, template<class> class State >
class StateInitializer {
 public:
  StateInitializer() {}
  virtual ~StateInitializer() {}
  virtual void Init(FactorGraph<Alphabet, State>& graph) const = 0;
};  // class StateInitializer

// Strategy that randomly samples states from training profiles.
template< class Alphabet, template<class> class State >
class SamplingStateInitializer : public StateInitializer<Alphabet, State> {
 public:
  typedef shared_ptr< CountProfile<Alphabet> > ProfilePtr;
  typedef typename std::vector<ProfilePtr> ProfileVec;
  typedef typename std::vector<ProfilePtr>::const_iterator ProfileIter;

  SamplingStateInitializer(ProfileVec profiles,
                           float sample_rate,
                           const Pseudocounts<Alphabet>* pc = NULL,
                           float pc_admixture = 1.0f);
  virtual ~SamplingStateInitializer() {};
  virtual void Init(FactorGraph<Alphabet, State>& graph) const;

 private:
  // Pool of full length sequence profiles to sample from.
  ProfileVec profiles_;
  // Fraction of profile windows sampled from each subject.
  float sample_rate_;
  // Pseudocount factory for state profiles.
  const Pseudocounts<Alphabet>* pc_;
  // Constant pseudocount admixture for state profiles.
  float pc_admixture_;
};  // class SamplingStateInitializer

// Strategy that uses context profiles from a profile library to initialize
// states in the factor graph.
template< class Alphabet, template<class> class State >
class LibraryStateInitializer : public StateInitializer<Alphabet, State> {
 public:
  LibraryStateInitializer(const ProfileLibrary<Alphabet>* lib);
  virtual ~LibraryStateInitializer() {};
  virtual void Init(FactorGraph<Alphabet, State>& graph) const;

 private:
  typedef shared_ptr< ContextProfile<Alphabet> > ContextProfilePtr;
  typedef std::vector<ContextProfilePtr> ContextProfileVec;
  typedef typename ContextProfileVec::const_iterator ContextProfileIter;

  // Profile library of context profiles.
  const ProfileLibrary<Alphabet>* lib_;
};  // class LibraryStateInitializer


// Abstract base class for transition initialization strategies.
template< class Alphabet, template<class> class State >
class TransitionInitializer {
 public:
  TransitionInitializer() {}
  virtual ~TransitionInitializer() {}
  virtual void Init(FactorGraph<Alphabet, State>& graph) const = 0;
};  // class TransitionInitializer

// Strategy that initializes transitions homogeneously.
template< class Alphabet, template<class> class State >
class HomogeneousTransitionInitializer
    : public TransitionInitializer<Alphabet, State> {
 public:
  HomogeneousTransitionInitializer() {}
  virtual ~HomogeneousTransitionInitializer() {}
  virtual void Init(FactorGraph<Alphabet, State>& graph) const;
};  // class HomogeneousTransitionInitializer

// Strategy that initializes transition probabilities at random.
template< class Alphabet, template<class> class State >
class RandomTransitionInitializer :
      public TransitionInitializer<Alphabet, State> {
 public:
  RandomTransitionInitializer() {}
  virtual ~RandomTransitionInitializer() {}
  virtual void Init(FactorGraph<Alphabet, State>& graph) const;
};  // class RandomTransitionInitializer

// Strategy that uses the co-emission scores of profile pairs to determine
// transition probabilities. No transition is added if the co-emission is below
// given threshold.
template< class Alphabet, template<class> class State >
class CoEmissionTransitionInitializer
    : public TransitionInitializer<Alphabet, State> {
 public:
  CoEmissionTransitionInitializer(const SubstitutionMatrix<Alphabet>* sm,
                                  float score_thresh);
  virtual ~CoEmissionTransitionInitializer() {}
  virtual void Init(FactorGraph<Alphabet, State>& graph) const;

 private:
  // Function object for calculation of co-emission scores
  CoEmission<Alphabet> co_emission_;
  // Minimal co-emission score for inclusion in transition set
  float score_thresh_;
}; // class CoEmissionTransitionInitializer

}  // namespace cs

#endif  // SRC_FACTOR_GRAPH_H_
