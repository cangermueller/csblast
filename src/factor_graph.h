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

template<class Alphabet>
class TransitionAdaptor;


// Factor graph base class for factor graphs and HMMs.
template< class Alphabet, template<class> class State >
class FactorGraph {
 public:
  // Public typedefs
  typedef std::vector< shared_ptr<State> > StateVec;
  typedef sparse_matrix<Transition> TransitionMatrix;
  typedef typename StateVec::iterator StateIter;
  typedef typename StateVec::const_iterator ConstStateIter;
  typedef typename TransitionMatrix::nonempty_iterator TransitionIter;
  typedef typename TransitionMatrix::const_nonempty_iterator ConstTransitionIter;

  // Constructs an empty factor graph without any states or transitions.
  explicit FactorGraph(int num_states);
  // Constructs a factor graph from serialized graph read from input stream.
  explicit FactorGraph(FILE* fin);

  virtual ~FactorGraph() {}

  // Adds a states to the factor graph initialized with given profile and
  // returns the index of the new state.
  virtual int AddState(const Profile<Alphabet>& profile) = 0;
  // Initializes states with the provided initializer.
  void InitStates(const StateInitializer<Alphabet>& st_init);
  // Initializes factor graph transitions with the provided initializer.
  void InitTransitions(const TransitionInitializer<Alphabet>& tr_init);
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
  TransitionAdaptor<Alphabet> operator() (int k, int l) {
    return TransitionAdaptor<Alphabet>(this, k, l);
  }
  float operator() (int k, int l) const {
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

 private:
  // Scaling factor for serialization of profile log values
  static const int kLogScale = 1000;
  // Buffer size for reading
  static const int kBufferSize = KB;

  // Initializes the factor graph.
  void Init();
  // Initializes the factor graph from a serialized factor graph read from stream.
  void Read(FILE* fin);
  // Writes header section for serialization.
  virtual ReadHeader(FILE* fout) const;
  // Writes states in serialization format.
  virtual ReadStates(FILE* fout) const;
  // Writes transitions in serialization format.
  virtual ReadTransitions(FILE* fout) const;
  // Writes header section for serialization.
  virtual WriteHeader(FILE* fout) const;
  // Writes states in serialization format.
  virtual WriteStates(FILE* fout) const;
  // Writes transitions in serialization format.
  virtual WriteTransitions(FILE* fout) const;
  // Prints the factor graph in human-readable format to output stream.
  void Print(std::ostream& out) const;
  // Returns serialization class identity.
  virtual const char* class_id() const = 0;

  // Number states in the fully assembled factor graph
  int num_states_;
  // Number of columns in each context state.
  int num_cols_;
  // Number of training iterations performed on this factor graph.
  int iterations_;
  // Factor graph states ordered by index.
  std::vector< shared_ptr< State<Alphabet> > > states_;
  // Sparse matrix with state transitions.
  sparse_matrix<Transition> transitions_;
  // Flag indicating if HMM transitions are log- or linspace
  bool transitions_logspace_;

  friend class TransitionAdaptor<Alphabet>;
};  // class FactorGraph

template<class Alphabet>
class TransitionAdaptor {
 public:
  TransitionAdaptor(FactorGraph<Alphabet>* g, int k, int l)
      : g_(g), k_(k), l_(l) {}
  TransitionAdaptor& operator= (float val) {
    g_->set_transition(k_, l_, val);
    return *this;
  }
  operator float() { return g_->transitions_.get(k_, l_).weight; }
  Transition* operator& () { return &hmm_->transitions_.mutating_get(k_, l_); }

 private:
  CRF<Alphabet>* g_;
  int k_;
  int l_;
};  // class TransitionAdapter

}  // namespace cs

#endif  // SRC_FACTOR_GRAPH_H_
