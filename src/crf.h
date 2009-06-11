// Copyright 2009, Andreas Biegert

#ifndef SRC_CRF_H_
#define SRC_CRF_H_

#include <cstdlib>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "context_weight_state.h"
#include "chain_graph-inl.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs {

// A conditional random field (CRF) that stores context information in
// form of states of context weights and state transition weights.
template<class Alphabet>
class CRF : public ChainGraph<Alphabet, ContextWeightState> {
 public:
  // Public typedefs
  typedef ContextWeightState<Alphabet> State;
  typedef shared_ptr<State> StatePtr;
  typedef std::vector<StatePtr> StateVec;
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
  // initializer. State profiles are initially set to lin-space.
  CRF(int num_states,
      int num_cols,
      const StateInitializer<Alphabet, ContextWeightState>& st_init,
      const TransitionInitializer<Alphabet, ContextWeightState>& tr_init);

  virtual ~CRF() {}

  // Adds a profile as state to the CRF and returns its state index.
  virtual int AddState(const Profile<Alphabet>& profile);

 protected:
  // Needed to access names in templatized ChainGraph base class
  using ChainGraph<Alphabet, ContextWeightState>::full;
  using ChainGraph<Alphabet, ContextWeightState>::num_cols;
  using ChainGraph<Alphabet, ContextWeightState>::num_states;
  using ChainGraph<Alphabet, ContextWeightState>::states_;
  using ChainGraph<Alphabet, ContextWeightState>::Read;

 private:
  // Class identifier
  static const char* kClassID;

  // Return serialization class identity.
  virtual const char* class_id() const { return kClassID; }
};  // class CRF


// Strategy that randomly samples states from training profiles.
template<class Alphabet>
class SamplingStateInitializerCRF :
      public SamplingStateInitializer<Alphabet, ContextWeightState> {
 public:
  typedef shared_ptr< CountProfile<Alphabet> > ProfilePtr;
  typedef typename std::vector<ProfilePtr> ProfileVec;

  SamplingStateInitializerCRF(ProfileVec profiles,
                              float sample_rate,
                              const Pseudocounts<Alphabet>* pc = NULL,
                              float pc_admixture = 1.0f);
};

// Strategy that uses context profiles from a profile library to initialize
// states in the factor graph.
template<class Alphabet>
class LibraryBasedStateInitializerCRF :
      public LibraryBasedStateInitializer<Alphabet, ContextWeightState> {
 public:
  LibraryBasedStateInitializerCRF(const ProfileLibrary<Alphabet>* lib);
};


// Strategy that initializes transitions homogeneously.
template<class Alphabet>
class HomogeneousTransitionInitializerCRF
    : public HomogeneousTransitionInitializer<Alphabet, ContextWeightState> {
 public:
  HomogeneousTransitionInitializerCRF() {}
};

// Strategy that initializes transition probabilities at random.
template<class Alphabet>
class RandomTransitionInitializerCRF :
      public RandomTransitionInitializer<Alphabet, ContextWeightState> {
 public:
  RandomTransitionInitializerCRF() {}
};


// Erases all transitions whose log-value is equal to or below a given threshold.
template<class Alphabet>
void Sparsify(CRF<Alphabet>& crf, float threshold);

}  // namespace cs

#endif  // SRC_CRF_H_
