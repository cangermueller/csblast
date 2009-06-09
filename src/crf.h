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
#include "factor_graph-inl.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs {

// A conditional random field (CRF) that stores context information in
// form of states of context weights and state transition weights.
template<class Alphabet>
class CRF : public FactorGraph<Alphabet, ContextWeightState> {
 public:
  // Public typedefs
  typedef shared_ptr< ContextWeightState<Alphabet> > StatePtr;
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
  // Needed to access names in templatized FactorGraph base class
  using FactorGraph<Alphabet, ContextWeightState>::full;
  using FactorGraph<Alphabet, ContextWeightState>::num_cols;
  using FactorGraph<Alphabet, ContextWeightState>::num_states;
  using FactorGraph<Alphabet, ContextWeightState>::Read;
  using FactorGraph<Alphabet, ContextWeightState>::ReadHeader;
  using FactorGraph<Alphabet, ContextWeightState>::WriteHeader;
  using FactorGraph<Alphabet, ContextWeightState>::states_;
  using FactorGraph<Alphabet, ContextWeightState>::states_begin;
  using FactorGraph<Alphabet, ContextWeightState>::states_end;
  using FactorGraph<Alphabet, ContextWeightState>::kBufferSize;

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
class LibraryStateInitializerCRF :
      public LibraryStateInitializer<Alphabet, ContextWeightState> {
 public:
  LibraryStateInitializerCRF(const ProfileLibrary<Alphabet>* lib);
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

}  // namespace cs

#endif  // SRC_CRF_H_
