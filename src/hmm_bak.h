// Copyright 2009, Andreas Biegert

#ifndef CS_HMM_H_
#define CS_HMM_H_

#include <stdlib.h>
#include <stdio.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "co_emission.h"
#include "context_profile_state.h"
#include "chain_graph-inl.h"
#include "profile.h"
#include "profile_library.h"
#include "pseudocounts-inl.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"

namespace cs {

// Forwards declarations
template< class Alphabet, template<class> class State >
class StateInitializer;

template< class Alphabet, template<class> class State >
class TransitionInitializer;

template< class Alphabet, template<class> class State >
class SamplingStateInitializer;

template< class Alphabet, template<class> class State >
class LibraryBasedStateInitializer;

template< class Alphabet, template<class> class State >
class HomogeneousTransitionInitializer;

template< class Alphabet, template<class> class State >
class RandomTransitionInitializer;

template< class Alphabet, template<class> class State >
class CoEmissionTransitionInitializer;


// A Hidden Markov Model that stores context information in form of states
// of context profiles and inter state transition probabilities.
template<class Alphabet>
class Hmm : public ChainGraph<Alphabet, ContextProfileState> {
 public:
  // Public typedefs
  typedef ContextProfileState<Alphabet> State;
  typedef shared_ptr<State> StatePtr;
  typedef std::vector<StatePtr> StateVec;
  typedef sparse_matrix<Transition> TransitionMatrix;
  typedef typename StateVec::iterator StateIter;
  typedef typename StateVec::const_iterator ConstStateIter;
  typedef typename TransitionMatrix::nonempty_iterator TransitionIter;
  typedef typename TransitionMatrix::const_nonempty_iterator ConstTransitionIter;

  // Needed to access names in templatized ChainGraph base class
  using ChainGraph<Alphabet, ContextProfileState>::full;
  using ChainGraph<Alphabet, ContextProfileState>::num_cols;
  using ChainGraph<Alphabet, ContextProfileState>::num_states;

  // Constructs an empty HMM of given size without any states or transitions.
  Hmm(int num_states, int num_cols);
  // Constructs context HMM from serialized HMM read from input stream.
  explicit Hmm(FILE* fin);
  // Constructs context HMM with the help of a state- and a transition-
  // initializer. State profiles are initially set to lin-space.
  Hmm(int num_states,
      int num_cols,
      const ::cs::StateInitializer<Alphabet, ContextProfileState>& st_init,
      const ::cs::TransitionInitializer<Alphabet, ContextProfileState>& tr_init);

  virtual ~Hmm() {}

  // Adds a profile as state to the HMM and returns its state index.
  virtual int AddState(const Profile<Alphabet>& profile);
  // Returns true if state profiles are in logspace.
  bool states_logspace() const { return states_logspace_; }
  // Transforms state profiles to logspace.
  void TransformStatesToLogSpace();
  // Transforms state profiles to linspace.
  void TransformStatesToLinSpace();
  // Writes header section for serialization.
  virtual void WriteHeader(FILE* fout) const;

 protected:
  // Needed to access names in templatized ChainGraph base class
  using ChainGraph<Alphabet, ContextProfileState>::Read;
  using ChainGraph<Alphabet, ContextProfileState>::ReadHeader;
  using ChainGraph<Alphabet, ContextProfileState>::WriteHeader;
  using ChainGraph<Alphabet, ContextProfileState>::states_;
  using ChainGraph<Alphabet, ContextProfileState>::states_begin;
  using ChainGraph<Alphabet, ContextProfileState>::states_end;
  using ChainGraph<Alphabet, ContextProfileState>::kBufferSize;

  // Initializes the HMM from a serialized HMM read from stream.
  virtual void ReadHeader(FILE* fin);

 private:
  // Class identifier
  static const char* kClassID;

  // Return serialization class identity.
  virtual const char* class_id() const { return kClassID; }

  // Flag indicating if HMM profile probabilities are in log- or linspace
  bool states_logspace_;
};  // class Hmm


// Strategy that randomly samples states from training profiles.
template<class Alphabet>
class SamplingStateInitializerHmm
    : public SamplingStateInitializer<Alphabet, ContextProfileState> {
 public:
  typedef shared_ptr< CountProfile<Alphabet> > ProfilePtr;
  typedef typename std::vector<ProfilePtr> ProfileVec;

  SamplingStateInitializerHmm(ProfileVec profiles,
                              float sample_rate,
                              const Pseudocounts<Alphabet>* pc = NULL,
                              float pc_admixture = 1.0f);
};

// Strategy that uses context profiles from a profile library to initialize
// states in the factor graph.
template<class Alphabet>
class LibraryBasedStateInitializerHmm
    : public LibraryBasedStateInitializer<Alphabet, ContextProfileState> {
 public:
  LibraryBasedStateInitializerHmm(const ProfileLibrary<Alphabet>* lib);
};


// Strategy that initializes transitions homogeneously.
template<class Alphabet>
class HomogeneousTransitionInitializerHmm
    : public HomogeneousTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  HomogeneousTransitionInitializerHmm() {}
};

// Strategy that initializes transition probabilities at random.
template<class Alphabet>
class RandomTransitionInitializerHmm :
      public RandomTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  RandomTransitionInitializerHmm() {}
};

// Strategy that uses the co-emission scores of profile pairs to determine
// transition probabilities. No transition is added if the co-emission is below
// given threshold.
template<class Alphabet>
class CoEmissionTransitionInitializerHmm
    : public CoEmissionTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  CoEmissionTransitionInitializerHmm(const SubstitutionMatrix<Alphabet>* sm,
                                     float score_thresh);
};

}  // namespace cs

#endif  // CS_HMM_H_
