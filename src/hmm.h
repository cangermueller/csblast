// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_H_
#define SRC_HMM_H_

#include <cstdlib>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <vector>

#include "globals.h"
#include "co_emission.h"
#include "context_profile_state.h"
#include "factor_graph-inl.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs {

// A Hidden Markov Model that stores context information in form of states
// of context profiles and inter state transition probabilities.
template<class Alphabet>
class HMM : public FactorGraph<Alphabet, ContextProfileState> {
 public:
  // Public typedefs
  typedef shared_ptr< ContextProfileState<Alphabet> > StatePtr;
  typedef std::vector<StatePtr> StateVec;
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
      const StateInitializer<Alphabet, ContextProfileState>& st_init,
      const TransitionInitializer<Alphabet, ContextProfileState>& tr_init);

  virtual ~HMM() {}

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
  // Needed to access names in templatized FactorGraph base class
  using FactorGraph<Alphabet, ContextProfileState>::full;
  using FactorGraph<Alphabet, ContextProfileState>::num_cols;
  using FactorGraph<Alphabet, ContextProfileState>::num_states;
  using FactorGraph<Alphabet, ContextProfileState>::Read;
  using FactorGraph<Alphabet, ContextProfileState>::ReadHeader;
  using FactorGraph<Alphabet, ContextProfileState>::WriteHeader;
  using FactorGraph<Alphabet, ContextProfileState>::states_;
  using FactorGraph<Alphabet, ContextProfileState>::states_begin;
  using FactorGraph<Alphabet, ContextProfileState>::states_end;
  using FactorGraph<Alphabet, ContextProfileState>::kBufferSize;

  // Initializes the HMM from a serialized HMM read from stream.
  virtual void ReadHeader(FILE* fin);

 private:
  // Class identifier
  static const char* kClassID;

  // Return serialization class identity.
  virtual const char* class_id() const { return kClassID; }

  // Flag indicating if HMM profile probabilities are in log- or linspace
  bool states_logspace_;
};  // class HMM


// Strategy that randomly samples states from training profiles.
template<class Alphabet>
class SamplingStateInitializerHMM :
      public SamplingStateInitializer<Alphabet, ContextProfileState> {
 public:
  typedef shared_ptr< CountProfile<Alphabet> > ProfilePtr;
  typedef typename std::vector<ProfilePtr> ProfileVec;

  SamplingStateInitializerHMM(ProfileVec profiles,
                              float sample_rate,
                              const Pseudocounts<Alphabet>* pc = NULL,
                              float pc_admixture = 1.0f);
};

// Strategy that uses context profiles from a profile library to initialize
// states in the factor graph.
template<class Alphabet>
class LibraryStateInitializerHMM :
      public LibraryStateInitializer<Alphabet, ContextProfileState> {
 public:
  LibraryStateInitializerHMM(const ProfileLibrary<Alphabet>* lib);
};


// Strategy that initializes transitions homogeneously.
template<class Alphabet>
class HomogeneousTransitionInitializerHMM
    : public HomogeneousTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  HomogeneousTransitionInitializerHMM() {}
};

// Strategy that initializes transition probabilities at random.
template<class Alphabet>
class RandomTransitionInitializerHMM :
      public RandomTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  RandomTransitionInitializerHMM() {}
};

// Strategy that uses the co-emission scores of profile pairs to determine
// transition probabilities. No transition is added if the co-emission is below
// given threshold.
template<class Alphabet>
class CoEmissionTransitionInitializerHMM
    : public CoEmissionTransitionInitializer<Alphabet, ContextProfileState> {
 public:
  CoEmissionTransitionInitializerHMM(const SubstitutionMatrix<Alphabet>* sm,
                                     float score_thresh);
};

}  // namespace cs

#endif  // SRC_HMM_H_
