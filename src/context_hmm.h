#ifndef CS_COUNTS_PROFILE_H
#define CS_COUNTS_PROFILE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A hidden Markov model that stores context information in the form of
// state-specific context profiles and state transition probabilities.

#include <list>
#include <vector>
#include <utility>

#include "context_profile.h"
#include "shared_ptr.h"

namespace cs
{

class ContextHMM
{
  public:
    struct Transition
    {
        // Transition probability.
        float probability;
        // Index of state from/to which the transition connects.
        int state;
    };

    class State
    {
      public:
        // Constructs HMM state from serialized state read from input stream.
        State(std::istream& in, const ContextHMM* hmm);

        // Index of state in HMM states vector.
        int index;
        // Context profile with residue emission probabilities.
        ContextProfile profile;
        // List of in-transitions.
        std::list<Transition> in_transitions;
        // List of out-transitions.
        std::list<Transition> out_transitions;
    };

    // Constructs context HMM from serialized HMM read from input stream.
    ContextHMM(std::istream& in, const SequenceAlphabet* alphabet);

    virtual ~ContextHMM() {}

    // Returns the number of states in the HMM (including BEGIN- and END-state)
    float size(int i) const { return states_.size(); }
    // Access methods for HMM state i.
    State& operator[](int i) { return *states_[i]; }
    const State& operator[](int i) const { return *states_[i]; }
    // Accessors for BEGIN state.
    State& begin_state() { return *states_.front(); }
    const State& begin_state() const { return *states_.front(); }
    // Accessors for END state.
    State& end_state() { return *states_.back(); }
    const State& end_state() const { return *states_.back(); }

  private:
    // The HMM states ordered by index (BEGIN=0, 1, 2, ..., K, END=K+1)
    std::vector< shared_ptr<State> > states_;
};  // ContextHMM

}  // cs

#endif
