#ifndef CS_CONTEXT_HMM_H
#define CS_CONTEXT_HMM_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A hidden Markov model that stores context information in the form of
// state-specific context profiles and state transition probabilities.

#include <algorithm>
#include <functional>
#include <list>
#include <vector>
#include <utility>

#include "context_profile.h"
#include "shared_ptr.h"

namespace cs
{

// Forward declaraions
class ContextHMM;
class ContextHMM::Transition;
class ContextHMM::State;

class ContextHMM
{
  public:
    // public typedefs
    typedef std::list<Transition> transition_list;
    typedef transition_list::iterator transition_iterator;
    typedef transition_list::const_iterator const_transition_iterator;

    struct Transition
    {
        // Index of state from/to which the transition connects.
        int state;
        // Transition probability.
        float probability;
    }; // Transition

    class State : public Profile
    {
      public:
        // Constructs HMM state from serialized state read from input stream.
        State(std::istream& in, shared_ptr<const ContextHMM> hmm);

        // Returns number of in-transitions.
        int num_in_transitions() const { return in_transitions.size(); }
        // Returns number of out-transitions.
        int num_out_transitions() const { return out_transitions.size(); }
        // Returns the transition probability from this state to state k.
        float transition_probability_to(int k);
        // Returns the transition probability from state k to this state.
        float transition_probability_from(int k);

        // Returns pair of in-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<transition_iterator, transition_iterator> in_transitions()
        { return makepair(in_transitions_.begin(), in_transitions_.end()); }
        // Returns pair of out-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<transition_iterator, transition_iterator> out_transitions()
        { return makepair(out_transitions_.begin(), out_transitions_.end()); }
        // Returns pair of const in-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<const_transition_iterator, const_transition_iterator> in_transitions() const
        { return makepair(in_transitions_.begin(), in_transitions_.end()); }
        // Returns pair of const out-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<const_transition_iterator, const_transition_iterator> out_transitions() const
        { return makepair(out_transitions_.begin(), out_transitions_.end()); }

      private:
        // Returns serialization class identity.
        virtual const std::string class_identity() const { static std::string id("State"); return id;}

        // Index of state in HMM states vector.
        int index_;
        // List of in-transitions.
        transition_list in_transitions_;
        // List of out-transitions.
        transition_list out_transitions_;
        // Pointer to parent HMM.
        shared_ptr<const ContextHMM> hmm_;
    };  // State


    // Constructs context HMM from serialized HMM read from input stream.
    ContextHMM(std::istream& in, const SequenceAlphabet* alphabet);

    virtual ~ContextHMM() {}

    // Returns the number of states in the HMM (not counting the BEGIN/END state)
    int size() const { return states_.size(); }
    // Access methods for state i, where i is from interval [1,K] (that is no access to BEGIN/END state).
    State& operator[](int i) { return *states_[i-1]; }
    const State& operator[](int i) const { return *states_[i-1]; }
    // Returns the transition probability from state k to state l (runtime O(1) for BEGIN/END state; O(K) for others).
    float transition_probability(int k, int l) const;
    // Returns the sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

  private:
    // HMM states ordered by index (1, 2, ..., K)
    std::vector< shared_ptr<State> > states_;
    // The underlying alphabet of the sequence.
    const SequenceAlphabet* alphabet_;
};  // ContextHMM



// Binary function for finding transitions in transition list.
struct transition_state_equal_to : public std::binary_function<ContextHMM::Transition, ContextHMM::Transition, bool>
{
    bool operator() (const ContextHMM::Transition& a, const ContextHMM::Transition& b) const
    { return a.state == b.state; }
};

inline float ContextHMM::transition_probability(int k, int l) const
{
    return states_[k-1]->transition_probability_to(l);
}

inline float ContextHMM::State::transition_probability_to(int k) const
{
    const_transition_iterator ti = find_if(out_transitions.begin(), out_transitions.end(), bind2nd(transition_state_equal_to(), k));
    return ti == out_transitions.end() ? 0.0f : ti->state;
}

inline float ContextHMM::State::transition_probability_from(int k) const
{
    const_transition_iterator ti = find_if(in_transitions_.begin(), in_transitions_.end(), bind2nd(transition_state_equal_to(), k));
    return ti == in_transitions_.end() ? 0.0f : ti->state;
}

}  // cs

#endif
