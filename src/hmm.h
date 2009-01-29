#ifndef CS_HMM_H
#define CS_HMM_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A hidden Markov model that stores context information in the form of
// state-specific context profiles and state transition probabilities.

#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <vector>
#include <utility>

//#include "initializer.h"
#include "profile.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"

namespace cs
{

class HMM
{
  public:
    // Forward declarations
    struct Transition;
    class State;

    // Public typedefs
    typedef std::list<Transition> transition_list;
    typedef transition_list::iterator transition_iterator;
    typedef transition_list::const_iterator const_transition_iterator;
    typedef std::vector< shared_ptr<State> > state_vector;
    typedef state_vector::iterator state_iterator;
    typedef state_vector::const_iterator const_state_iterator;

    struct Transition
    {
        // Index of state from/to which the transition connects.
        int state;
        // Transition probability.
        float probability;
    };

    class State : public Profile
    {
      public:
        // Constructs HMM state from serialized state read from input stream.
        State(std::istream& in, const SequenceAlphabet* alphabet) : Profile(alphabet), index_(0) {}
        virtual ~State() {}

        // Returns the index of this state.
        int index() const { return index_; }
        // Returns number of in-transitions.
        int num_in_transitions() const { return in_transitions_.size(); }
        // Returns number of out-transitions.
        int num_out_transitions() const { return out_transitions_.size(); }
        // Returns the transition probability from this state to state k.
        float transition_probability_to(int k) const;
        // Returns the transition probability from state k to this state.
        float transition_probability_from(int k) const;
        // Returns pair of in-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<transition_iterator, transition_iterator> in_transitions()
        { return make_pair(in_transitions_.begin(), in_transitions_.end()); }
        // Returns pair of out-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<transition_iterator, transition_iterator> out_transitions()
        { return make_pair(out_transitions_.begin(), out_transitions_.end()); }
        // Returns pair of const in-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<const_transition_iterator, const_transition_iterator> in_transitions() const
        { return make_pair(in_transitions_.begin(), in_transitions_.end()); }
        // Returns pair of const out-transition iterators (the first iterator points to the "beginning" the second "past the end")
        std::pair<const_transition_iterator, const_transition_iterator> out_transitions() const
        { return make_pair(out_transitions_.begin(), out_transitions_.end()); }

      protected:
        // // Reads and initializes serialized scalar data members from stream.
        // virtual void read_header(std::istream& in);
        // // Reads and initializes array data members from stream.
        // virtual void read_body(std::istream& in);
        // // Writes serialized scalar data members to stream.
        // virtual void write_header(std::ostream& out) const;
        // // Writes serialized array data members to stream.
        // virtual void write_body(std::ostream& out) const;
        // // Prints the profile in human-readable format to output stream.
        // virtual void print(std::ostream& out) const;

      private:
        // Returns serialization class identity.
        virtual const std::string class_identity() const { static std::string id("State"); return id;}

        // Index of state in HMM states vector.
        int index_;
        // List of in-transitions.
        transition_list in_transitions_;
        // List of out-transitions.
        transition_list out_transitions_;
    };  // State


    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in, const SequenceAlphabet* alphabet) : alphabet_(alphabet) { read(in); }
    virtual ~HMM() {}

    // Returns the number of states in the HMM (not counting the BEGIN/END state)
    int size() const { return states_.size(); }
    // Access methods for state i, where i is from interval [1,K] (that is no access to BEGIN/END state).
    State& operator[](int i) { return *states_[i-1]; }
    const State& operator[](int i) const { return *states_[i-1]; }
    // Returns the transition probability from state k to state l (runtime O(1) for BEGIN/END state; O(K) for others).
    float transition_probability(int k, int l) const;
    // Returns pair of state iterators (the first iterator points to the "beginning" the second "past the end")
    std::pair<state_iterator, state_iterator> states() { return make_pair(states_.begin(), states_.end()); }
    // Returns pair of const state iterators (the first iterator points to the "beginning" the second "past the end")
    std::pair<const_state_iterator, const_state_iterator> states() const { return make_pair(states_.begin(), states_.end()); }
    // Initializes the HMM from a serialized HMM read from stream.
    void read(std::istream& in);
    // Writes the HMM in serialization format to output stream.
    void write(std::ostream& out) const;
    // Returns the sequence alphabet.
    const SequenceAlphabet* alphabet() const { return alphabet_; }

    friend std::ostream& operator<< (std::ostream& out, const HMM& hmm);

  private:
    // Class identity keyword for serialization.
    static const char kClassIdentity[];

    // Prints the HMM in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

    // HMM states ordered by index (1, 2, ..., K)
    state_vector states_;
    // // Pointer to .
    // shared_ptr<Initializer> initializer_;
    // The underlying alphabet of the sequence.
    const SequenceAlphabet* alphabet_;
};  // HMM



// Returns the transition probability from state k to state l (runtime O(1) for BEGIN/END state; O(K) for others).
inline float HMM::transition_probability(int k, int l) const
{
    return (*states_[k-1]).transition_probability_to(l);
}

// Returns the transition probability from this state to state k.
inline float HMM::State::transition_probability_to(int k) const
{
    for (const_transition_iterator ti = out_transitions_.begin(); ti != out_transitions_.end(); ++ti)
        if (ti->state == k) return ti->probability;
    return 0.0f;
}

// Returns the transition probability from state k to this state.
inline float HMM::State::transition_probability_from(int k) const
{
    for (const_transition_iterator ti = in_transitions_.begin(); ti != in_transitions_.end(); ++ti)
        if (ti->state == k) return ti->probability;
    return 0.0f;
}

// Prints HMM in human-readable format for debugging.
inline std::ostream& operator<< (std::ostream& out, const HMM& hmm)
{
    hmm.print(out);
    return out;
}

}  // cs

#endif
