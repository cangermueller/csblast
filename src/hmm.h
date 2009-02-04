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
#include "context_profile.h"
#include "exception.h"
#include "context_profile.h"
#include "shared_ptr.h"

namespace cs
{

// Forward declarations
template<class AlphabetType>
class Transition;
template<class AlphabetType>
class State;

template<class AlphabetType>
class HMM
{
  public:
    // Public typedefs
    typedef std::vector< shared_ptr<State> > state_vector;
    typedef state_vector::iterator state_iterator;
    typedef state_vector::const_iterator const_state_iterator;

    // BEGIN and END states have both index 0, this is because is all transitions to state 0
    // go to the END state and all transitions from state 0 start from the BEGIN. There are no
    // in-transitions for the BEGIN state and no out-transitions for the END state.
    static const int BEGIN_END_STATE = 0;

    // Cosntructs an empty HMM without any states or transitions.
    HMM() {}
    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in) : { read(in); }
    // Constructs context HMM with the help of an initializer.
    //HMM(shared_ptr<Initializer> initializer) : alphabet_(initializer.alphabet()) { initializer.init(*this); }

    virtual ~HMM() {}

    // Returns the number of states in the HMM (not counting the BEGIN/END state)
    int size() const { return states_.size(); }
    // Access method for state i, where i is from interval [1,K] (that is no access to BEGIN/END state).
    const State<AlphabetType>& operator[](int i) const { return *states_[i-1]; }
    // Returns the transition probability from state k to state l in runtime O(#transitions).
    float transition_probability(int k, int l) const;
    // Sets the transition probability from state k to state l in runtime O(#transitions).
    void set_transition_probability(int k, int l, float prob);
    // Removes all transitions with probability below or equal to given threshold.
    void sparsify(float threshold);
    // Clears all states and transitions.
    void clear() { states_.clear(); }
    // Clears all transitions but leaves profile of states untouched.
    void clear_transitions();
    // Adds the given state to the HMM and returns the index of the newly added state.
    int add_state(const State<AlphabetType>& state);
    // Adds the given profile as state to the HMM and returns the index of the newly added state.
    int add_profile(const ContextProfile<AlphabetType>& profile);
    // Returns pair of state iterators (the first iterator points to the "beginning" the second "past the end")
    std::pair<state_iterator, state_iterator> states() { return make_pair(states_.begin(), states_.end()); }
    // Returns pair of const state iterators (the first iterator points to the "beginning" the second "past the end")
    std::pair<const_state_iterator, const_state_iterator> states() const { return make_pair(states_.begin(), states_.end()); }
    // Initializes the HMM from a serialized HMM read from stream.
    void read(std::istream& in);
    // Writes the HMM in serialization format to output stream.
    void write(std::ostream& out) const;

    // Prints HMM in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const HMM& hmm)
    {
        hmm.print(out);
        return out;
    }

  private:
    // Class identity keyword for serialization.
    static const char CLASS_IDENTITY[];

    // Prints the HMM in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

    // HMM states ordered by index (1, 2, ..., K)
    state_vector states_;
    // // Pointer to .
    // shared_ptr<Initializer> initializer_;
};  // HMM



template<class AlphabetType>
class State : public ContextProfile<AlphabetType>
{
  public:
    typedef std::list<Transition> transition_list;
    typedef transition_list::iterator transition_iterator;
    typedef transition_list::const_iterator const_transition_iterator;

    // Constructs HMM state from serialized state read from input stream.
    State(std::istream& in) : index_(0) {} //read(in) }
    // Constructs HMM state with given profile and all transitions initialized to zero.
    State(const ContextProfile<AlphabetType>& profile);

    virtual ~State() {}

    // Returns the index of this state.
    int index() const { return index_; }
    // Returns number of in-transitions.
    int num_in_transitions() const { return in_transitions_.size(); }
    // Returns number of out-transitions.
    int num_out_transitions() const { return out_transitions_.size(); }
    // Returns the transition probability from this state to state k in runtime O(#transitions).
    float transition_probability_to(int k) const;
    // Returns the transition probability from state k to this state in runtime O(#transitions).
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

    friend class HMM;

  protected:
    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes array data members from stream.
    virtual void read_body(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized array data members to stream.
    virtual void write_body(std::ostream& out) const;
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

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



class Transition
{
  public:
    // Simple Constructors
    Transition() : state(0), probability(0.0f) {}
    Transition(int s, float p) : state(s), probability(p) {}
    ~Transition() {}

    // Index of state from/to which the transition connects.
    int state;
    // Transition probability.
    float probability;
};  // Transition




// Functor predicate that returns true if a transition has a transition probability below given threshold.
class FailsProbabilityThreshold : public std::unary_function<Transition, bool>
{
  public:
    FailsProbabilityThreshold(float threshold) : threshold_(threshold) {}

    bool operator() (const Transition& transition) const
    {
        return transition.probability < threshold_;
    }

  private:
    const float threshold_;
};

// Functor predicate that returns true if a transition transits from or to a given state.
class TransitsState : public std::unary_function<Transition, bool>
{
  public:
    TransitsState(int state) : state_(state) {}

    bool operator() (const Transition& transition) const
    {
        return transition.state == state_;
    }

  private:
    const int state_;
};




const char HMM<AlphabetType>::CLASS_IDENTITY[] = "HMM";

// Returns the transition probability from state k to state l in runtime O(#transitions).
template<class AlphabetType>
inline float HMM<AlphabetType>::transition_probability(int k, int l) const
{
    if (k > BEGIN_END_STATE)
        return states_[k-1]->transition_probability_to(l);
    else
        return (*states_[l-1]).transition_probability_from(k);
}

// Sets the transition probability from state k to state l in runtime O(#transitions).
template<class AlphabetType>
inline void HMM<AlphabetType>::set_transition_probability(int k, int l, float prob)
{
    if (prob == 0.0f) {
        // Remove transitions from transition lists if probability is zero
        if (k > BEGIN_END_STATE) (*states_[k-1]).out_transitions_.remove_if( TransitsState(l) );
        if (l > BEGIN_END_STATE) (*states_[l-1]).in_transitions_.remove_if( TransitsState(k) );

    } else {
        // Find correct in- and out-transitions and set their probability to new value.
        if (k > BEGIN_END_STATE) {
            // set out-transition of state k
            std::pair<transition_iterator, transition_iterator> tp = (*states_[k-1]).out_transitions();
            transition_iterator ti = find_if(tp.first, tp.second, TransitsState(l));
            if (ti == tp.second)  // transition does not yet exist
                (*states_[k-1]).out_transitions_.push_back(Transition(l, prob));
            else  // modify existing transition
                ti->probability = prob;
        }

        if (l > BEGIN_END_STATE) {
            // set out-transition of state k
            std::pair<transition_iterator, transition_iterator> tp = (*states_[l-1]).in_transitions();
            transition_iterator ti = find_if(tp.first, tp.second, TransitsState(k));
            if (ti == tp.second)  // transition does not yet exist
                (*states_[l-1]).in_transitions_.push_back(Transition(k, prob));
            else  // modify existing transition
                ti->probability = prob;
        }
    }
}

template<class AlphabetType>
void HMM<AlphabetType>::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading HMM from stream ...";

    // Check if stream actually contains a serialized HMM
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find(CLASS_IDENTITY) == std::string::npos)
        throw Exception("Bad format: serialized HMM does not start with '%s' keyword!", CLASS_IDENTITY);

    // Read number of states
    int size = 0;
    if (getline(in, tmp) && tmp.find("size") != std::string::npos)
        size = atoi(tmp.c_str() + 4);
    else
        throw Exception("Bad format: serialized profile does not contain 'size' record!");

    states_.clear();
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr<State> state(new State(in, alphabet_));
        states_.push_back(state);
    }

    LOG(DEBUG1) << *this;
}

template<class AlphabetType>
void HMM<AlphabetType>::write(std::ostream& out) const
{
    out << CLASS_IDENTITY << std::endl;
    out << "size\t" << size() << std::endl;
    for (std::pair<const_state_iterator, const_state_iterator> sp = states(); sp.first != sp.second; ++sp.first)
        (**sp.first).write(out);
}

template<class AlphabetType>
void HMM<AlphabetType>::sparsify(float threshold)
{
    std::pair<state_iterator, state_iterator> sp = states();
    for (std::pair<state_iterator, state_iterator> sp = states(); sp.first != sp.second; ++sp.first) {
        (*sp.first)->in_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
        (*sp.first)->out_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
    }
}



// Returns the transition probability from this state to state k.
template<class AlphabetType>
inline float State<AlphabetType>::transition_probability_to(int k) const
{
    const_transition_iterator ti = find_if(out_transitions_.begin(), out_transitions_.end(), TransitsState(k));
    return ti == out_transitions_.end() ? 0.0f : ti->probability;
}

// Returns the transition probability from state k to this state.
template<class AlphabetType>
inline float State<AlphabetType>::transition_probability_from(int k) const
{
    const_transition_iterator ti = find_if(in_transitions_.begin(), in_transitions_.end(), TransitsState(k));
    return ti == out_transitions_.end() ? 0.0f : ti->probability;
}

template<class AlphabetType>
void State<AlphabetType>::read_header(std::istream& in)
{
    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("index") != std::string::npos)
        index_ = atoi(tmp.c_str() + 5);

    Profile::read_header(in);
}

template<class AlphabetType>
void State<AlphabetType>::read_body(std::istream& in)
{
    // TODO
}

template<class AlphabetType>
void State<AlphabetType>::write_header(std::ostream& out) const
{
    out << "index\t" << index_ << std::endl;
    Profile::write_header(out);
}

template<class AlphabetType>
void State<AlphabetType>::write_body(std::ostream& out) const
{
    // TODO
}

template<class AlphabetType>
void State<AlphabetType>::print(std::ostream& out) const
{
    out << "State " << index_ << ":" << std::endl;
    Profile<AlphabetType>::print(out);
}



// Comparison functor that compares the index of two states.
struct StateIndexCompare : public std::binary_function<State, State, bool>
{
    bool operator()(const State& lhs, const State& rhs) const
    {
        return lhs.index() < rhs.index();
    }
};

// Comparison functor that compares the number of in-transitions of twho states.
struct NumInTransitionsCompare : public std::binary_function<State, State, bool>
{
    bool operator()(const State& lhs, const State& rhs) const
    {
        return lhs.num_in_transitions() < rhs.num_in_transitions();
    }
};

// Comparison functor that compares the number of out-transitions of twho states.
struct NumOutTransitionsCompare : public std::binary_function<State, State, bool>
{
    bool operator()(const State& lhs, const State& rhs) const
    {
        return lhs.num_out_transitions() < rhs.num_out_transitions();
    }
};

}  // cs

#endif
