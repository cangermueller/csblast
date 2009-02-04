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
#include "sparse_matrix.h"
#include "state.h"

namespace cs
{

struct Transition
{
    // Simple Constructors
    Transition() : from(0), to(0), probability(0.0f) {}
    Transition(int f, int t, float p) : from(f), to(t), probability(p) {}
    ~Transition() {}

    // Index of state from which the transition connects.
    int from;
    // Index of state to which the transition connects.
    int to;
    // Transition probability.
    float probability;
};  // Transition

template<class Alphabet_T>
class HMM
{
  public:
    // Public typedefs
    typedef std::vector< shared_ptr<State> >::iterator state_iterator;
    typedef std::vector< shared_ptr<State> >::const_iterator const_state_iterator;
    typedef sparse_matrix<Transition>::iterator transition_iterator;
    typedef sparse_matrix<Transition>::const_iterator const_transition_iterator;

    // BEGIN and END states have both index 0, this is because is all transitions to state 0
    // go to the END state and all transitions from state 0 start from the BEGIN. There are no
    // in-transitions for the BEGIN state and no out-transitions for the END state.
    static const int BEGIN_END_STATE = 0;

    // Cosntructs an empty HMM of given size without any states or transitions.
    HMM(int size) : size_(size), transitions_(size, size) {}
    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in) : { read(in); }
    // Constructs context HMM with the help of an initializer.
    //HMM(shared_ptr<Initializer> initializer) : alphabet_(initializer.alphabet()) { initializer.init(*this); }

    virtual ~HMM() {}

    // Returns the number of states in the HMM (not counting the BEGIN/END state)
    int size() const { return states_.size(); }
    // Returns the number of states in the HMM (not counting the BEGIN/END state)
    int num_states() const { return states_.size(); }
    // Access method for state i, where i is from interval [1,K] (that is no access to BEGIN/END state).
    const State<Alphabet_T>& operator[](int i) const { return *states_[i-1]; }
    // Returns the transition probability from state k to state l in runtime O(#transitions).
    float transition_probability(int k, int l) const { transitions_.get(k,l).probability; }
    // Sets the transition probability from state k to state l in runtime O(#transitions).
    void set_transition_probability(int k, int l, float prob) { transitions_.set(k, l, Transition(k, l, prob)); }
    // Removes all transitions with probability below or equal to given threshold.
    void sparsify(float threshold);
    // Clears all states and transitions.
    void clear() { states_.clear(); transitions_.clear(); }
    // Clears all transitions but leaves profile of states untouched.
    void clear_transitions();
    // Adds the given profile as state to the HMM and returns its state index. Note that number of columns must be odd!
    int add_profile(const Profile<Alphabet_T>& profile);
    // Returns an iterator to a list of pointers of states.
    state_iterator states_begin() { return states_.begin(); }
    // Returns an iterator pointing past the end of a list of pointers of states.
    state_iterator states_end() { return states_.end(); }
     // Returns a const iterator to a list of pointers of states.
    const_state_iterator states_begin() const { return states_.begin(); }
    // Returns a const iterator pointing past the end of a list of pointers of states.
    const_state_iterator states_end() const { return states_.end(); }
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
    // Initializes the HMM with the mandatory BEGIN/END state.
    void init();

    // Final number of states in the HMM (excl. BEGIN/END state)
    const int num_states_;
    // HMM states ordered by index (1, 2, ..., size)
    std::vector< shared_ptr<State> > states_;
    // Sparse matrix with state transitions
    sparse_matrix<Transition> transitions_;
    // // Pointer to .
    // shared_ptr<Initializer> initializer_;
};  // HMM







const char HMM<Alphabet_T>::CLASS_IDENTITY[] = "HMM";

HMM(int size) : size_(size), transitions_(size, size) {}

// Returns the transition probability from state k to state l in runtime O(#transitions).
template<class Alphabet_T>
inline float HMM<Alphabet_T>::transition_probability(int k, int l) const
{
    if (k > BEGIN_END_STATE)
        return states_[k-1]->transition_probability_to(l);
    else
        return (*states_[l-1]).transition_probability_from(k);
}

// Sets the transition probability from state k to state l in runtime O(#transitions).
template<class Alphabet_T>
inline void HMM<Alphabet_T>::set_transition_probability(int k, int l, float prob)
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

template<class Alphabet_T>
void HMM<Alphabet_T>::read(std::istream& in)
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

template<class Alphabet_T>
void HMM<Alphabet_T>::write(std::ostream& out) const
{
    out << CLASS_IDENTITY << std::endl;
    out << "size\t" << size() << std::endl;
    for (std::pair<const_state_iterator, const_state_iterator> sp = states(); sp.first != sp.second; ++sp.first)
        (**sp.first).write(out);
}

template<class Alphabet_T>
void HMM<Alphabet_T>::sparsify(float threshold)
{
    std::pair<state_iterator, state_iterator> sp = states();
    for (std::pair<state_iterator, state_iterator> sp = states(); sp.first != sp.second; ++sp.first) {
        (*sp.first)->in_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
        (*sp.first)->out_transitions_.remove_if( FailsProbabilityThreshold(threshold) );
    }
}

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
        return transition.from == state_;
    }

  private:
    const int state_;
};

}  // cs

#endif
