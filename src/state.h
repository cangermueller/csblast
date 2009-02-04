#ifndef CS_STATE_H
#define CS_STATE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The state object of a context HMM.

#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <vector>
#include <utility>

#include "context_profile.h"
#include "exception.h"
#include "context_profile.h"
#include "shared_ptr.h"
#include "transition.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class HMM;

template<class Alphabet_T>
class State : public ContextProfile<Alphabet_T>
{
  public:
    typedef std::list<Transition> transition_list;
    typedef transition_list::iterator transition_iterator;
    typedef transition_list::const_iterator const_transition_iterator;

    // Constructs HMM state from serialized state read from input stream.
    State(std::istream& in) : index_(0) {} //read(in) }
    // Constructs HMM state with given profile and all transitions initialized to zero.
    State(const ContextProfile<Alphabet_T>& profile);

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



// Returns the transition probability from this state to state k.
template<class Alphabet_T>
inline float State<Alphabet_T>::transition_probability_to(int k) const
{
    const_transition_iterator ti = find_if(out_transitions_.begin(), out_transitions_.end(), TransitsState(k));
    return ti == out_transitions_.end() ? 0.0f : ti->probability;
}

// Returns the transition probability from state k to this state.
template<class Alphabet_T>
inline float State<Alphabet_T>::transition_probability_from(int k) const
{
    const_transition_iterator ti = find_if(in_transitions_.begin(), in_transitions_.end(), TransitsState(k));
    return ti == out_transitions_.end() ? 0.0f : ti->probability;
}

template<class Alphabet_T>
void State<Alphabet_T>::read_header(std::istream& in)
{
    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("index") != std::string::npos)
        index_ = atoi(tmp.c_str() + 5);

    Profile::read_header(in);
}

template<class Alphabet_T>
void State<Alphabet_T>::read_body(std::istream& in)
{
    // TODO
}

template<class Alphabet_T>
void State<Alphabet_T>::write_header(std::ostream& out) const
{
    out << "index\t" << index_ << std::endl;
    Profile::write_header(out);
}

template<class Alphabet_T>
void State<Alphabet_T>::write_body(std::ostream& out) const
{
    // TODO
}

template<class Alphabet_T>
void State<Alphabet_T>::print(std::ostream& out) const
{
    out << "State " << index_ << ":" << std::endl;
    Profile<Alphabet_T>::print(out);
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
