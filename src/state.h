#ifndef CS_STATE_H
#define CS_STATE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The state object of a context HMM.

#include <google/sparsetable>

#include "context_profile.h"
#include "hmm.h"
#include "profile.h"

using google::sparsetable;

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class HMM;

template<class Alphabet_T>
class State : public ContextProfile<Alphabet_T>
{
  public:
    typedef sparsetable<Transition*>::nonempty_iterator transition_iterator;
    typedef sparsetable<Transition*>::const_nonempty_iterator const_transition_iterator;

    // Constructs dummy state for BEGIN/END state in HMM.
    State();
    // Constructs HMM state from serialized state read from input stream.
    explicit State(std::istream& in);
    // Constructs HMM state with given profile and all transitions initialized to zero.
    explicit State(int index, const Profile<Alphabet_T>& profile);

    virtual ~State() {}

    // Returns the index of this state.
    int index() const { return index_; }
    // Returns number of in-transitions.
    int num_in_transitions() const { return in_transitions_.num_nonempty(); }
    // Returns number of out-transitions.
    int num_out_transitions() const { return out_transitions_.num_nonempty(); }
    // Returns the transition probability from this state to state k in runtime O(#transitions).
    float transition_probability_to(int k) const { return out_transitions_.get(k)->probability; }
    // Returns the transition probability from state k to this state in runtime O(#transitions).
    float transition_probability_from(int k) const { return in_transitions_.get(k)->probability; }
    // Returns an iterator to start of list with non-null in-transition pointers.
    transition_iterator in_transitions_begin() { return in_transitions_.nonempty_begin(); }
    // Returns an iterator past the end of list with non-null in-transition pointers.
    transition_iterator in_transitions_end() { return in_transitions_.nonempty_end(); }
    // Returns an iterator to start of list with non-null out-transition pointers.
    transition_iterator out_transitions_begin() { return out_transitions_.nonempty_begin(); }
    // Returns an iterator past the end of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_end() { return out_transitions_.nonempty_end(); }
    // Returns a const iterator to start of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_begin() const { return in_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_end() const { return in_transitions_.nonempty_end(); }
    // Returns a const iterator to start of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_begin() const { return out_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_end() const { return out_transitions_.nonempty_end(); }

    // HMM needs access to transition tables.
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
    sparsetable<Transition*> in_transitions_;
    // List of out-transitions.
    sparsetable<Transition*> out_transitions_;
};  // State



template<class Alphabet_T>
State<Alphabet_T>::State()
        : index_(index),
          in_transition_(0),
          out_transition_(0)
{}

template<class Alphabet_T>
State<Alphabet_T>::State(std::istream& in)
        : index_(0),
          in_transition_(0),
          out_transition_(0)
{
    // read(in);
}

template<class Alphabet_T>
State<Alphabet_T>::State(int index, const Profile<Alphabet_T>& profile)
        : index_(index),
          ContextProfile<Alphabet_T>(profile)
{}

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
    out << "index\t\t" << index_ << std::endl;
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
