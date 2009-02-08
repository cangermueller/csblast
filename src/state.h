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
#include "profile.h"
#include "transition.h"

using google::sparsetable;

namespace cs
{

template<class Alphabet_T>
class HMM;

template<class Alphabet_T>
class State : public ContextProfile<Alphabet_T>
{
  public:
    typedef typename sparsetable<AnchoredTransition>::const_nonempty_iterator const_transition_iterator;

    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::read;
    using Profile<Alphabet_T>::write;

    // Constructs dummy state for BEGIN/END state in HMM.
    explicit State(int hmm_size);
    // Constructs HMM state from serialized state read from input stream.
    explicit State(std::istream& in);
    // Constructs HMM state with given profile and all transitions initialized to zero.
    State(int index, const Profile<Alphabet_T>& profile, int hmm_size);

    virtual ~State() {}

    // Returns the index of this state.
    int index() const { return index_; }
    // Returns number of in-transitions.
    int num_in_transitions() const { return in_transitions_.num_nonempty(); }
    // Returns number of out-transitions.
    int num_out_transitions() const { return out_transitions_.num_nonempty(); }
    // Returns the transition probability from this state to state k.
    float to(int k) const;
    // Returns the transition probability from state k to this state.
    float from(int k) const;
    // Clears all in- and out-transitions.
    void clear_transitions();
    // Resizes the transition tables to new HMM size.
    void resize(int hmm_size);

    // Returns a const iterator to start of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_begin() const { return in_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_end() const { return in_transitions_.nonempty_end(); }
    // Returns a const iterator to start of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_begin() const { return out_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_end() const { return out_transitions_.nonempty_end(); }

  protected:
    // HMM needs access to transition tables.
    friend class HMM<Alphabet_T>;

    // Needed to access names in templatized Profile base class
    using Profile<Alphabet_T>::read_header;
    using Profile<Alphabet_T>::read_body;
    using Profile<Alphabet_T>::write_header;
    using Profile<Alphabet_T>::write_body;
    using Profile<Alphabet_T>::print;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Prints the profile in human-readable format to output stream.
    virtual void print(std::ostream& out) const;

  private:
    // Returns serialization class identity.
    virtual const std::string class_identity() const { static std::string id("State"); return id;}

    // Index of state in HMM states vector.
    int index_;
    // Size of HMM to which the state belongs.
    int hmm_size_;
    // List of in-transitions.
    sparsetable<AnchoredTransition> in_transitions_;
    // List of out-transitions.
    sparsetable<AnchoredTransition> out_transitions_;
};  // State



template<class Alphabet_T>
State<Alphabet_T>::State(int hmm_size)
        : ContextProfile<Alphabet_T>(),
          index_(0),
          hmm_size_(hmm_size),
          in_transitions_(hmm_size + 1),
          out_transitions_(hmm_size + 1)
{}

template<class Alphabet_T>
State<Alphabet_T>::State(std::istream& in)
        : ContextProfile<Alphabet_T>(),
          index_(0),
          hmm_size_(0),
          in_transitions_(0),
          out_transitions_(0)
{
    read(in);
}

template<class Alphabet_T>
State<Alphabet_T>::State(int index, const Profile<Alphabet_T>& profile, int hmm_size)
        : ContextProfile<Alphabet_T>(profile),
          index_(index),
          hmm_size_(hmm_size),
          in_transitions_(hmm_size + 1),
          out_transitions_(hmm_size + 1)
{}

template<class Alphabet_T>
float State<Alphabet_T>::to(int k) const
{
    return out_transitions_.test(k) ? out_transitions_.get(k).probability : 0.0f;
}

template<class Alphabet_T>
float State<Alphabet_T>::from(int k) const
{
    return in_transitions_.test(k) ? in_transitions_.get(k).probability : 0.0f;
}

template<class Alphabet_T>
void State<Alphabet_T>::clear_transitions()
{
    in_transitions_.clear();
    out_transitions_.clear();
}

template<class Alphabet_T>
void State<Alphabet_T>::resize(int hmm_size)
{
    clear_transitions();
    in_transitions_.resize(hmm_size + 1);
    out_transitions_.resize(hmm_size + 1);
}

template<class Alphabet_T>
void State<Alphabet_T>::read_header(std::istream& in)
{
    // Read has_counts
    std::string tmp;
    if (getline(in, tmp) && tmp.find("index") != std::string::npos)
        index_ = atoi(tmp.c_str() + 5);
    // Read HMM size
    if (getline(in, tmp) && tmp.find("hmm_size") != std::string::npos)
        hmm_size_ = atoi(tmp.c_str() + 8);
    in_transitions_.resize(hmm_size_ + 1);
    out_transitions_.resize(hmm_size_ + 1);

    Profile<Alphabet_T>::read_header(in);
}

template<class Alphabet_T>
void State<Alphabet_T>::write_header(std::ostream& out) const
{
    out << "index\t\t" << index_ << std::endl;
    out << "hmm_size\t" << hmm_size_ << std::endl;
    Profile<Alphabet_T>::write_header(out);
}

template<class Alphabet_T>
void State<Alphabet_T>::print(std::ostream& out) const
{
    out << "State " << index_ << ":" << std::endl;
    Profile<Alphabet_T>::print(out);
}



// Comparison functor that compares the index of two states.
template<class Alphabet_T>
struct StateIndexCompare : public std::binary_function< State<Alphabet_T>, State<Alphabet_T>, bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.index() < rhs.index();
    }
};

// Comparison functor that compares the number of in-transitions of twho states.
template<class Alphabet_T>
struct NumInTransitionsCompare : public std::binary_function< State<Alphabet_T>, State<Alphabet_T>, bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.num_in_transitions() < rhs.num_in_transitions();
    }
};

// Comparison functor that compares the number of out-transitions of twho states.
template<class Alphabet_T>
struct NumOutTransitionsCompare : public std::binary_function< State<Alphabet_T>, State<Alphabet_T>, bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.num_out_transitions() < rhs.num_out_transitions();
    }
};

}  // cs

#endif
