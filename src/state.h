#ifndef CS_STATE_H
#define CS_STATE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// The state object of a context HMM.

#include <cstdio>
#include <cstring>

#include <google/sparsetable>

#include "context_profile.h"
#include "profile.h"
#include "transition.h"

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
    typedef typename sparsetable<AnchoredTransition>::const_nonempty_iterator
    const_transition_iterator;

    // Needed to access names in templatized Profile base class
    using ContextProfile<Alphabet_T>::num_cols;
    using ContextProfile<Alphabet_T>::alphabet_size;
    using ContextProfile<Alphabet_T>::logspace;
    using ContextProfile<Alphabet_T>::prior;
    using ContextProfile<Alphabet_T>::set_prior;
    using ContextProfile<Alphabet_T>::index;
    using ContextProfile<Alphabet_T>::set_index;
    using ContextProfile<Alphabet_T>::read;

    // Constructs HMM state from serialized state read from input stream.
    explicit State(std::istream& in);
    // Constructs HMM state from serialized state read from input stream.
    explicit State(FILE* fin);
    // Constructs HMM state with given profile and all transitions initialized to zero.
    State(int index, const Profile<Alphabet_T>& profile, int num_states);

    virtual ~State() { }

    // Returns number of in-transitions.
    int num_in_transitions() const
    { return in_transitions_.num_nonempty(); }
    // Returns number of out-transitions.
    int num_out_transitions() const
    { return out_transitions_.num_nonempty(); }
    // Returns the transition probability from this state to state k.
    float to(int k) const;
    // Returns the transition probability from state k to this state.
    float from(int k) const;
    // Clears all in- and out-transitions.
    void clear_transitions();
    // Resizes the transition tables to new HMM size.
    void resize(int num_states);

    // Returns a const iterator to start of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_begin() const
    { return in_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null in-transition pointers.
    const_transition_iterator in_transitions_end() const
    { return in_transitions_.nonempty_end(); }
    // Returns a const iterator to start of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_begin() const
    { return out_transitions_.nonempty_begin(); }
    // Returns a const iterator past the end of list with non-null out-transition pointers.
    const_transition_iterator out_transitions_end() const
    { return out_transitions_.nonempty_end(); }

  protected:
    // HMM needs access to transition tables.
    friend class HMM<Alphabet_T>;

    // Needed to access names in templatized Profile base class
    using ContextProfile<Alphabet_T>::BUFFER_SIZE;
    using ContextProfile<Alphabet_T>::read_header;
    using ContextProfile<Alphabet_T>::read_body;
    using ContextProfile<Alphabet_T>::write_header;
    using ContextProfile<Alphabet_T>::write_body;

    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(std::istream& in);
    // Reads and initializes serialized scalar data members from stream.
    virtual void read_header(FILE* in);
    // Writes serialized scalar data members to stream.
    virtual void write_header(std::ostream& out) const;
    // Writes serialized scalar data members to stream.
    virtual void write_header(FILE* fout) const;

  private:
    // Class identifier
    static const char* CLASS_ID;

    // Returns serialization class identity.
    virtual const std::string class_identity() const
    { static std::string id("State"); return id; }
    virtual const char* class_id() const { return CLASS_ID; }

    // Size of HMM to which the state belongs.
    int num_states_;
    // List of in-transitions.
    sparsetable<AnchoredTransition> in_transitions_;
    // List of out-transitions.
    sparsetable<AnchoredTransition> out_transitions_;
};  // State



template<class Alphabet_T>
const char* State<Alphabet_T>::CLASS_ID = "State";

template<class Alphabet_T>
inline State<Alphabet_T>::State(std::istream& in)
        : ContextProfile<Alphabet_T>(),
          num_states_(0),
          in_transitions_(0),
          out_transitions_(0)
{
    read(in);
}

template<class Alphabet_T>
inline State<Alphabet_T>::State(FILE* fin)
        : ContextProfile<Alphabet_T>(),
          num_states_(0),
          in_transitions_(0),
          out_transitions_(0)
{
    read(fin);
}

template<class Alphabet_T>
inline State<Alphabet_T>::State(int index, const Profile<Alphabet_T>& profile, int num_states)
        : ContextProfile<Alphabet_T>(index, profile),
          num_states_(num_states),
          in_transitions_(num_states),
          out_transitions_(num_states)
{ }

template<class Alphabet_T>
inline float State<Alphabet_T>::to(int k) const
{
    return out_transitions_.test(k) ? out_transitions_.get(k).probability : 0.0f;
}

template<class Alphabet_T>
inline float State<Alphabet_T>::from(int k) const
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
void State<Alphabet_T>::resize(int num_states)
{
    clear_transitions();
    in_transitions_.resize(num_states);
    out_transitions_.resize(num_states);
}

template<class Alphabet_T>
void State<Alphabet_T>::read_header(std::istream& in)
{
    ContextProfile<Alphabet_T>::read_header(in);

    // Read HMM size
    std::string tmp;
    if (getline(in, tmp) && tmp.find("num_states") != std::string::npos)
        num_states_ = atoi(tmp.c_str() + 10);

    in_transitions_.resize(num_states_);
    out_transitions_.resize(num_states_);
}

template<class Alphabet_T>
void State<Alphabet_T>::read_header(FILE* fin)
{
    ContextProfile<Alphabet_T>::read_header(fin);

    // Read HMM size
    char buffer[BUFFER_SIZE];
    const char* ptr = buffer;
    if (fgetline(buffer, BUFFER_SIZE, fin) && strstr(buffer, "num_states")) {
        ptr = buffer;
        num_states_ = strtoi(ptr);
    } else {
        throw Exception("Bad format: profile does not contain 'num_states' record!");
    }

    in_transitions_.resize(num_states_);
    out_transitions_.resize(num_states_);
}

template<class Alphabet_T>
void State<Alphabet_T>::write_header(std::ostream& out) const
{
    ContextProfile<Alphabet_T>::write_header(out);
    out << "num_states\t" << num_states_ << std::endl;
}

template<class Alphabet_T>
void State<Alphabet_T>::write_header(FILE* fout) const
{
    ContextProfile<Alphabet_T>::write_header(fout);

    fprintf(fout, "num_states\t%i\n", num_states_);
}



// Comparison functor that compares the index of two states.
template<class Alphabet_T>
struct StateIndexCompare : public std::binary_function< State<Alphabet_T>,
                                                        State<Alphabet_T>,
                                                        bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.index() < rhs.index();
    }
};

// Comparison functor that compares the number of in-transitions of twho states.
template<class Alphabet_T>
struct NumInTransitionsCompare : public std::binary_function< State<Alphabet_T>,
                                                              State<Alphabet_T>,
                                                              bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.num_in_transitions() < rhs.num_in_transitions();
    }
};

// Comparison functor that compares the number of out-transitions of twho states.
template<class Alphabet_T>
struct NumOutTransitionsCompare : public std::binary_function< State<Alphabet_T>,
                                                               State<Alphabet_T>,
                                                               bool >
{
    bool operator()(const State<Alphabet_T>& lhs, const State<Alphabet_T>& rhs) const
    {
        return lhs.num_out_transitions() < rhs.num_out_transitions();
    }
};

}  // cs

#endif
