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

//#include "initializer.h"
#include "context_profile.h"
#include "exception.h"
#include "context_profile.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"
#include "state.h"
#include "transition.h"
#include "util.h"

namespace cs
{

template<class Alphabet_T>
class HMM
{
  public:
    // Public typedefs
    typedef typename std::vector< shared_ptr< State<Alphabet_T> > >::const_iterator const_state_iterator;
    typedef typename sparse_matrix<Transition>::iterator transition_iterator;
    typedef typename sparse_matrix<Transition>::const_iterator const_transition_iterator;

    // BEGIN and END states have both index 0, this is because is all transitions to state 0
    // go to the END state and all transitions from state 0 start from the BEGIN. There are no
    // in-transitions for the BEGIN state and no out-transitions for the END state.
    static const int BEGIN_END_STATE = 0;

    // Constructs an empty HMM of given size without any states or transitions.
    HMM(int size);
    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in);
    // Constructs context HMM with the help of an initializer.
    //HMM(shared_ptr<Initializer> initializer) : alphabet_(initializer.alphabet()) { initializer.init(*this); }

    virtual ~HMM() {}

    // Returns the number of states in the fully assembled HMM (not counting the BEGIN/END state)
    int size() const { return size_; }
    // Returns the current number of states in the HMM (not counting the BEGIN/END state)
    int num_states() const { return states_.size() - 1; }
    // Returns the number of non-null transitions in the HMM.
    int num_transitions() const { return transitions_.num_nonempty(); }
    // Accessor methods for state i, where i is from interval [0,num_states].
    const State<Alphabet_T>& operator[](int i) const { return *states_[i]; }
    State<Alphabet_T>& operator[](int i) { return *states_[i]; }
    // Returns the transition probability from state k to state l.
    float transition_probability(int k, int l) const { return transitions_.get(k,l).probability; }
    // Sets the transition probability from state k to state l.
    void set_transition(int k, int l, float prob);
    // Removes the transition between state k and state l from the HMM.
    void erase_transition(int k, int l);
    // Removes all transitions with probability below or equal to given threshold.
    void sparsify(float threshold);
    // Clears all states and transitions.
    void clear();
    // Clears all transitions but leaves profile of states untouched.
    void clear_transitions();
    // Adds the given profile as state to the HMM and returns its state index. Note that number of profile columns must be odd!
    int add_state(const Profile<Alphabet_T>& profile);
    // Returns a const iterator to a list of pointers of states.
    const_state_iterator states_begin() const { return states_.begin(); }
    // Returns a const iterator pointing past the end of a list of pointers of states.
    const_state_iterator states_end() const { return states_.end(); }
    // Returns an iterator to a list of transitions.
    transition_iterator transitions_begin() { return transitions_.nonempty_begin(); }
    // Returns an iterator pointing past the end of a list of transitions.
    transition_iterator transitions_end() { return transitions_.nonempty_end(); }
    // Returns a const iterator to a list of transitions.
    const_transition_iterator transitions_begin() const { return transitions_.nonempty_begin(); }
    // Returns a const iterator pointing past the end of a list of transitions.
    const_transition_iterator transitions_end() const { return transitions_.nonempty_end(); }
    // Writes the HMM in serialization format to output stream.
    void write(std::ostream& out) const;
    // Returns true if the HMM transitions and state profiles are in logspace
    bool logspace() const { return logspace_; }

    // Prints HMM in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const HMM& hmm)
    {
        hmm.write(out);
        return out;
    }

  private:
    // Scaling factor for serialization of profile log values
    static const int SCALE_FACTOR = 1000;
    // Class identity keyword for serialization.
    static const char CLASS_IDENTITY[];

    // Prints the HMM in human-readable format to output stream. TODO!!!
    // virtual void print(std::ostream& out) const;

    // Initializes the HMM from a serialized HMM read from stream.
    void read(std::istream& in);
    // Initializes the HMM with the mandatory BEGIN/END state.
    void init();

    // Number states in the fully assembled HMM (excl. BEGIN/END state)
    const int size_;
    // HMM states ordered by index (1, 2, ..., size)
    std::vector< shared_ptr< State<Alphabet_T> > > states_;
    // Sparse matrix with state transitions
    sparse_matrix<Transition> transitions_;
    // Flag indicating if HMM is in log- or linspace
    bool logspace_;
    // // Pointer to .
    // shared_ptr<Initializer> initializer_;
};  // HMM


template<class Alphabet_T>
const char HMM<Alphabet_T>::CLASS_IDENTITY[] = "HMM";

template<class Alphabet_T>
HMM<Alphabet_T>::HMM(int size)
  : size_(size),
    states_(),  // we add states with push_back
    transitions_(size + 1, size + 1),
    logspace_(false)
{
    init();
}

template<class Alphabet_T>
HMM<Alphabet_T>::HMM(std::istream& in)
  : size_(0),
    states_(),
    transitions_(),
    logspace_(false)
{
    read(in);
}

template<class Alphabet_T>
void HMM<Alphabet_T>::init()
{
    states_.reserve(size() + 1);
    states_.push_back( shared_ptr< State<Alphabet_T> >(new State<Alphabet_T>(size() + 1)) );
    transitions_.resize(size() + 1, size() + 1);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::set_transition(int k, int l, float prob)
{
    const Transition* tr_ptr = &transitions_.set(k, l, Transition(k, l, prob));
    states_[k].out_transitions_.set(l, tr_ptr);
    states_[l].in_transitions_.set(k, tr_ptr);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::erase_transition(int k, int l)
{
    transitions_.erase(k,l);
    states_[k].out_transitions_.erase(l);
    states_[l].in_transitions_.erase(k);
}

template<class Alphabet_T>
void HMM<Alphabet_T>::sparsify(float threshold)
{
    transitions_.clear();
    for (int k = 0; k < size_; ++k)
        for (int l = 0; l < size_; ++l)
            if (transition_probability(k,l) <= threshold)
                erase_transition(k,l);
}

template<class Alphabet_T>
void HMM<Alphabet_T>::clear()
{
    states_.clear();
    transitions_.clear();
    init();
}

template<class Alphabet_T>
void HMM<Alphabet_T>::clear_transitions()
{
    transitions_.clear();
    for (const_state_iterator si = states_begin(); si != states_end(); ++si)
        si->clear_transitions();
}

template<class Alphabet_T>
int HMM<Alphabet_T>::add_state(const Profile<Alphabet_T>& profile)
{
    if (num_states() >= size())
        throw Exception("Unable to add state: the HMM contains already %i states!", size());

    shared_ptr< Profile<Alphabet_T> > state_ptr(new State<Alphabet_T>(num_states(), profile, size() + 1));
    states_.push_back(state_ptr);
    return num_states();
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
    if (getline(in, tmp) && tmp.find("size") != std::string::npos)
        size_ = atoi(tmp.c_str() + 4);
    else
        throw Exception("Bad format: serialized profile does not contain 'size' record!");

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    // Read state records
    init();
    while (in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr< State<Alphabet_T> > state_ptr(new State<Alphabet_T>(in));
        states_.push_back(state_ptr);
    }

    // Read all transitions
    std::vector<std::string> tokens;
    getline(in, tmp);  // skip description line
    while (getline(in, tmp)) {
        if (tmp.empty()) continue;
        if (tmp.length() > 1 && tmp[0] == '/' && tmp[1] == '/') break;

        split(tmp, '\t', tokens);
        float log_p = *tokens.back().begin() == '*' ? std::numeric_limits<int>::max() : atoi(tokens.back().c_str());
        set_transition( atoi(tokens[0].c_str()),
                        atoi(tokens[1].c_str()),
                        logspace_ ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR) );
        tokens.clear();
    }
    if (num_states() != size())
        throw Exception("Error while reading HMM: number of states is %i but should be %i!", num_states(), size());

    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void HMM<Alphabet_T>::write(std::ostream& out) const
{
    // Write header
    out << CLASS_IDENTITY << std::endl;
    out << "size\t\t" << size() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;

    // Write states (excl. BEGIN/END state)
    for (const_state_iterator si = ++states_begin(); si != states_end(); ++si)
        si->write(out);

    // Write transitions
    out << "transitions" << std::endl;
    for (const_transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti) {
        out << ti->from << "\t" << ti->to << "\t";
        float logval = logspace_ ? ti->probability : log2(ti->probability);
        if (-logval == std::numeric_limits<float>::infinity())
            out << "*" << std::endl;
        else
            out << -iround(logval * SCALE_FACTOR) << std::endl;
    }
    out << "//" << std::endl;
}

}  // cs

#endif
