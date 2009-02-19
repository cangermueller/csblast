#ifndef CS_HMM_H
#define CS_HMM_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A hidden Markov model that stores context information in the form of
// state-specific context profiles and state transition probabilities.

#include <cstdlib>
#include <ctime>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <vector>

#include "exception.h"
#include "profile.h"
#include "context_profile.h"
#include "counts_profile.h"
#include "log.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "sparse_matrix.h"
#include "state.h"
#include "transition.h"
#include "util.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class TransitionAdaptor;

template<class Alphabet_T>
class StateInitializer
{
  public:
    StateInitializer() {}
    virtual ~StateInitializer() {};
    virtual void init(HMM<Alphabet_T>& hmm) const = 0;
};

template<class Alphabet_T>
class TransitionInitializer
{
  public:
    TransitionInitializer() {}
    virtual ~TransitionInitializer() {};
    virtual void init(HMM<Alphabet_T>& hmm) const = 0;
};

template<class Alphabet_T>
class HMM
{
  public:
    // Public typedefs
    typedef typename std::vector< shared_ptr< State<Alphabet_T> > >::iterator state_iterator;
    typedef typename std::vector< shared_ptr< State<Alphabet_T> > >::const_iterator const_state_iterator;
    typedef typename sparse_matrix<Transition>::nonempty_iterator transition_iterator;
    typedef typename sparse_matrix<Transition>::const_nonempty_iterator const_transition_iterator;

    // Constructs an empty HMM of given size without any states or transitions.
    HMM(int num_states);
    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in);
    // Constructs context HMM with the help of a state- and a transition-initializer.
    HMM(int num_states, const StateInitializer<Alphabet_T>& st_init, const TransitionInitializer<Alphabet_T>& tr_init);

    virtual ~HMM() {}

    // Initializes HMM states with the provided initializer (previous states are lost).
    void init_states(const StateInitializer<Alphabet_T>& st_init);
    // Initializes HMM transitions with the provided initializer (previous transitions are lost).
    void init_transitions(const TransitionInitializer<Alphabet_T>& tr_init);
    // Returns true if all states have been fully assembled.
    bool full() const { return static_cast<int>(states_.size()) == num_states_; }
    // Returns the number of states in the HMM
    int num_states() const { return num_states_; }
    // Returns the number of non-null transitions in the HMM.
    int num_transitions() const { return transitions_.num_nonempty(); }
    // Accessor methods for state i, where i is from interval [0,num_states].
    State<Alphabet_T>& operator[](int i) { return *states_[i]; }
    const State<Alphabet_T>& operator[](int i) const { return *states_[i]; }
    // Accessor methods for transition probability (k,l)
    TransitionAdaptor<Alphabet_T> operator() (int k, int l) { return TransitionAdaptor<Alphabet_T>(this, k, l); }
    float operator() (int k, int l) const { return transitions_.get(k,l).probability; }
    // Sets the transition probability from state k to state l.
    void set_transition(int k, int l, float prob);
    // Returns the transition probability from state k to state l.
    float transition_probability(int k, int l) const { return transitions_.get(k,l).probability; }
    // Removes the transition between state k and state l from the HMM.
    void erase_transition(int k, int l);
     // Returns true if there is a transition between state k and state l.
    bool test_transition(int k, int l) const { return transitions_.test(k,l); }
    // Clears all states and transitions.
    void clear();
    // Clears all transitions but leaves profile of states untouched.
    void clear_transitions();
    // Adds the given profile as state to the HMM and returns its state index. Note that number of profile columns must be odd!
    int add_profile(const Profile<Alphabet_T>& profile);
    // Returns an iterator to a list of pointers of states.
    state_iterator states_begin() { return states_.begin(); }
    // Returns an iterator pointing past the end of a list of pointers of states.
    state_iterator states_end() { return states_.end(); }
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
    // Returns true if transitions are in logspace.
    bool transitions_logspace() const { return transitions_logspace_; }
    // Returns true if state profiles are in logspace.
    bool states_logspace() const { return states_logspace_; }
    // Transforms transitions to logspace.
    void transform_transitions_to_logspace();
    // Transforms transitions to linspace.
    void transform_transitions_to_linspace();
    // Transforms state profiles to logspace.
    void transform_states_to_logspace();
    // Transforms state profiles to linspace.
    void transform_states_to_linspace();

    // Prints HMM in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const HMM& hmm)
    {
        hmm.print(out);
        return out;
    }

  private:
    // Scaling factor for serialization of profile log values
    static const int SCALE_FACTOR = 1000;

    // Prints the HMM in human-readable format to output stream. TODO!!!
    void print(std::ostream& out) const;
    // Initializes the HMM from a serialized HMM read from stream.
    void read(std::istream& in);
    // Initializes the HMM.
    void init();

    // Number states in the fully assembled HMM)
    int num_states_;
    // HMM states ordered by index.
    std::vector< shared_ptr< State<Alphabet_T> > > states_;
    // Sparse matrix with state transitions.
    sparse_matrix<Transition> transitions_;
    // Flag indicating if HMM transitions are log- or linspace
    bool transitions_logspace_;
    // Flag indicating if HMM profile probabilities are in log- or linspace
    bool states_logspace_;
};  // HMM



template<class Alphabet_T>
HMM<Alphabet_T>::HMM(int num_states)
        : num_states_(num_states),
          states_(),  // we add states with push_back
          transitions_(num_states, num_states),
          transitions_logspace_(false),
          states_logspace_(false)
{
    init();
}

template<class Alphabet_T>
HMM<Alphabet_T>::HMM(std::istream& in)
        : num_states_(0),
          states_(),
          transitions_(),
          transitions_logspace_(false),
          states_logspace_(false)
{
    read(in);
}

template<class Alphabet_T>
HMM<Alphabet_T>::HMM(int num_states,
                     const StateInitializer<Alphabet_T>& st_init,
                     const TransitionInitializer<Alphabet_T>& tr_init)
        : num_states_(num_states),
          states_(),
          transitions_(num_states, num_states),
          transitions_logspace_(false),
          states_logspace_(false)
{
    init();
    st_init.init(*this);
    tr_init.init(*this);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::init()
{
    states_.reserve(num_states());
    transitions_.resize(num_states(), num_states());
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::init_states(const StateInitializer<Alphabet_T>& st_init)
{
    clear();
    st_init.init(*this);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::init_transitions(const TransitionInitializer<Alphabet_T>& tr_init)
{
    clear_transitions();
    tr_init.init(*this);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::set_transition(int k, int l, float prob)
{
    if (transitions_.test(k,l)) {
        // transitions already set -> modify in place
        Transition* tr = &transitions_[k][l];
        tr->probability = prob;
        AnchoredTransition* out_tr = &states_[k]->out_transitions_[l];
        out_tr->probability = prob;
        AnchoredTransition* in_tr = &states_[l]->in_transitions_[k];
        in_tr->probability = prob;
    } else {
        // transitions unset -> insert into matrix and tables
        transitions_.set(k, l, Transition(k, l, prob));
        states_[k]->out_transitions_.set(l, AnchoredTransition(l, prob));
        states_[l]->in_transitions_.set(k, AnchoredTransition(k, prob));
    }
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::erase_transition(int k, int l)
{
    transitions_.erase(k,l);
    states_[k]->out_transitions_.erase(l);
    states_[l]->in_transitions_.erase(k);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::clear()
{
    states_.clear();
    transitions_.clear();
    init();
}

template<class Alphabet_T>
void HMM<Alphabet_T>::clear_transitions()
{
    transitions_.clear();
    for (state_iterator si = states_begin(); si != states_end(); ++si)
        (*si)->clear_transitions();
}

template<class Alphabet_T>
inline int HMM<Alphabet_T>::add_profile(const Profile<Alphabet_T>& profile)
{
    if (full())
        throw Exception("Unable to add state: the HMM contains already %i states!", num_states());

    shared_ptr< State<Alphabet_T> > state_ptr(new State<Alphabet_T>(states_.size(),
                                                                    profile,
                                                                    num_states()));
    state_ptr->set_prior(1.0f / num_states());  // start with unbiased prior probabilities
    states_.push_back(state_ptr);
    return states_.size();
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::transform_transitions_to_logspace()
{
    if (!transitions_logspace()) {
        for (transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti)
            ti->probability = log2(ti->probability);
        transitions_logspace_ = true;
    }
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::transform_transitions_to_linspace()
{
    if (transitions_logspace()) {
        for (transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti)
            ti->probability = pow(2.0, ti->probability);
        transitions_logspace_ = false;
    }
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::transform_states_to_logspace()
{
    if (!states_logspace()) {
        for (state_iterator si = states_begin(); si != states_end(); ++si)
            (*si)->transform_to_logspace();
        states_logspace_ = true;
    }
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::transform_states_to_linspace()
{
    if (states_logspace()) {
        for (state_iterator si = states_begin(); si != states_end(); ++si)
            (*si)->transform_to_linspace();
        states_logspace_ = false;
    }
}

template<class Alphabet_T>
void HMM<Alphabet_T>::read(std::istream& in)
{
    LOG(DEBUG1) << "Reading HMM from stream ...";

    // Check if stream actually contains a serialized HMM
    std::string tmp;
    while (getline(in, tmp) && tmp.empty()) continue;
    if (tmp.find("HMM") == std::string::npos)
        throw Exception("Bad format: serialized HMM does not start with 'HMM' keyword!");
    // Read number of states
    if (getline(in, tmp) && tmp.find("num_states") != std::string::npos)
        num_states_ = atoi(tmp.c_str() + 10);
    else
        throw Exception("Bad format: serialized profile does not contain 'size' record!");
    // Read number of transitions
    int ntr = 0;
    if (getline(in, tmp) && tmp.find("num_transitions") != std::string::npos)
        ntr = atoi(tmp.c_str() + 15);
    else
        throw Exception("Bad format: serialized profile does not contain 'num_transitions' record!");
    // Read transitions_logspace
    if (getline(in, tmp) && tmp.find("transitions_logspace") != std::string::npos)
        transitions_logspace_ = atoi(tmp.c_str() + 20) == 1;
    // Read states_logspace
    if (getline(in, tmp) && tmp.find("states_logspace") != std::string::npos)
        states_logspace_ = atoi(tmp.c_str() + 15) == 1;
    // Read state records
    init();
    while (!full() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr< State<Alphabet_T> > state_ptr(new State<Alphabet_T>(in));
        states_.push_back(state_ptr);
    }
    if (!full())
        throw Exception("Error while reading HMM: number of state records is %i but should be %i!",
                        states_.size(), num_states());

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
                        transitions_logspace() ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR) );
        tokens.clear();
    }
    if (num_transitions() != ntr)
        throw Exception("Serialized HMM has %i transition records but should have %i!", num_transitions(), ntr);
    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void HMM<Alphabet_T>::write(std::ostream& out) const
{
    // Write header
    out << "HMM" << std::endl;
    out << "num_states\t\t" << num_states() << std::endl;
    out << "num_transitions\t\t" << num_transitions() << std::endl;
    out << "transitions_logspace\t" << (transitions_logspace() ? 1 : 0) << std::endl;
    out << "states_logspace\t\t" << (states_logspace() ? 1 : 0) << std::endl;

    // Write states (excl. BEGIN/END state)
    for (const_state_iterator si = states_begin(); si != states_end(); ++si)
        (*si)->write(out);

    // Write transitions
    out << "transitions" << std::endl;
    for (const_transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti) {
        out << ti->from << "\t" << ti->to << "\t";
        float logval = transitions_logspace() ? ti->probability : log2(ti->probability);
        if (-logval == std::numeric_limits<float>::infinity())
            out << "*" << std::endl;
        else
            out << -iround(logval * SCALE_FACTOR) << std::endl;
    }
    out << "//" << std::endl;
}

template<class Alphabet_T>
void HMM<Alphabet_T>::print(std::ostream& out) const
{
    std::ios_base::fmtflags flags = out.flags();  // save flags

    out << "HMM" << std::endl;
    out << "Total number of states:        " << num_states() << std::endl;
    out << "Total number of transitions:   " << num_transitions() << std::endl;
    out << "Average number of transitions: " << iround(static_cast<float>(num_transitions()) / num_states()) << std::endl;

    for (const_state_iterator si = states_begin(); si != states_end(); ++si)
        (*si)->print(out);

    out << "Transition matrix:" << std::endl;
    out << "    ";
    for (int l = 0; l < num_states(); ++l) out << strprintf("%6i  ", l);
    out << std::endl;

    for (int k = 0; k < num_states(); ++k) {
        out << strprintf("%-4i", k);
        for (int l = 0; l < num_states(); ++l) {
            if (test_transition(k,l))
                out << strprintf("%6.2f  ", 100.0f * ( transitions_logspace() ?
                                                       pow(2.0, transition_probability(k,l)) : transition_probability(k,l) ));
            else
                out << "     *  ";
        }
        out << std::endl;
    }

    out.flags(flags);
}



// Normalizes transition probabilities to one.
template<class Alphabet_T>
void normalize_transitions(HMM<Alphabet_T>& hmm)
{
    const bool logspace = hmm.transitions_logspace();
    if (logspace) hmm.transform_transitions_to_linspace();

    for (int k = 0; k < hmm.num_states(); ++k) {
        float sum = 0.0f;
        for (int l = 0; l < hmm.num_states(); ++l)
            if (hmm.test_transition(k,l)) sum += hmm(k,l);

        if (sum != 0.0f) {
            float fac = 1.0f / sum;
            for (int l = 0; l < hmm.num_states(); ++l)
                if (hmm.test_transition(k,l)) hmm(k,l) = hmm(k,l) * fac;
        } else {
            throw Exception("Unable to normalize HMM transitions: state %i has no out-transitions!", k);
        }
    }

    if (logspace) hmm.transform_transitions_to_logspace();
}

// Removes all transitions with probability below or equal to given threshold.
template<class Alphabet_T>
void sparsify(HMM<Alphabet_T>& hmm, float threshold)
{
    const bool logspace = hmm.transitions_logspace();
    if (logspace) hmm.transform_transitions_to_linspace();

    for (int k = 0; k < hmm.num_states(); ++k)
        for (int l = 0; l < hmm.num_states(); ++l)
            if (hmm.test_transition(k,l) && hmm(k,l) <= threshold)
                hmm.erase_transition(k,l);

    normalize_transitions(hmm);

    if (logspace) hmm.transform_transitions_to_logspace();
}



template<class Alphabet_T>
class TransitionAdaptor {
  public:
    TransitionAdaptor(HMM<Alphabet_T>* hmm, int k, int l)
            : hmm_(hmm), k_(k), l_(l)
    {}

    TransitionAdaptor& operator= (float val)
    {
        hmm_->set_transition(k_, l_, val);
        return *this;
    }

    operator float() { return hmm_->transition_probability(k_, l_); }   // we look like a value

 private:
    HMM<Alphabet_T>* hmm_;
    int k_;
    int l_;
};



template<class Alphabet_T>
class RandomSampleStateInitializer : public StateInitializer<Alphabet_T>
{
  public:
    typedef typename std::vector< shared_ptr< CountsProfile<Alphabet_T> > > profile_vector;

    RandomSampleStateInitializer(profile_vector profiles,
                                 int num_cols,
                                 float sample_rate,
                                 const Pseudocounts<Alphabet_T>& pc)
            : profiles_(profiles),
              num_cols_(num_cols),
              sample_rate_(sample_rate),
              pc_(pc)
    {
        random_shuffle(profiles_.begin(), profiles_.end());
    }

    virtual ~RandomSampleStateInitializer() {};

    virtual void init(HMM<Alphabet_T>& hmm) const
    {
        LOG(INFO) << "Initializing HMM with " << hmm.num_states() << " profile windows randomly sampled from "
                  << profiles_.size() << " training profiles ...";

        // Iterate over randomly shuffled profiles; from each profile sample a fraction of profile windows until HMM is full.
        typedef typename profile_vector::const_iterator const_profile_iterator;
        for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end() && !hmm.full(); ++pi) {
            if ((*pi)->num_cols() < num_cols_) continue;

            LOG(DEBUG) << "Processing next training profile ...";
            LOG(DEBUG) << **pi;

            // Prepare sample of indices
            std::vector<int> idx;
            for (int i = 0; i <= (*pi)->num_cols() - num_cols_; ++i) idx.push_back(i);
            LOG(DEBUG) << "Available column indices:";
            LOG(DEBUG) << stringify_container(idx);

            random_shuffle(idx.begin(), idx.end());
            LOG(DEBUG) << "Shuffled column indices:";
            LOG(DEBUG) << stringify_container(idx);

            const int sample_size = iround(sample_rate_ * idx.size());
            idx.erase(idx.begin() + sample_size, idx.end());  // sample only a fraction of the profile indices.
            LOG(DEBUG) << "Sampled column indicices to be actually used::";
            LOG(DEBUG) << stringify_container(idx);

            // Add sub-profiles at sampled indices to HMM
            for (std::vector<int>::const_iterator i = idx.begin(); i != idx.end() && !hmm.full(); ++i) {
                CountsProfile<Alphabet_T> p(**pi, *i, num_cols_);
                LOG(DEBUG) << "Extracted profile window at position " << *i << ":";
                p.convert_to_frequencies(); // make sure that profile contains frequencies not counts
                pc_.add_to_profile(p);
                std::cerr << hmm.add_profile(p);
            }
        }
        if (!hmm.full())
            throw Exception("Could not fully initialize %i HMM states. Maybe too few training profiles provided?",
                            hmm.num_states());
        LOG(DEBUG) << "HMM after full state assembly:";
        LOG(DEBUG) << hmm;
    }

  private:
    // Pool of full length sequence profiles (possibly with pseudocounts) to be sampled from.
    profile_vector profiles_;
    // Number of columns per state profile.
    const int num_cols_;
    // Fraction of profile windows to be sampled from each training profile.
    const float sample_rate_;
    // Pseudocounts to be added to sampled profiles.
    const Pseudocounts<Alphabet_T>& pc_;
};  // RandomSampleStateInitializer

template<class Alphabet_T>
class HomogeneousTransitionInitializer : public TransitionInitializer<Alphabet_T>
{
  public:
    HomogeneousTransitionInitializer() {}
    virtual ~HomogeneousTransitionInitializer() {};

    virtual void init(HMM<Alphabet_T>& hmm) const
    {
        float prob = 1.0f / hmm.num_states();
        for (int k = 0; k < hmm.num_states(); ++k) {
            for (int l = 0; l < hmm.num_states(); ++l) {
                hmm(k,l) = prob;
            }
        }
    }
};

template<class Alphabet_T>
class RandomTransitionInitializer : public TransitionInitializer<Alphabet_T>
{
  public:
    RandomTransitionInitializer() {}
    virtual ~RandomTransitionInitializer() {};

    virtual void init(HMM<Alphabet_T>& hmm) const
    {
        srand(static_cast<unsigned int>(clock()));
        for (int k = 0; k < hmm.num_states(); ++k)
            for (int l = 0; l < hmm.num_states(); ++l)
                hmm(k,l) = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
        normalize_transitions(hmm);
    }
};

}  // cs

#endif
