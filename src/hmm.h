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
    typedef typename std::vector< shared_ptr< State<Alphabet_T> > >::const_iterator const_state_iterator;
    typedef typename sparse_matrix<Transition>::nonempty_iterator transition_iterator;
    typedef typename sparse_matrix<Transition>::const_nonempty_iterator const_transition_iterator;

    // Constructs an empty HMM of given size without any states or transitions.
    HMM(int size);
    // Constructs context HMM from serialized HMM read from input stream.
    HMM(std::istream& in);
    // Constructs context HMM with the help of a state- and a transition-initializer.
    HMM(int size, const StateInitializer<Alphabet_T>& st_init, const TransitionInitializer<Alphabet_T>& tr_init);

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
    void set_transition_probability(int k, int l, float prob);
    // Removes the transition between state k and state l from the HMM.
    void erase_transition(int k, int l);
     // Returns true if there is a transition between state k and state l.
    bool test_transition(int k, int l) { return transitions_.test(k,l); }
    // Clears all states and transitions.
    void clear();
    // Clears all transitions but leaves profile of states untouched.
    void clear_transitions();
    // Adds the given profile as state to the HMM and returns its state index. Note that number of profile columns must be odd!
    int add_profile(const Profile<Alphabet_T>& profile);
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
    // Transforms state profiles and transitions to logspace
    void transform_to_logspace();
    // Transforms state profiles and transitions to linspace
    void transform_to_linspace();

    // Prints HMM in human-readable format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const HMM& hmm)
    {
        hmm.write(out);
        return out;
    }

  private:
    // Scaling factor for serialization of profile log values
    static const int SCALE_FACTOR = 1000;

    // Prints the HMM in human-readable format to output stream. TODO!!!
    // virtual void print(std::ostream& out) const;

    // Initializes the HMM from a serialized HMM read from stream.
    void read(std::istream& in);
    // Initializes the HMM with the mandatory BEGIN/END state.
    void init();

    // Number states in the fully assembled HMM (excl. BEGIN/END state)
    int size_;
    // HMM states ordered by index (1, 2, ..., size)
    std::vector< shared_ptr< State<Alphabet_T> > > states_;
    // Sparse matrix with state transitions
    sparse_matrix<Transition> transitions_;
    // Flag indicating if HMM is in log- or linspace
    bool logspace_;
};  // HMM



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
HMM<Alphabet_T>::HMM(int size,
                     const StateInitializer<Alphabet_T>& st_init,
                     const TransitionInitializer<Alphabet_T>& tr_init)
        : size_(size),
          states_(),
          transitions_(size + 1, size + 1),
          logspace_(false)
{
    init();
    st_init.init(*this);
    tr_init.init(*this);
}

template<class Alphabet_T>
void HMM<Alphabet_T>::init()
{
    states_.reserve(size() + 1);
    states_.push_back( shared_ptr< State<Alphabet_T> >(new State<Alphabet_T>(size())) );
    transitions_.resize(size() + 1, size() + 1);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::set_transition_probability(int k, int l, float prob)
{
    const Transition* tr_ptr = &transitions_.set(k, l, Transition(k, l, prob));
    states_[k]->out_transitions_.set(l, tr_ptr);
    states_[l]->in_transitions_.set(k, tr_ptr);
}

template<class Alphabet_T>
inline void HMM<Alphabet_T>::erase_transition(int k, int l)
{
    transitions_.erase(k,l);
    states_[k]->out_transitions_.erase(l);
    states_[l]->in_transitions_.erase(k);
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
inline int HMM<Alphabet_T>::add_profile(const Profile<Alphabet_T>& profile)
{
    if (num_states() >= size())
        throw Exception("Unable to add state: the HMM contains already %i states!", size());

    shared_ptr< State<Alphabet_T> > state_ptr(new State<Alphabet_T>(num_states() + 1,
                                                                    profile,
                                                                    size()));
    states_.push_back(state_ptr);
    return num_states();
}

template<class Alphabet_T>
void HMM<Alphabet_T>::transform_to_logspace()
{
    for (const_state_iterator si = states_begin(); si != states_end(); ++si)
        (*si)->transform_to_logspace();

    for (transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti) {
        float prob = ti->probability;
        ti->probability = prob == 0.0f ? -std::numeric_limits<float>::infinity() : log2(prob);
    }

    logspace_ = true;
}

template<class Alphabet_T>
void HMM<Alphabet_T>::transform_to_linspace()
{
    for (const_state_iterator si = states_begin(); si != states_end(); ++si)
        (*si)->transform_to_linspace();

    for (transition_iterator ti = transitions_begin(); ti != transitions_end(); ++ti) {
        ti->probability = pow(2.0, ti->probability);
    }

    logspace_ = false;
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
    if (getline(in, tmp) && tmp.find("size") != std::string::npos)
        size_ = atoi(tmp.c_str() + 4);
    else
        throw Exception("Bad format: serialized profile does not contain 'size' record!");

    // Read logspace
    if (getline(in, tmp) && tmp.find("logspace") != std::string::npos)
        logspace_ = atoi(tmp.c_str() + 8) == 1;

    // Read state records
    init();
    while (num_states() < size() && in.peek() && in.good()) {  // peek first to make sure that we don't read beyond '//'
        shared_ptr< State<Alphabet_T> > state_ptr(new State<Alphabet_T>(in));
        states_.push_back(state_ptr);
    }
    if (num_states() != size())
        throw Exception("Error while reading HMM: number of states is %i but should be %i!", num_states(), size());

    // Read all transitions
    std::vector<std::string> tokens;
    getline(in, tmp);  // skip description line
    while (getline(in, tmp)) {
        if (tmp.empty()) continue;
        if (tmp.length() > 1 && tmp[0] == '/' && tmp[1] == '/') break;

        split(tmp, '\t', tokens);
        float log_p = *tokens.back().begin() == '*' ? std::numeric_limits<int>::max() : atoi(tokens.back().c_str());
        set_transition_probability( atoi(tokens[0].c_str()),
                                    atoi(tokens[1].c_str()),
                                    logspace_ ? -log_p / SCALE_FACTOR : pow(2.0, -log_p / SCALE_FACTOR) );
        tokens.clear();
    }

    LOG(DEBUG1) << *this;
}

template<class Alphabet_T>
void HMM<Alphabet_T>::write(std::ostream& out) const
{
    // Write header
    out << "HMM" << std::endl;
    out << "size\t\t" << size() << std::endl;
    out << "logspace\t" << (logspace() ? 1 : 0) << std::endl;

    // Write states (excl. BEGIN/END state)
    for (const_state_iterator si = ++states_begin(); si != states_end(); ++si)
        (*si)->write(out);

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



// Normalizes transition probabilities to one.
template<class Alphabet_T>
void normalize_transitions(HMM<Alphabet_T>& hmm)
{
    for (int k = 0; k <= hmm.size(); ++k) {
        float sum = 0.0f;
        for (int l = 0; l <= hmm.size(); ++l)
            if (hmm.test_transition(k,l)) sum += hmm.transition_probability(k,l);
        if (sum != 0.0f) {
            float fac = 1.0f / sum;
            for (int l = 0; l <= hmm.size(); ++l)
                if (hmm.test_transition(k,l))
                    hmm.set_transition_probability(k, l, fac * hmm.transition_probability(k,l));
        } else {
            // no out-transitions for this state -> connect to END-state
            hmm.set_transition_probability(k, 0, 1.0f);
        }
    }
}

// Removes all transitions with probability below or equal to given threshold.
template<class Alphabet_T>
void sparsify(HMM<Alphabet_T>& hmm, float threshold)
{
    for (int k = 0; k <= hmm.size(); ++k)
        for (int l = 0; l <= hmm.size(); ++l)
            if (hmm.test_transition(k,l) && hmm.transition_probability(k,l) <= threshold)
                hmm.erase_transition(k,l);
    normalize_transitions(hmm);
}



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
        LOG(INFO) << "Initializing HMM with " << hmm.size() << "profile windows randomly sampled from "
                  << profiles_.size() << " training profiles ...";

        // Iterate over randomly shuffled profiles; from each profile sample a fraction of profile windows until HMM is full.
        typedef typename profile_vector::const_iterator const_profile_iterator;
        for (const_profile_iterator pi = profiles_.begin(); pi != profiles_.end() && hmm.num_states() < hmm.size(); ++pi) {
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
            for (std::vector<int>::const_iterator i = idx.begin(); i != idx.end() && hmm.num_states() < hmm.size(); ++i) {
                CountsProfile<Alphabet_T> p(**pi, *i, num_cols_);
                LOG(DEBUG) << "Extracted profile window at position " << *i << ":";
                p.convert_to_frequencies(); // make sure that profile contains frequencies not counts
                pc_.add_to_profile(p);
                hmm.add_profile(p);
            }
        }
        if (hmm.num_states() < hmm.size())
            throw Exception("Initialized only %i out of %i states in HMM. Maybe too few training profiles provided?",
                            hmm.num_states(), hmm.size());
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
class ConstantTransitionInitializer : public TransitionInitializer<Alphabet_T>
{
  public:
    ConstantTransitionInitializer() {}
    virtual ~ConstantTransitionInitializer() {};

    virtual void init(HMM<Alphabet_T>& hmm) const
    {
        float prob = 1.0f / (hmm.size() + 1);
        for (int k = 0; k <= hmm.size(); ++k)
            for (int l = 0; l <= hmm.size(); ++l)
                hmm.set_transition_probability(k, l, prob);
        normalize_transitions(hmm);
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
        for (int k = 0; k <= hmm.size(); ++k)
            for (int l = 0; l <= hmm.size(); ++l) {
                float r = static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
                hmm.set_transition_probability(k, l, r);
            }
        normalize_transitions(hmm);
    }
};

}  // cs

#endif
