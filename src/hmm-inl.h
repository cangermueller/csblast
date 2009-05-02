// Copyright 2009, Andreas Biegert

#ifndef SRC_HMM_INL_H_
#define SRC_HMM_INL_H_

#include "hmm.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"

#include "log.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "pseudocounts.h"
#include "shared_ptr.h"
#include "state-inl.h"
#include "transition.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
HMM<Alphabet>::HMM(int num_states, int num_cols)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      states_(),
      transitions_(num_states, num_states),
      transitions_logspace_(false),
      states_logspace_(false)  {
  init();
}

template<class Alphabet>
HMM<Alphabet>::HMM(FILE* fin)
    : num_states_(0),
      num_cols_(0),
      iterations_(0),
      states_(),
      transitions_(),
      transitions_logspace_(false),
      states_logspace_(false) {
  read(fin);
}

template<class Alphabet>
HMM<Alphabet>::HMM(int num_states,
                   int num_cols,
                   const StateInitializer<Alphabet>& st_init,
                   const TransitionInitializer<Alphabet>& tr_init)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      states_(),
      transitions_(num_states, num_states),
      transitions_logspace_(false),
      states_logspace_(false) {
  init();
  st_init.init(*this);
  tr_init.init(*this);
}

template<class Alphabet>
void HMM<Alphabet>::init() {
  states_.reserve(num_states());
  transitions_.resize(num_states(), num_states());
}

template<class Alphabet>
void HMM<Alphabet>::init_states(const StateInitializer<Alphabet>& st_init) {
  clear();
  st_init.init(*this);
}

template<class Alphabet>
void HMM<Alphabet>::init_transitions(
    const TransitionInitializer<Alphabet>& tr_init) {
  clear_transitions();
  tr_init.init(*this);
}

template<class Alphabet>
inline void HMM<Alphabet>::set_transition(int k, int l, float prob) {
  if (transitions_.test(k,l)) {
    // Transitions already set -> modify in place
    Transition* tr = &transitions_[k][l];
    tr->probability = prob;
    AnchoredTransition* out_tr = &states_[k]->out_transitions_[l];
    out_tr->probability = prob;
    AnchoredTransition* in_tr = &states_[l]->in_transitions_[k];
    in_tr->probability = prob;
  } else {
    // Transitions unset -> insert into matrix and tables
    transitions_.set(k, l, Transition(k, l, prob));
    states_[k]->out_transitions_.set(l, AnchoredTransition(l, prob));
    states_[l]->in_transitions_.set(k, AnchoredTransition(k, prob));
  }
}

template<class Alphabet>
inline void HMM<Alphabet>::erase_transition(int k, int l) {
  transitions_.erase(k,l);
  states_[k]->out_transitions_.erase(l);
  states_[l]->in_transitions_.erase(k);
}

template<class Alphabet>
inline void HMM<Alphabet>::clear() {
  states_.clear();
  transitions_.clear();
  init();
}

template<class Alphabet>
void HMM<Alphabet>::clear_transitions() {
  transitions_.clear();
  for (state_iterator si = states_begin(); si != states_end(); ++si)
    (*si)->clear_transitions();
}

template<class Alphabet>
inline int HMM<Alphabet>::AddState(const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("Unable to add state: the HMM contains already %i states!",
                    num_states());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  shared_ptr< State<Alphabet> > state_ptr(new State<Alphabet>(states_.size(),
                                                              profile,
                                                              num_states()));
  state_ptr->set_prior(1.0f / num_states());
  states_.push_back(state_ptr);

  return states_.size() - 1;
}

template<class Alphabet>
inline int HMM<Alphabet>::AddState(const ContextProfile<Alphabet>& profile) {
  if (full())
    throw Exception("Unable to add state: the HMM contains already %i states!",
                    num_states());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  shared_ptr< State<Alphabet> > state(new State<Alphabet>(states_.size(),
                                                          profile,
                                                          num_states()));
  states_.push_back(state);

  return states_.size() - 1;
}

template<class Alphabet>
inline void HMM<Alphabet>::transform_transitions_to_logspace() {
  if (!transitions_logspace()) {
    for (transition_iterator ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->probability = fast_log2(ti->probability);
    transitions_logspace_ = true;
  }
}

template<class Alphabet>
inline void HMM<Alphabet>::transform_transitions_to_linspace() {
  if (transitions_logspace()) {
    for (transition_iterator ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->probability = fast_pow2(ti->probability);
    transitions_logspace_ = false;
  }
}

template<class Alphabet>
inline void HMM<Alphabet>::transform_states_to_logspace() {
  if (!states_logspace()) {
    for (state_iterator si = states_begin(); si != states_end(); ++si)
      (*si)->transform_to_logspace();
    states_logspace_ = true;
  }
}

template<class Alphabet>
inline void HMM<Alphabet>::transform_states_to_linspace() {
  if (states_logspace()) {
    for (state_iterator si = states_begin(); si != states_end(); ++si)
      (*si)->transform_to_linspace();
    states_logspace_ = false;
  }
}

template<class Alphabet>
void HMM<Alphabet>::read(FILE* fin) {
  LOG(DEBUG1) << "Reading HMM from stream ...";

  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Check if stream actually contains a serialized HMM
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, "HMM"))
    throw Exception("HMM does not start with 'HMM' keyword!");

  // Read number of states
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'NSTATES' record!");
  }
  // Read number of transitions
  int ntr = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NTRANS")) {
    ptr = buffer;
    ntr = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'NTRANS' record!");
  }
  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NCOLS")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'NCOLS' record!");
  }
  // Read number of iterations
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ITERS")) {
    ptr = buffer;
    iterations_ = strtoi(ptr);
  } else {
    throw Exception("HMM does not contain 'ITERS' record!");
  }
  // Read transitions logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "TRLOG")) {
    ptr = buffer;
    transitions_logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("hMM does not contain 'TRLOG' record!");
  }
  // Read states logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "STLOG")) {
    ptr = buffer;
    states_logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("Bad format: HMM does not contain 'STLOG' record!");
  }

  init();

  // Read HMM states
  while (!full() && !feof(fin)) {
    shared_ptr< State<Alphabet> > state_ptr(new State<Alphabet>(fin));
    states_.push_back(state_ptr);
  }
  if (!full())
    throw Exception("HMM has %i states but should have %i!",
                    states_.size(), num_states());

  // Read HMM transitions
  int k, l;
  float tr_prob;
  fgetline(buffer, kBufferSize, fin);  // skip description line
  while (fgetline(buffer, kBufferSize, fin)
         && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    k = strtoi(ptr);
    l = strtoi(ptr);
    if (transitions_logspace())
      tr_prob = static_cast<float>(-strtoi_ast(ptr)) / kLogScale;
    else
      tr_prob = fast_pow2(static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
    set_transition(k, l, tr_prob);
  }
  if (num_transitions() != ntr)
    throw Exception("HMM has %i transition records but should have %i!",
                    num_transitions(), ntr);

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void HMM<Alphabet>::write(FILE* fout) const {
  // Write header
  fputs("HMM\n", fout);
  fprintf(fout, "NSTATES\t%i\n", num_states());
  fprintf(fout, "NTRANS\t%i\n", num_transitions());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());
  fprintf(fout, "TRLOG\t%i\n", transitions_logspace() ? 1 : 0);
  fprintf(fout, "STLOG\t%i\n", states_logspace() ? 1 : 0);

  // Write states
  for (const_state_iterator si = states_begin(); si != states_end(); ++si)
    (*si)->write(fout);

  // Write transitions
  fputs("TRANS\n", fout);
  for (const_transition_iterator ti = transitions_begin();
       ti != transitions_end(); ++ti) {
    fprintf(fout, "%i\t%i\t",
            static_cast<int>(ti->from), static_cast<int>(ti->to));
    float log_p =
      transitions_logspace() ? ti->probability : fast_log2(ti->probability);
    if (log_p == -INFINITY)
      fputs("*\n", fout);
    else
      fprintf(fout, "%i\n", -iround(log_p * kLogScale));
  }
  fputs("//\n", fout);
}

template<class Alphabet>
void HMM<Alphabet>::print(std::ostream& out) const {
  out << "HMM" << std::endl;
  out << "Total number of states:      " << num_states() << std::endl;
  out << "Total number of transitions: " << num_transitions() << std::endl;
  out << "Average connectivity:        " << strprintf("%-7.1f", connectivity())
      << std::endl;
  out << "Context profile columns:     " << num_cols() << std::endl;
  out << "Training iterations:         " << iterations() << std::endl;

  for (const_state_iterator si = states_begin(); si != states_end(); ++si)
    out << **si;

  out << "Transition matrix:" << std::endl;
  out << "    ";
  for (int l = 0; l < num_states(); ++l) out << strprintf("%6i  ", l);
  out << std::endl;

  for (int k = 0; k < num_states(); ++k) {
    out << strprintf("%-4i", k);
    for (int l = 0; l < num_states(); ++l) {
      if (test_transition(k,l)) {
        if (transitions_logspace())
          out << strprintf("%6.2f  ",
                           100.0f * fast_pow2(transition_probability(k,l)));
        else
          out << strprintf("%6.2f  ", 100.0f * transition_probability(k,l));
      } else {
        out << "     *  ";
      }
    }
    out << std::endl;
  }
}

// Normalizes transition probabilities to one.
template<class Alphabet>
void normalize_transitions(HMM<Alphabet>& hmm) {
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
      throw Exception("Unable to normalize: state %i has no out-transitions!", k);
    }
  }

  if (logspace) hmm.transform_transitions_to_logspace();
}

// Removes all transitions with probability below or equal to given threshold.
template<class Alphabet>
void sparsify(HMM<Alphabet>& hmm, float threshold) {
  const bool logspace = hmm.transitions_logspace();
  if (logspace) hmm.transform_transitions_to_linspace();

  for (int k = 0; k < hmm.num_states(); ++k)
    for (int l = 0; l < hmm.num_states(); ++l)
      if (hmm.test_transition(k,l) && hmm(k,l) <= threshold)
        hmm.erase_transition(k,l);

  normalize_transitions(hmm);

  if (logspace) hmm.transform_transitions_to_logspace();
}


template<class Alphabet>
void SamplingStateInitializer<Alphabet>::init(HMM<Alphabet>& hmm) const {
  LOG(DEBUG) << "Initializing HMM with " << hmm.num_states()
             << " profile windows randomly sampled from "
             << profiles_.size() << " training profiles ...";

  // Iterate over randomly shuffled profiles; from each profile we sample a
  // fraction of profile windows.
  for (profile_iterator pi = profiles_.begin();
       pi != profiles_.end() && !hmm.full(); ++pi) {
    if ((*pi)->num_cols() < hmm.num_cols()) continue;

    LOG(DEBUG1) << "Processing next training profile ...";
    LOG(DEBUG1) << **pi;

    // Prepare sample of indices
    std::vector<int> idx;
    for (int i = 0; i <= (*pi)->num_cols() - hmm.num_cols(); ++i)
      idx.push_back(i);
    LOG(DEBUG2) << "Available column indices:";
    LOG(DEBUG2) << stringify_container(idx);

    random_shuffle(idx.begin(), idx.end());
    LOG(DEBUG2) << "Shuffled column indices:";
    LOG(DEBUG2) << stringify_container(idx);

    const int sample_size = iround(sample_rate_ * idx.size());
    // sample only a fraction of the profile indices.
    idx.erase(idx.begin() + sample_size, idx.end());
    LOG(DEBUG2) << "Sampled column indicices to be actually used::";
    LOG(DEBUG2) << stringify_container(idx);

    // Add sub-profiles at sampled indices to HMM
    for (std::vector<int>::const_iterator i = idx.begin();
         i != idx.end() && !hmm.full(); ++i) {
      CountProfile<Alphabet> p(**pi, *i, hmm.num_cols());
      LOG(DEBUG1) << "Extracted profile window at position " << *i << ":";
      if (pc_) pc_->add_to_profile(ConstantAdmixture(pc_admixture_), &p);
      hmm.AddState(p);
    }
  }
  if (!hmm.full())
    throw Exception("Could not fully initialize all %i HMM states. "
                    "Maybe too few training profiles provided?",
                    hmm.num_states());

  LOG(DEBUG) << "HMM after state initialization:";
  LOG(DEBUG) << hmm;
}

template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs) {
  return lhs->prior() > rhs->prior();
}

template<class Alphabet>
void LibraryStateInitializer<Alphabet>::init(HMM<Alphabet>& hmm) const {
  assert(lib_->num_cols() == hmm.num_cols());
  LOG(DEBUG) << "Initializing HMM states with profile library ...";

  typedef std::vector< shared_ptr< ContextProfile<Alphabet> > > ContextProfiles;
  typedef typename ContextProfiles::const_iterator ContextProfileIter;
  ContextProfiles profiles(lib_->begin(), lib_->end());
  sort(profiles.begin(), profiles.end(), PriorCompare<Alphabet>);

  for (ContextProfileIter it = profiles.begin(); it != profiles.end() &&
         !hmm.full(); ++it) {
    hmm.AddState(**it);
  }
  hmm.set_states_logspace(lib_->logspace());
  hmm.transform_states_to_linspace();

  if (!hmm.full())
    throw Exception("Could not fully initialize all %i HMM states. "
                    "Context library contains too few profiles!",
                    hmm.num_states());

  LOG(DEBUG) << "HMM after state initialization:";
  LOG(DEBUG) << hmm;
}

}  // namespace cs

#endif  // SRC_HMM_INL_H_
