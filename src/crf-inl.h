// Copyright 2009, Andreas Biegert

#ifndef SRC_CRF_INL_H_
#define SRC_CRF_INL_H_

#include "crf.h"

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
#include "crf_state-inl.h"
#include "transition.h"
#include "utils-inl.h"

namespace cs {

template<class Alphabet>
CRF<Alphabet>::CRF(int num_states, int num_cols)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      transitions_(num_states, num_states),
      transitions_logspace_(false) {
  Init();
}

template<class Alphabet>
CRF<Alphabet>::CRF(FILE* fin)
    : num_states_(0),
      num_cols_(0),
      iterations_(0),
      transitions_logspace_(false) {
  Read(fin);
}

template<class Alphabet>
CRF<Alphabet>::CRF(int num_states,
                   int num_cols,
                   const CRFStateInitializer<Alphabet>& st_init,
                   const CRFTransitionInitializer<Alphabet>& tr_init)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      transitions_(num_states, num_states),
      transitions_logspace_(false) {
  Init();
  st_init.Init(this);
  tr_init.Init(this);
}

template<class Alphabet>
void CRF<Alphabet>::Init() {
  states_.reserve(num_states());
  transitions_.resize(num_states(), num_states());
}

template<class Alphabet>
void CRF<Alphabet>::InitStates(const CRFStateInitializer<Alphabet>& st_init) {
  Clear();
  st_init.Init(*this);
}

template<class Alphabet>
void CRF<Alphabet>::InitTransitions(
    const CRFTransitionInitializer<Alphabet>& tr_init) {
  ClearTransitions();
  tr_init.Init(*this);
}

template<class Alphabet>
inline void CRF<Alphabet>::set_transition(int k, int l, float w) {
  if (transitions_.test(k,l)) {
    // Transition is already set -> modify in place
    (&transitions_[k][l])->weight = w;
    (&states_[k]->out_transitions_[l])->weight = w;
    (&states_[l]->in_transitions_[k])->weight = w;
  } else {
    // Transitions unset -> insert into matrix and tables
    transitions_.set(k, l, Transition(k, l, w));
    states_[k]->out_transitions_.set(l, AnchoredTransition(l, w));
    states_[l]->in_transitions_.set(k, AnchoredTransition(k, w));
  }
}

template<class Alphabet>
inline void CRF<Alphabet>::erase_transition(int k, int l) {
  transitions_.erase(k,l);
  states_[k]->out_transitions_.erase(l);
  states_[l]->in_transitions_.erase(k);
}

template<class Alphabet>
inline void CRF<Alphabet>::Clear() {
  states_.clear();
  transitions_.clear();
  Init();
}

template<class Alphabet>
void CRF<Alphabet>::ClearTransitions() {
  transitions_.clear();
  for (StateIter si = states_begin(); si != states_end(); ++si)
    (*si)->ClearTransitions();
}

template<class Alphabet>
int CRF<Alphabet>::AddState(const Profile<Alphabet>& profile) {
  if (full())
    throw Exception("Unable to add state: the CRF contains already %i states!",
                    num_states());
  if (profile.num_cols() != num_cols())
    throw Exception("Profile to add as state has %i columns but should have %i!",
                    profile.num_cols(), num_cols());

  shared_ptr< CRFState<Alphabet> > state_ptr(
      new CRFState<Alphabet>(states_.size(), num_states(), profile));
  states_.push_back(state_ptr);

  return states_.size() - 1;
}

template<class Alphabet>
void CRF<Alphabet>::Read(FILE* fin) {
  LOG(DEBUG1) << "Reading CRF from stream ...";

  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Check if stream actually contains a serialized CRF
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, "CRF"))
    throw Exception("CRF does not start with 'CRF' keyword!");

  // Read number of states
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states_ = strtoi(ptr);
  } else {
    throw Exception("CRF does not contain 'NSTATES' record!");
  }
  // Read number of transitions
  int ntr = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NTRANS")) {
    ptr = buffer;
    ntr = strtoi(ptr);
  } else {
    throw Exception("CRF does not contain 'NTRANS' record!");
  }
  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NCOLS")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("CRF does not contain 'NCOLS' record!");
  }
  // Read number of iterations
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ITERS")) {
    ptr = buffer;
    iterations_ = strtoi(ptr);
  } else {
    throw Exception("CRF does not contain 'ITERS' record!");
  }

  Init();

  // Read CRF states
  while (!full() && !feof(fin)) {
    shared_ptr< CRFState<Alphabet> > state_ptr(new CRFState<Alphabet>(fin));
    states_.push_back(state_ptr);
  }
  if (!full())
    throw Exception("CRF has %i states but should have %i!",
                    states_.size(), num_states());

  // Read CRF transitions
  int k, l;
  fgetline(buffer, kBufferSize, fin);  // skip description line
  while (fgetline(buffer, kBufferSize, fin)
         && buffer[0] != '/' && buffer[1] != '/') {
    ptr = buffer;
    k = strtoi(ptr);
    l = strtoi(ptr);
    set_transition(k, l, static_cast<float>(-strtoi_ast(ptr)) / kLogScale);
  }
  if (num_transitions() != ntr)
    throw Exception("CRF has %i transition records but should have %i!",
                    num_transitions(), ntr);

  LOG(DEBUG1) << *this;
}

template<class Alphabet>
void CRF<Alphabet>::Write(FILE* fout) const {
  // Write header
  fputs("CRF\n", fout);
  fprintf(fout, "NSTATES\t%i\n", num_states());
  fprintf(fout, "NTRANS\t%i\n", num_transitions());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());

  // Write states
  for (ConstStateIter si = states_begin(); si != states_end(); ++si)
    (*si)->Write(fout);

  // Write transitions
  fputs("TRANS\n", fout);
  for (ConstTransitionIter ti = transitions_begin();
       ti != transitions_end(); ++ti) {
    fprintf(fout, "%i\t%i\t",
            static_cast<int>(ti->source), static_cast<int>(ti->target));
    float log_p =
      transitions_logspace() ? ti->weight : fast_log2(ti->weight);
    if (log_p == -INFINITY)
      fputs("*\n", fout);
    else
      fprintf(fout, "%i\n", -iround(log_p * kLogScale));
  }
  fputs("//\n", fout);
}

template<class Alphabet>
void CRF<Alphabet>::Print(std::ostream& out) const {
  out << "CRF" << std::endl;
  out << "Total number of states:      " << num_states() << std::endl;
  out << "Total number of transitions: " << num_transitions() << std::endl;
  out << "Average connectivity:        " << strprintf("%-7.1f", connectivity())
      << std::endl;
  out << "Context profile columns:     " << num_cols() << std::endl;
  out << "Training iterations:         " << iterations() << std::endl;

  for (ConstStateIter si = states_begin(); si != states_end(); ++si)
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
                           100.0f * fast_pow2(tr(k,l)));
        else
          out << strprintf("%6.2f  ", 100.0f * tr(k,l));
      } else {
        out << "     *  ";
      }
    }
    out << std::endl;
  }
}

// Normalizes transition probabilities to one.
template<class Alphabet>
void normalize_transitions(CRF<Alphabet>& crf) {
  const bool logspace = crf.transitions_logspace();
  if (logspace) crf.TransformTransitionsToLinSpace();

  for (int k = 0; k < crf.num_states(); ++k) {
    float sum = 0.0f;
    for (int l = 0; l < crf.num_states(); ++l)
      if (crf.test_transition(k,l)) sum += crf(k,l);

    if (sum != 0.0f) {
      float fac = 1.0f / sum;
      for (int l = 0; l < crf.num_states(); ++l)
        if (crf.test_transition(k,l)) crf(k,l) = crf(k,l) * fac;
    } else {
      throw Exception("Unable to normalize: state %i has no out-transitions!", k);
    }
  }

  if (logspace) crf.TransformTransitionsToLogSpace();
}

// Removes all transitions with probability below or equal to given threshold.
template<class Alphabet>
void sparsify(CRF<Alphabet>& crf, float threshold) {
  const bool logspace = crf.transitions_logspace();
  if (logspace) crf.TransformTransitionsToLinSpace();

  for (int k = 0; k < crf.num_states(); ++k)
    for (int l = 0; l < crf.num_states(); ++l)
      if (crf.test_transition(k,l) && crf(k,l) <= threshold)
        crf.erase_transition(k,l);

  normalize_transitions(crf);

  if (logspace) crf.TransformTransitionsToLogSpace();
}


template<class Alphabet>
void SamplingCRFStateInitializer<Alphabet>::Init(CRF<Alphabet>& crf) const {
  LOG(DEBUG) << "Initializing CRF with " << crf.num_states()
             << " profile windows randomly sampled from "
             << profiles_.size() << " training profiles ...";

  // Iterate over randomly shuffled profiles; from each profile we sample a
  // fraction of profile windows.
  for (profile_iterator pi = profiles_.begin();
       pi != profiles_.end() && !crf.full(); ++pi) {
    if ((*pi)->num_cols() < crf.num_cols()) continue;

    LOG(DEBUG1) << "Processing next training profile ...";
    LOG(DEBUG1) << **pi;

    // Prepare sample of indices
    std::vector<int> idx;
    for (int i = 0; i <= (*pi)->num_cols() - crf.num_cols(); ++i)
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

    // Add sub-profiles at sampled indices to CRF
    for (std::vector<int>::const_iterator i = idx.begin();
         i != idx.end() && !crf.full(); ++i) {
      CountProfile<Alphabet> p(**pi, *i, crf.num_cols());
      LOG(DEBUG1) << "Extracted profile window at position " << *i << ":";
      if (pc_) pc_->add_to_profile(ConstantAdmixture(pc_admixture_), &p);
      crf.AddState(p);
    }
  }
  if (!crf.full())
    throw Exception("Could not fully initialize all %i CRF states. "
                    "Maybe too few training profiles provided?",
                    crf.num_states());

  LOG(DEBUG) << "CRF after state initialization:";
  LOG(DEBUG) << crf;
}

template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs) {
  return lhs->prior() > rhs->prior();
}

template<class Alphabet>
void LibraryCRFStateInitializer<Alphabet>::Init(CRF<Alphabet>& crf) const {
  assert(lib_->num_cols() == crf.num_cols());
  LOG(DEBUG) << "Initializing CRF states with profile library ...";

  typedef std::vector< shared_ptr< ContextProfile<Alphabet> > > ContextProfiles;
  typedef typename ContextProfiles::const_iterator ContextProfileIter;
  ContextProfiles profiles(lib_->begin(), lib_->end());
  sort(profiles.begin(), profiles.end(), PriorCompare<Alphabet>);

  for (ContextProfileIter it = profiles.begin(); it != profiles.end() &&
         !crf.full(); ++it) {
    crf.AddState(**it);
  }
  crf.set_states_logspace(lib_->logspace());
  crf.TransformStatesToLinSpace();

  if (!crf.full())
    throw Exception("Could not fully initialize all %i CRF states. "
                    "Context library contains too few profiles!",
                    crf.num_states());

  LOG(DEBUG) << "CRF after state initialization:";
  LOG(DEBUG) << crf;
}

template<class Alphabet>
void CoEmissionCRFTransitionInitializer<Alphabet>::Init(CRF<Alphabet>& crf) const {
  const int ncols = crf.num_cols() - 1;

  for (int k = 0; k < crf.num_states(); ++k) {
    for (int l = 0; l < crf.num_states(); ++l) {
      float score = co_emission_(crf[k], crf[l], 1, 0, ncols);
      if (score > score_thresh_)
        crf(k,l) = score - score_thresh_;
    }
  }
  normalize_transitions(crf);
}

}  // namespace cs

#endif  // SRC_CRF_INL_H_
