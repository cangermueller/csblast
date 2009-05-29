// Copyright 2009, Andreas Biegert

#ifndef SRC_INITIALIZER_INL_H_
#define SRC_INITIALIZER_INL_H_

#include "initializer.h"

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
#include "utils-inl.h"

namespace cs {

template< class Alphabet, template<class> class Graph >
void SamplingStateInitializer<Alphabet>::Init(HMM<Alphabet>& hmm) const {
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
void LibraryStateInitializer<Alphabet>::Init(HMM<Alphabet>& hmm) const {
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
  hmm.states_logspace_ = lib_->logspace();
  hmm.TransformStatesToLinSpace();

  if (!hmm.full())
    throw Exception("Could not fully initialize all %i HMM states. "
                    "Context library contains too few profiles!",
                    hmm.num_states());

  LOG(DEBUG) << "HMM after state initialization:";
  LOG(DEBUG) << hmm;
}

template<class Alphabet>
void CoEmissionTransitionInitializer<Alphabet>::Init(HMM<Alphabet>& hmm) const {
  const int ncols = hmm.num_cols() - 1;

  for (int k = 0; k < hmm.num_states(); ++k) {
    for (int l = 0; l < hmm.num_states(); ++l) {
      float score = co_emission_(hmm[k], hmm[l], 1, 0, ncols);
      if (score > score_thresh_)
        hmm(k,l) = score - score_thresh_;
    }
  }
  NormalizeTransitions(hmm);
}

// Normalizes transition probabilities to one.
template< class Alphabet,
          template<class A> class Graph >
void NormalizeTransitions(Graph<Alphabet>* graph) {
  const bool logspace = graph->transitions_logspace();
  if (logspace) graph->TransformTransitionsToLinSpace();

  for (int k = 0; k < graph->num_states(); ++k) {
    double sum = 0.0;
    for (int l = 0; l < graph->num_states(); ++l)
      if (graph->test_transition(k,l)) sum += graph->tr(k,l);

    if (sum != 0.0) {
      float fac = 1.0 / sum;
      for (int l = 0; l < graph->num_states(); ++l)
        if (graph->test_transition(k,l))
          graph->tr(k,l) = graph->tr(k,l) * fac;
    }
  }

  if (logspace) graph->TransformTransitionsToLogSpace();
}

}  // namespace cs

#endif  // SRC_INITIALIZER_INL_H_
