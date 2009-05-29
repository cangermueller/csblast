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

template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs) {
  return lhs->prior() > rhs->prior();
}

template< class Alphabet, template<class> class Graph >
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


template< class Alphabet, template<class> class Graph >
SamplingStateInitializer<Alphabet>::SamplingStateInitializer(
    profile_vector profiles,
    float sample_rate,
    const Pseudocounts<Alphabet>* pc = NULL,
    float pc_admixture = 1.0f)
    : profiles_(profiles),
      sample_rate_(sample_rate),
      pc_(pc),
      pc_admixture_(pc_admixture) {
  random_shuffle(profiles_.begin(), profiles_.end());
}

template< class Alphabet, template<class> class Graph >
void SamplingStateInitializer<Alphabet>::Init(Graph<Alphabet>& graph) const {
  // Iterate over randomly shuffled profiles; from each profile we sample a
  // fraction of profile windows.
  for (profile_iterator pi = profiles_.begin();
       pi != profiles_.end() && !graph.full(); ++pi) {
    if ((*pi)->num_cols() < graph.num_cols()) continue;

    // Prepare sample of indices
    std::vector<int> idx;
    for (int i = 0; i <= (*pi)->num_cols() - graph.num_cols(); ++i)
      idx.push_back(i);

    random_shuffle(idx.begin(), idx.end());

    const int sample_size = iround(sample_rate_ * idx.size());
    // sample only a fraction of the profile indices.
    idx.erase(idx.begin() + sample_size, idx.end());

    // Add sub-profiles at sampled indices to graph
    for (std::vector<int>::const_iterator i = idx.begin();
         i != idx.end() && !graph.full(); ++i) {
      CountProfile<Alphabet> p(**pi, *i, graph.num_cols());
      if (pc_) pc_->add_to_profile(ConstantAdmixture(pc_admixture_), &p);
      graph.AddState(p);
    }
  }
  if (!graph.full())
    throw Exception("Could not fully initialize all %i states. "
                    "Maybe too few training profiles provided?",
                    graph.num_states());
}

template< class Alphabet, template<class> class Graph >
LibraryStateInitializer<Alphabet>::LibraryStateInitializer(
    const ProfileLibrary<Alphabet>* lib)  : lib_(lib) {}

template< class Alphabet, template<class> class Graph >
void LibraryStateInitializer<Alphabet>::Init(Graph<Alphabet>& graph) const {
  assert(lib_->num_cols() == graph.num_cols());

  typedef std::vector< shared_ptr< ContextProfile<Alphabet> > > ContextProfiles;
  typedef typename ContextProfiles::const_iterator ContextProfileIter;
  ContextProfiles profiles(lib_->begin(), lib_->end());
  sort(profiles.begin(), profiles.end(), PriorCompare<Alphabet>);

  for (ContextProfileIter it = profiles.begin(); it != profiles.end() &&
         !graph.full(); ++it) {
    graph.AddState(**it);
  }
  graph.states_logspace_ = lib_->logspace();
  graph.TransformStatesToLinSpace();

  if (!graph.full())
    throw Exception("Could not fully initialize all %i states. "
                    "Context library contains too few profiles!",
                    graph.num_states());
}


template< class Alphabet, template<class> class Graph >
void HomogeneousTransitionInitializer<Alphabet>::Init(Graph<Alphabet>& graph) const {
  float w = 1.0f / graph.num_states();
  for (int k = 0; k < graph.num_states(); ++k) {
    for (int l = 0; l < graph.num_states(); ++l) {
      graph(k,l) = w;
    }
  }
}

template< class Alphabet, template<class> class Graph >
void RandomTransitionInitializer<Alphabet>::Init(Graph<Alphabet>& graph) const {
  srand(static_cast<unsigned int>(clock()));
  for (int k = 0; k < graph.num_states(); ++k)
    for (int l = 0; l < graph.num_states(); ++l)
      graph(k,l) =
        static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
  NormalizeTransitions(&graph);
}

template< class Alphabet, template<class> class Graph >
CoEmissionTransitionInitializer<Alphabet>::CoEmissionTransitionInitializer(
    const SubstitutionMatrix<Alphabet>* sm, float score_thresh)
    : co_emission_(sm), score_thresh_(score_thresh) {}

template< class Alphabet, template<class> class Graph >
void CoEmissionTransitionInitializer<Alphabet>::Init(Graph<Alphabet>& graph) const {
  const int ncols = graph.num_cols() - 1;

  for (int k = 0; k < graph.num_states(); ++k) {
    for (int l = 0; l < graph.num_states(); ++l) {
      float score = co_emission_(graph[k], graph[l], 1, 0, ncols);
      if (score > score_thresh_)
        graph(k,l) = score - score_thresh_;
    }
  }
  NormalizeTransitions(graph);
}

}  // namespace cs

#endif  // SRC_INITIALIZER_INL_H_
