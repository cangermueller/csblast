// Copyright 2009, Andreas Biegert

#ifndef SRC_CHAIN_GRAPH_INL_H_
#define SRC_CHAIN_GRAPH_INL_H_

#include "chain_graph.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "context_profile-inl.h"
#include "count_profile-inl.h"
#include "profile-inl.h"
#include "profile_library-inl.h"
#include "pseudocounts.h"
#include "log.h"
#include "profile.h"
#include "shared_ptr.h"
#include "substitution_matrix-inl.h"
#include "transition.h"

namespace cs {

template< class Alphabet, template<class> class State >
ChainGraph<Alphabet, State>::ChainGraph()
    : num_states_(0),
      num_cols_(0),
      iterations_(0),
      transitions_logspace_(false) {
  Init();
}

template< class Alphabet, template<class> class State >
ChainGraph<Alphabet, State>::ChainGraph(int num_states, int num_cols)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      transitions_(num_states, num_states),
      transitions_logspace_(false) {
  Init();
}

template< class Alphabet, template<class> class State >
ChainGraph<Alphabet, State>::ChainGraph(FILE* fin)
    : num_states_(0),
      num_cols_(0),
      iterations_(0),
      transitions_logspace_(false) {
  Read(fin);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::Init() {
  states_.reserve(num_states());
  transitions_.resize(num_states(), num_states());
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::InitStates(
    const StateInitializer<Alphabet, State>& st) {
  Clear();
  st.Init(*this);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::InitTransitions(
    const TransitionInitializer<Alphabet, State>& tr) {
  ClearTransitions();
  tr.Init(*this);
}

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::set_transition(int k, int l, float w) {
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

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::erase_transition(int k, int l) {
  transitions_.erase(k,l);
  states_[k]->out_transitions_.erase(l);
  states_[l]->in_transitions_.erase(k);
}

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::Clear() {
  states_.clear();
  transitions_.clear();
  Init();
}

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::ClearTransitions() {
  transitions_.clear();
  for (StateIter si = states_begin(); si != states_end(); ++si)
    (*si)->ClearTransitions();
}

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::TransformTransitionsToLogSpace() {
  if (!transitions_logspace()) {
    for (TransitionIter ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->weight = fast_log2(ti->weight);
    transitions_logspace_ = true;
  }
}

template< class Alphabet, template<class> class State >
inline void ChainGraph<Alphabet, State>::TransformTransitionsToLinSpace() {
  if (transitions_logspace()) {
    for (TransitionIter ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->weight = fast_pow2(ti->weight);
    transitions_logspace_ = false;
  }
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::Read(FILE* fin) {
  LOG(DEBUG1) << "Reading HMM from stream ...";

  ReadHeader(fin);
  Init();
  ReadStates(fin);
  ReadTransitions(fin);

  LOG(DEBUG1) << *this;
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::ReadHeader(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;

  // Check if stream actually contains a serialized HMM
  while (fgetline(buffer, kBufferSize, fin))
    if (strscn(buffer)) break;
  if (!strstr(buffer, class_id()))
    throw Exception("%s does not start with '%s' keyword!", class_id(),
                    class_id());

  // Read number of states
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NSTATES")) {
    ptr = buffer;
    num_states_ = strtoi(ptr);
  } else {
    throw Exception("%s does not contain 'NSTATES' record!", class_id());
  }

  // Read number of transitions
  int ntr = 0;
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NTRANS")) {
    ptr = buffer;
    ntr = strtoi(ptr);
  } else {
    throw Exception("%s does not contain 'NTRANS' record!", class_id());
  }

  // Read number of columns
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "NCOLS")) {
    ptr = buffer;
    num_cols_ = strtoi(ptr);
  } else {
    throw Exception("%s does not contain 'NCOLS' record!", class_id());
  }

  // Read number of iterations
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "ITERS")) {
    ptr = buffer;
    iterations_ = strtoi(ptr);
  } else {
    throw Exception("%s does not contain 'ITERS' record!", class_id());
  }

  // Read transitions logspace
  if (fgetline(buffer, kBufferSize, fin) && strstr(buffer, "TRLOG")) {
    ptr = buffer;
    transitions_logspace_ = strtoi(ptr) == 1;
  } else {
    throw Exception("%s does not contain 'TRLOG' record!", class_id());
  }
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::ReadStates(FILE* fin) {
  while (!full() && !feof(fin)) {
    shared_ptr< State<Alphabet> > state_ptr(new State<Alphabet>(fin));
    states_.push_back(state_ptr);
  }

  if (!full())
    throw Exception("%s has %i states but should have %i!", class_id(),
                    states_.size(), num_states());
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::ReadTransitions(FILE* fin) {
  char buffer[kBufferSize];
  const char* ptr = buffer;
  int k, l;
  float w;

  fgetline(buffer, kBufferSize, fin);  // skip description line
  while (fgetline(buffer, kBufferSize, fin) && buffer[0] != '/' &&
         buffer[1] != '/') {
    ptr = buffer;
    k = strtoi(ptr);
    l = strtoi(ptr);

    if (transitions_logspace())
      w = static_cast<float>(-strtoi_ast(ptr)) / kLogScale;
    else
      w = fast_pow2(static_cast<float>(-strtoi_ast(ptr)) / kLogScale);

    set_transition(k, l, w);
  }
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::Write(FILE* fout) const {
  WriteHeader(fout);
  WriteStates(fout);
  WriteTransitions(fout);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::WriteHeader(FILE* fout) const {
  fprintf(fout, "%s\n", class_id());
  fprintf(fout, "NSTATES\t%i\n", num_states());
  fprintf(fout, "NTRANS\t%i\n", num_transitions());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());
  fprintf(fout, "TRLOG\t%i\n", transitions_logspace() ? 1 : 0);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::WriteStates(FILE* fout) const {
  for (ConstStateIter si = states_begin(); si != states_end(); ++si)
    (*si)->Write(fout);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::WriteTransitions(FILE* fout) const {
  fputs("TRANS\n", fout);

  for (ConstTransitionIter ti = transitions_begin();
       ti != transitions_end(); ++ti) {
    fprintf(fout, "%i\t%i\t",
            static_cast<int>(ti->source), static_cast<int>(ti->target));
    float w = transitions_logspace() ? ti->weight : fast_log2(ti->weight);
    if (w == -INFINITY)
      fputs("*\n", fout);
    else
      fprintf(fout, "%i\n", -iround(w * kLogScale));
  }

  fputs("//\n", fout);
}

template< class Alphabet, template<class> class State >
void ChainGraph<Alphabet, State>::Print(std::ostream& out) const {
  out << class_id() << std::endl;
  out << "Number of states:      " << num_states() << std::endl;
  out << "Number of transitions: " << num_transitions() << std::endl;
  out << "Mean connectivity:     " << strprintf("%-7.1f\n", connectivity());
  out << "Number of columns:     " << num_cols() << std::endl;
  out << "Training iterations:   " << iterations() << std::endl;

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
          out << strprintf("%6.2f  ", 100.0f * fast_pow2(tr(k,l)));
        else
          out << strprintf("%6.2f  ", 100.0f * tr(k,l));
      } else {
        out << "     *  ";
      }
    }
    out << std::endl;
  }
}


template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs) {
  return lhs->prior() > rhs->prior();
}

template< class Alphabet, template<class> class State >
void NormalizeTransitions(ChainGraph<Alphabet, State>* graph, float f) {
  const bool logspace = graph->transitions_logspace();
  if (logspace) graph->TransformTransitionsToLinSpace();

  for (int k = 0; k < graph->num_states(); ++k) {
    double sum = 0.0;
    bool is_zero = true;
    for (int l = 0; l < graph->num_states(); ++l)
      if (graph->test_transition(k,l)) {
        sum += (*graph)(k,l);
        is_zero = false;
      }

    if (!is_zero) {
      float fac = f / sum;
      for (int l = 0; l < graph->num_states(); ++l)
        if (graph->test_transition(k,l)) {
          (*graph)(k,l) = (*graph)(k,l) * fac;
        }
    }
  }

  if (logspace) graph->TransformTransitionsToLogSpace();
}


template< class Alphabet, template<class> class State >
SamplingStateInitializer<Alphabet, State>::SamplingStateInitializer(
    ProfileVec profiles,
    float sample_rate,
    const Pseudocounts<Alphabet>* pc,
    float pc_admixture)
    : profiles_(profiles),
      sample_rate_(sample_rate),
      pc_(pc),
      pc_admixture_(pc_admixture) {
  random_shuffle(profiles_.begin(), profiles_.end());
}

template< class Alphabet, template<class> class State >
void SamplingStateInitializer<Alphabet, State>::Init(
    ChainGraph<Alphabet, State>& graph) const {
  // Iterate over randomly shuffled profiles; from each profile we sample a
  // fraction of profile windows.
  for (ProfileIter pi = profiles_.begin(); pi != profiles_.end() &&
         !graph.full(); ++pi) {
    assert(!(*pi)->logspace());
    if ((*pi)->num_cols() < graph.num_cols()) continue;

    // Prepare sample of indices
    std::vector<int> idx;
    for (int i = 0; i <= (*pi)->num_cols() - graph.num_cols(); ++i)
      idx.push_back(i);
    random_shuffle(idx.begin(), idx.end());
    const int sample_size = iround(sample_rate_ * idx.size());
    idx.erase(idx.begin() + sample_size, idx.end());

    // Add sub-profiles at sampled indices to graph
    for (std::vector<int>::const_iterator i = idx.begin(); i != idx.end() &&
           !graph.full(); ++i) {
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

template< class Alphabet, template<class> class State >
LibraryBasedStateInitializer<Alphabet, State>::LibraryBasedStateInitializer(
    const ProfileLibrary<Alphabet>* lib)  : lib_(lib) {}

template< class Alphabet, template<class> class State >
void LibraryBasedStateInitializer<Alphabet, State>::Init(
    ChainGraph<Alphabet, State>& graph) const {
  assert(lib_->num_cols() == graph.num_cols());
  assert(!lib_->logspace());

  ContextProfileVec profiles(lib_->begin(), lib_->end());
  sort(profiles.begin(), profiles.end(), PriorCompare<Alphabet>);

  for (ContextProfileIter it = profiles.begin(); it != profiles.end() &&
         !graph.full(); ++it) {
    graph.AddState(**it);
  }

  if (!graph.full())
    throw Exception("Could not fully initialize all %i states. "
                    "Context library contains too few profiles!",
                    graph.num_states());
}


template< class Alphabet, template<class> class State >
void HomogeneousTransitionInitializer<Alphabet, State>::Init(
    ChainGraph<Alphabet, State>& graph) const {
  float w = 1.0f / graph.num_states();
  for (int k = 0; k < graph.num_states(); ++k) {
    for (int l = 0; l < graph.num_states(); ++l) {
      graph(k,l) = w;
    }
  }
}

template< class Alphabet, template<class> class State >
void RandomTransitionInitializer<Alphabet, State>::Init(
    ChainGraph<Alphabet, State>& graph) const {
  srand(static_cast<unsigned int>(clock()));

  for (int k = 0; k < graph.num_states(); ++k)
    for (int l = 0; l < graph.num_states(); ++l)
      graph(k,l) = static_cast<float>(rand()) / (1.0f + RAND_MAX);

  NormalizeTransitions(&graph);
}

template< class Alphabet, template<class> class State >
CoEmissionTransitionInitializer<Alphabet, State>::CoEmissionTransitionInitializer(
    const SubstitutionMatrix<Alphabet>* sm, float score_thresh)
    : co_emission_(sm), score_thresh_(score_thresh) {}

template< class Alphabet, template<class> class State >
void CoEmissionTransitionInitializer<Alphabet, State>::Init(
    ChainGraph<Alphabet, State>& graph) const {
  const int ncols = graph.num_cols() - 1;

  for (int k = 0; k < graph.num_states(); ++k) {
    for (int l = 0; l < graph.num_states(); ++l) {
      float score = co_emission_(graph[k], graph[l], 1, 0, ncols);
      if (score > score_thresh_)
        graph(k,l) = score - score_thresh_;
    }
  }

  NormalizeTransitions(&graph);
}

}  // namespace cs

#endif  // SRC_CHAIN_GRAPH_INL_H_
