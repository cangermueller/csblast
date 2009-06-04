// Copyright 2009, Andreas Biegert

#ifndef SRC_FACTOR_GRAPH_INL_H_
#define SRC_FACTOR_GRAPH_INL_H_

#include "factor_graph.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <iostream>
#include <vector>

#include "exception.h"
#include "initializer-inl.h"
#include "log.h"
#include "profile.h"
#include "shared_ptr.h"

namespace cs {

template< class Alphabet, template<class> class State >
FactorGraph<Alphabet, State>::FactorGraph(int num_states, int num_cols)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      states_(),
      transitions_(num_states, num_states),
      transitions_logspace_(false) {
  Init();
}

template< class Alphabet, template<class> class State >
FactorGraph<Alphabet, State>::FactorGraph(FILE* fin)
    : num_states_(0),
      num_cols_(0),
      iterations_(0),
      transitions_logspace_(false) {
  Read(fin);
}

template< class Alphabet, template<class> class State >
FactorGraph<Alphabet, State>::FactorGraph(
    int num_states,
    int num_cols,
    const StateInitializer<Alphabet, State>& st_init,
    const TransitionInitializer<Alphabet, State>& tr_init)
    : num_states_(num_states),
      num_cols_(num_cols),
      iterations_(0),
      transitions_(num_states, num_states),
      transitions_logspace_(false) {
  Init();
  st_init.Init(*this);
  tr_init.Init(*this);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::Init() {
  states_.reserve(num_states());
  transitions_.resize(num_states(), num_states());
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::InitStates(
    const StateInitializer<Alphabet>& st_init) {
  Clear();
  st_init.Init(*this);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::InitTransitions(
    const TransitionInitializer<Alphabet>& tr_init) {
  ClearTransitions();
  tr_init.Init(*this);
}

template< class Alphabet, template<class> class State >
inline void FactorGraph<Alphabet, State>::set_transition(int k, int l, float w) {
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
inline void FactorGraph<Alphabet, State>::erase_transition(int k, int l) {
  transitions_.erase(k,l);
  states_[k]->out_transitions_.erase(l);
  states_[l]->in_transitions_.erase(k);
}

template< class Alphabet, template<class> class State >
inline void FactorGraph<Alphabet, State>::Clear() {
  states_.clear();
  transitions_.clear();
  Init();
}

template< class Alphabet, template<class> class State >
inline void FactorGraph<Alphabet, State>::ClearTransitions() {
  transitions_.clear();
  for (StateIter si = states_begin(); si != states_end(); ++si)
    (*si)->ClearTransitions();
}

template< class Alphabet, template<class> class State >
inline void FactorGraph<Alphabet, State>::TransformTransitionsToLogSpace() {
  if (!transitions_logspace()) {
    for (TransitionIter ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->weight = fast_log2(ti->weight);
    transitions_logspace_ = true;
  }
}

template< class Alphabet, template<class> class State >
inline void FactorGraph<Alphabet, State>::TransformTransitionsToLinSpace() {
  if (transitions_logspace()) {
    for (TransitionIter ti = transitions_begin();
         ti != transitions_end(); ++ti)
      ti->weight = fast_pow2(ti->weight);
    transitions_logspace_ = false;
  }
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::Read(FILE* fin) {
  LOG(DEBUG1) << "Reading HMM from stream ...";

  ReadHeader(fin);
  ReadStates(fin);
  ReadTransitions(fin);

  LOG(DEBUG1) << *this;
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::ReadHeader(FILE* fin) {
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

  Init();
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::ReadStates(FILE* fin) {
  while (!full() && !feof(fin)) {
    shared_ptr< HMMState<Alphabet> > state_ptr(new HMMState<Alphabet>(fin));
    states_.push_back(state_ptr);
  }

  if (!full())
    throw Exception("%s has %i states but should have %i!", class_id(),
                    states_.size(), num_states());
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::ReadTransitions(FILE* fin) {
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

  if (num_transitions() != ntr)
    throw Exception("%s has %i transition records but should have %i!",
                    class_id(), num_transitions(), ntr);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::Write(FILE* fout) const {
  WriteHeader(fout);
  WriteStates(fout);
  WriteTransitions(fout);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::WriteHeader(FILE* fout) const {
  fprintf(fout, "%s\n", class_id());
  fprintf(fout, "NSTATES\t%i\n", num_states());
  fprintf(fout, "NTRANS\t%i\n", num_transitions());
  fprintf(fout, "NCOLS\t%i\n", num_cols());
  fprintf(fout, "ITERS\t%i\n", iterations());
  fprintf(fout, "TRLOG\t%i\n", transitions_logspace() ? 1 : 0);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::WriteStates(FILE* fout) const {
  for (ConstStateIter si = states_begin(); si != states_end(); ++si)
    (*si)->Write(fout);
}

template< class Alphabet, template<class> class State >
void FactorGraph<Alphabet, State>::WriteTransitions(FILE* fout) const {
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
void FactorGraph<Alphabet, State>::Print(std::ostream& out) const {
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

}  // namespace cs

#endif  // SRC_FACTOR_GRAPH_INL_H_
