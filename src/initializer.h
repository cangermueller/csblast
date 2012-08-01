/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_INITIALIZER_H_
#define CS_INITIALIZER_H_

#include <vector>

#include "profile_library.h"
#include "pseudocounts-inl.h"

namespace cs {

template< class Alphabet, template<class> class Graph >
class StateInitializer {
 public:
  StateInitializer() {}
  virtual ~StateInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const = 0;
};  // class StateInitializer

template< class Alphabet, template<class> class Graph >
class SamplingStateInitializer : public StateInitializer<Alphabet, Graph> {
 public:
  typedef typename
  std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;

  SamplingStateInitializer(profile_vector profiles,
                           float sample_rate,
                           const Pseudocounts<Alphabet>* pc = NULL,
                           float pc_admixture = 1.0f);
  virtual ~SamplingStateInitializer() {};
  virtual void Init(Graph<Alphabet>& graph) const;

 private:
  // Pool of full length sequence profiles to sample from.
  profile_vector profiles_;
  // Fraction of profile windows sampled from each subject.
  float sample_rate_;
  // Pseudocount factory for state profiles.
  const Pseudocounts<Alphabet>* pc_;
  // Constant pseudocount admixture for state profiles.
  float pc_admixture_;
};  // class SamplingStateInitializer

template< class Alphabet, template<class> class Graph >
class LibraryBasedStateInitializer : public StateInitializer<Alphabet, Graph> {
 public:
  LibraryBasedStateInitializer(const ProfileLibrary<Alphabet>* lib);
  virtual ~LibraryBasedStateInitializer() {};
  virtual void Init(Graph<Alphabet>& graph) const;

 private:
  // Profile library of context profiles.
  const ProfileLibrary<Alphabet>* lib_;
};  // class LibraryBasedStateInitializer


template< class Alphabet, template<class> class Graph >
class TransitionInitializer {
 public:
  TransitionInitializer() {}
  virtual ~TransitionInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const = 0;
};  // class TransitionInitializer

template< class Alphabet, template<class> class Graph >
class HomogeneousTransitionInitializer
    : public TransitionInitializer<Alphabet, Graph> {
 public:
  HomogeneousTransitionInitializer() {}
  virtual ~HomogeneousTransitionInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const;
};  // class HomogeneousTransitionInitializer

template< class Alphabet, template<class> class Graph >
class RandomTransitionInitializer : public TransitionInitializer<Alphabet, Graph> {
 public:
  RandomTransitionInitializer() {}
  virtual ~RandomTransitionInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const;
};  // class RandomTransitionInitializer

template< class Alphabet, template<class> class Graph >
class CoEmissionTransitionInitializer
    : public TransitionInitializer<Alphabet, Graph> {
 public:
  CoEmissionTransitionInitializer(const SubstitutionMatrix<Alphabet>* sm,
                                  float score_thresh);
  virtual ~CoEmissionTransitionInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const;

 private:
  // Function object for calculation of co-emission scores
  CoEmission<Alphabet> co_emission_;
  // Minimal co-emission score for inclusion in transition set
  float score_thresh_;
}; // class CoEmissionTransitionInitializer


// Normalizes transition probabilities to one.
template< class Alphabet, template<class> class Graph >
void NormalizeTransitions(Graph<Alphabet>* graph);

// Compare function to sort profiles in descending prior probability.
template<class Alphabet>
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs);

}  // namespace cs

#endif  // CS_INITIALIZER_H_
