// Copyright 2009, Andreas Biegert

#ifndef SRC_INITIALIZER_H_
#define SRC_INITIALIZER_H_

namespace cs {

template< class Alphabet, template<class> class Graph >
class StateInitializer {
 public:
  StateInitializer() {}
  virtual ~StateInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const = 0;
};  // class StateInitializer

template< class Alphabet, template<class> class Graph >
class SamplingStateInitializer : public StateInitializer<Alphabet> {
 public:
  typedef typename
  std::vector< shared_ptr< CountProfile<Alphabet> > > profile_vector;
  typedef typename profile_vector::const_iterator profile_iterator;

  SamplingStateInitializer(profile_vector profiles,
                           float sample_rate,
                           const Pseudocounts<Alphabet>* pc = NULL,
                           float pc_admixture = 1.0f)
      : profiles_(profiles),
        sample_rate_(sample_rate),
        pc_(pc),
        pc_admixture_(pc_admixture) {
    random_shuffle(profiles_.begin(), profiles_.end());
  }

  virtual ~SamplingStateInitializer() {};
  virtual void Init(<Alphabet>& graph) const;

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

// Compare function to sort states in descending prior probability.
template< class Alphabet, template<class> class Graph >
bool PriorCompare(const shared_ptr< ContextProfile<Alphabet> >& lhs,
                  const shared_ptr< ContextProfile<Alphabet> >& rhs);

template< class Alphabet, template<class> class Graph >
class LibraryStateInitializer : public StateInitializer<Alphabet> {
 public:
  LibraryStateInitializer(const ProfileLibrary<Alphabet>* lib)
      : lib_(lib) {}

  virtual ~LibraryStateInitializer() {};
  virtual void Init(Graph<Alphabet>& graph) const;

 private:
  // Profile library of context profiles.
  const ProfileLibrary<Alphabet>* lib_;
};  // class LibraryStateInitializer



template< class Alphabet, template<class> class Graph >
class TransitionInitializer {
 public:
  TransitionInitializer() {}
  virtual ~TransitionInitializer() {}
  virtual void Init(Graph<Alphabet>& graph) const = 0;
};  // class TransitionInitializer

template< class Alphabet, template<class> class Graph >
class HomogeneousTransitionInitializer
    : public TransitionInitializer<Alphabet> {
 public:
  HomogeneousTransitionInitializer() {}
  virtual ~HomogeneousTransitionInitializer() {}

  virtual void Init(Graph<Alphabet>& graph) const {
    float prob = 1.0f / graph.num_states();
    for (int k = 0; k < graph.num_states(); ++k) {
      for (int l = 0; l < graph.num_states(); ++l) {
        graph(k,l) = prob;
      }
    }
  }
};  // class HomogeneousTransitionInitializer

template< class Alphabet, template<class> class Graph >
class RandomTransitionInitializer : public TransitionInitializer<Alphabet> {
 public:
  RandomTransitionInitializer() {}
  virtual ~RandomTransitionInitializer() {}

  virtual void Init(Graph<Alphabet>& graph) const {
    srand(static_cast<unsigned int>(clock()));
    for (int k = 0; k < graph.num_states(); ++k)
      for (int l = 0; l < graph.num_states(); ++l)
        graph(k,l) =
          static_cast<float>(rand()) / (static_cast<float>(RAND_MAX) + 1.0f);
    NormalizeTransitions(&graph);
  }
};  // class RandomTransitionInitializer

template< class Alphabet, template<class> class Graph >
class CoEmissionTransitionInitializer
    : public TransitionInitializer<Alphabet> {
 public:
  CoEmissionTransitionInitializer(const SubstitutionMatrix<Alphabet>* sm,
                                  float score_thresh)
      : co_emission_(sm), score_thresh_(score_thresh) {}
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

}  // namespace cs

#endif  // SRC_INITIALIZER_H_
