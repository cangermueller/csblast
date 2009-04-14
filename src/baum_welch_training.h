// Copyright 2009, Andreas Biegert

#ifndef SRC_BAUM_WELCH_TRAINING_H_
#define SRC_BAUM_WELCH_TRAINING_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <vector>

#include "context_profile-inl.h"
#include "emitter-inl.h"
#include "expectation_maximization-inl.h"
#include "forward_backward_algorithm.h"
#include "hmm-inl.h"
#include "progress_table.h"
#include "shared_ptr.h"

namespace cs {

// Forwrad declarations
template< class Alphabet,
          template<class A> class Subject >
class BaumWelchProgressTable;

// Parameter wrapper for Baum-Welch training.
struct BaumWelchOptions : public EmissionOptions,
                         public ExpectationMaximizationOptions {
  BaumWelchOptions()
      : EmissionOptions(),
        ExpectationMaximizationOptions(),
        transition_pseudocounts(1.0f),
        max_connectivity(0) { }

  BaumWelchOptions(const BaumWelchOptions& opts)
      : EmissionOptions(opts),
        ExpectationMaximizationOptions(opts),
        transition_pseudocounts(opts.transition_pseudocounts),
        max_connectivity(opts.max_connectivity) { }

  // Pseudocounts added to transitions (values below 1 enforce sparsity).
  float transition_pseudocounts;
  // Maximum average connectivity for convergence.
  int max_connectivity;
};

// Encapsulation of Baum-Welch training for HMMs.
template< class Alphabet,
          template<class A> class Subject >
class BaumWelchTraining : public ExpectationMaximization<Alphabet, Subject> {
 public:
  typedef typename std::vector< shared_ptr< Subject<Alphabet> > > data_vector;
  typedef typename std::vector< shared_ptr< ContextProfile<Alphabet> > > profiles_vector;
  typedef typename HMM<Alphabet>::const_transition_iterator const_transition_iterator;

  // Needed to access names in templatized base class
  using ExpectationMaximization<Alphabet, Subject>::log_likelihood_change;

  // Initializes a new training object without output.
  BaumWelchTraining(const BaumWelchOptions& opts,
                    const data_vector& data,
                    HMM<Alphabet>& hmm);
  // Initializes a new training object with output.
  BaumWelchTraining(const BaumWelchOptions& opts,
                    const data_vector& data,
                    HMM<Alphabet>& hmm,
                    FILE* fout);

  virtual ~BaumWelchTraining();

 protected:
  // Needed to access names in templatized base class
  using ExpectationMaximization<Alphabet, Subject>::progress_table_;
  using ExpectationMaximization<Alphabet, Subject>::log_likelihood_;
  using ExpectationMaximization<Alphabet, Subject>::num_eff_cols_;
  using ExpectationMaximization<Alphabet, Subject>::epsilon_;
  using ExpectationMaximization<Alphabet, Subject>::data_;
  using ExpectationMaximization<Alphabet, Subject>::scan_;

  // Runs forward backward algorithm on provided data.
  virtual void expectation_step(const data_vector& block);
  // Calculates and assigns new HMM parameters by Maxmimum-Likelihood estimation.
  virtual void maximization_step();
  // Prepares all members for HMM training.
  virtual void init();
  // Returns true if any termination condition is fullfilled.
  virtual bool terminate() const;
  // Returns parameter wrapper
  virtual const BaumWelchOptions& opts() const { return opts_; }
  // Adds the contribution of a subject's forward-backward matrices to prior
  // probabilities of states.
  void add_contribution_to_priors(const ForwardBackwardMatrices& m);
  // Adds the contribution of a subject's forward-backward matrices to
  // transition counts.
  void add_contribution_to_transitions(const ForwardBackwardMatrices& m);
  // Adds the contribution of a count profile's forward-backward matrices to
  // emission counts.
  void add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                     const CountProfile<Alphabet>& c);
  // Adds the contribution of a sequence's forward-backward matrices to emission
  // counts.
  void add_contribution_to_emissions(const ForwardBackwardMatrices& m,
                                     const Sequence<Alphabet>& s);
  // Updates global sufficient statistics with sufficient statistics calculated
  // on current block.
  void update_sufficient_statistics();

  // Parameter wrapper for clustering.
  const BaumWelchOptions& opts_;
  // HMM to be trained
  HMM<Alphabet>& hmm_;
  // Profile matcher for calculation of emission probabilities.
  Emitter<Alphabet> emitter_;
  // Global expected sufficient statistics for transitions
  sparse_matrix<float> transition_stats_;
  // Global expeted sufficient statistics for emissions and state priors
  profiles_vector profile_stats_;
  // Expected sufficient statistics for transitions based on current block
  sparse_matrix<float> transition_stats_block_;
  // Expeted sufficient statistics for emissions and state priors based on
  // current block
  profiles_vector profile_stats_block_;

  friend class BaumWelchProgressTable<Alphabet, Subject>;
};


template< class Alphabet,
          template<class A> class Subject >
class BaumWelchProgressTable : public ProgressTable {
 public:
  BaumWelchProgressTable(const BaumWelchTraining<Alphabet, Subject>* training,
                         FILE* fout = stdout,
                         int width = 30);

  virtual ~BaumWelchProgressTable() { }

  virtual void print_header();
  virtual void print_row_begin();
  virtual void print_row_end();

 protected:
  const BaumWelchTraining<Alphabet, Subject>* training_;
};

}  // namespace cs

#endif  // SRC_BAUM_WELCH_TRAINING_H_
