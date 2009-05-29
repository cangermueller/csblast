// Copyright 2009, Andreas Biegert

#ifndef SRC_CLUSTERING_H_
#define SRC_CLUSTERING_H_

#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <valarray>
#include <vector>

#include "context_profile-inl.h"
#include "mult_emission-inl.h"
#include "expectation_maximization-inl.h"
#include "profile_library-inl.h"
#include "progress_table.h"
#include "shared_ptr.h"

namespace cs {

// Parameter wrapper for expectation-maximization clustering.
struct ClusteringOptions : public ExpectationMaximizationOptions {
  ClusteringOptions()
      : ExpectationMaximizationOptions(),
        weight_center(1.3f),
        weight_decay(0.9f) {}

  ClusteringOptions(const ClusteringOptions& opts)
      : ExpectationMaximizationOptions(opts),
        weight_center(opts.weight_center),
        weight_decay(opts.weight_decay) {}

  // Weight of central column in multinomial emission
  float weight_center;
  // Exponential decay of window weights
  float weight_decay;
};

// Encapsulation of expectation-maximization clustering.
template< class Alphabet, template<class> class Subject >
class Clustering : public ExpectationMaximization<Alphabet, Subject> {
 public:
  typedef typename std::vector< shared_ptr< Subject<Alphabet> > > data_vector;
  typedef typename
  std::vector< shared_ptr< ContextProfile<Alphabet> > > profiles_vector;

  // Initializes a new clustering object without output.
  Clustering(const ClusteringOptions& opts,
             const data_vector& data,
             ProfileLibrary<Alphabet>& lib);
  // Initializes a new clustering object with output.
  Clustering(const ClusteringOptions& opts,
             const data_vector& data,
             ProfileLibrary<Alphabet>& lib,
             FILE* fout);

  virtual ~Clustering();

 protected:
  // Needed to access names in templatized base class
  using ExpectationMaximization<Alphabet, Subject>::progress_table_;
  using ExpectationMaximization<Alphabet, Subject>::log_likelihood_;
  using ExpectationMaximization<Alphabet, Subject>::num_eff_cols_;
  using ExpectationMaximization<Alphabet, Subject>::epsilon_;
  using ExpectationMaximization<Alphabet, Subject>::data_;

  // Evaluates the responsibilities using the current parameter values.
  virtual void ExpectationStep(const data_vector& block);
  // Reestimate teh parameters using the current responsibilities.
  virtual void MaximizationStep();
  // Prepares all members for clustering.
  virtual void Init();
  // Returns parameter wrapper
  virtual const ClusteringOptions& opts() const { return opts_; }
  // Adds the contribution of the responsibilities for a subject to sufficient
  // statistics for priors.
  void add_contribution_to_priors(const std::valarray<double>& p_zn);
  // Adds the contribution of the responsibilities for a counts profile to
  // sufficient statistics for emissions.
  void add_contribution_to_emissions(const std::valarray<double>& p_zn,
                                     const CountProfile<Alphabet>& c);
  // Adds the contribution of the responsibilities for a sequence to sufficient
  // statistics for emissions.
  void add_contribution_to_emissions(const std::valarray<double>& p_zn,
                                     const Sequence<Alphabet>& s);
  // Updates global sufficient statistics with sufficient statistics calculated
  // on current block.
  void update_sufficient_statistics();

  // Parameter wrapper for clustering.
  const ClusteringOptions& opts_;
  // Profile library with context profiles
  ProfileLibrary<Alphabet> lib_;
  // Calculation of emission probabilities.
  MultEmission<Alphabet> emission_;
  // Global expected sufficient statistics for emissions and state priors.
  profiles_vector profile_stats_;
  // Expected sufficient statistics for emissions and state priors based on
  // current block.
  profiles_vector profile_stats_block_;
};


template< class Alphabet, template<class> class Subject >
class ClusteringProgressTable : public ProgressTable {
 public:
  ClusteringProgressTable(const Clustering<Alphabet, Subject>* clustering,
                          FILE* fout = stdout,
                          int width = 30);

  virtual ~ClusteringProgressTable() {}

  virtual void PrintHeader();
  virtual void PrintRowBegin();
  virtual void PrintRowEnd();

 protected:
  const Clustering<Alphabet, Subject>* clustering_;
};

}  // namespace cs

#endif  // SRC_CLUSTERING_H_
