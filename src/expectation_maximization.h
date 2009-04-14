// Copyright 2009, Andreas Biegert

#ifndef SRC_EXPECTATION_MAXIMIZATION_H_
#define SRC_EXPECTATION_MAXIMIZATION_H_

#include <vector>

#include "progress_table.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

// Abstract base class for expectation maximization algorithms.
struct ExpectationMaximizationOptions {
  ExpectationMaximizationOptions()
      : min_scans(10),
        max_scans(500),
        log_likelihood_change(2e-4f),
        num_blocks(0),
        epsilon_null(0.5f),
        beta(0.2f),
        epsilon_batch(0.05f) { }

  ExpectationMaximizationOptions(const ExpectationMaximizationOptions& opts)
      : min_scans(opts.min_scans),
        max_scans(opts.max_scans),
        log_likelihood_change(opts.log_likelihood_change),
        num_blocks(opts.num_blocks),
        epsilon_null(opts.epsilon_null),
        beta(opts.beta),
        epsilon_batch(opts.epsilon_batch) { }

  // Minimal number of data scans.
  int min_scans;
  // Maximum number of data scans.
  int max_scans;
  // Log-likelihood change per column for convergence.
  float log_likelihood_change;
  // Number of blocks into which the training data are divided
  // (default:  B=N^(3/8)).
  int num_blocks;
  // Initial value for learning rate epsilon (1-epsilon is preserved fraction of
  // sufficient statistics).
  float epsilon_null;
  // Paramter governing exponential decay of epsilon per scan.
  float beta;
  // Learning rate epsilon below which training is switched to batch mode.
  float epsilon_batch;
};

// Abstract base class for expectation maximization algorithms.
template< class Alphabet,
          template<class A> class Subject >
class ExpectationMaximization {
 public:
  typedef typename std::vector< shared_ptr< Subject<Alphabet> > > data_vector;
  typedef typename std::vector<data_vector> blocks_vector;

  // Runs EM algorithm until termination criterion is met.
  void run();
  // Returns number of current scan.
  int scan() const { return scan_; }
  // Returns number of completed EM iterations.
  int iterations() const { return iterations_; }
  // Returns number of training blocks in current scan.
  int num_blocks() const { return blocks_.size(); }
  // Returns learning rate epsilon in current scan.
  float epsilon() const { return epsilon_; }
  // Returns most recent log-likelihood.
  float log_likelihood() const { return log_likelihood_; }
  // Calculates the change between current log-likelihood and the
  // log-likelihood in previous scan.
  float log_likelihood_change() const {
    return log_likelihood_ - log_likelihood_prev_;
  }

 protected:
  // Constructs a new EM object.
  ExpectationMaximization(const data_vector& data);

  virtual ~ExpectationMaximization() { }

  // Evaluates the responsibilities using the current parameter values.
  virtual void expectation_step(const data_vector& block) = 0;
  // Reestimate teh parameters using the current responsibilities.
  virtual void maximization_step() = 0;
  // Initializes members for running the EM algorithm.
  virtual void init() = 0;
  // Returns parameter wrapper
  virtual const ExpectationMaximizationOptions& opts() const = 0;
  // Returns true if any termination condition is fullfilled.
  virtual bool terminate() const;
  // Fills the blocks vector with training data.
  void setup_blocks(bool force_batch = false);

  // Training data (either sequences or counts profiles)
  const data_vector& data_;
  // Progress table object for output.
  ProgressTable* progress_table_;
  // Effective number of training columns.
  int num_eff_cols_;
  // Blocks of training data.
  blocks_vector blocks_;
  // Incremental likelihood in current iteration
  float log_likelihood_;
  // Complete likelihood after last complete scan of training data
  float log_likelihood_prev_;
  // Number of traning iterations performed so far.
  int iterations_;
  // Number of current data scan.
  int scan_;
  // Current learning rate epsilon.
  float epsilon_;
};

}  // namespace cs

#endif  // SRC_EXPECTATION_MAXIMIZATION_H_
