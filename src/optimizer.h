// Copyright 2009, Andreas Biegert

#ifndef SRC_OPTIMIZER_H_
#define SRC_OPTIMIZER_H_

#include <vector>

#include "progress_table.h"
#include "scoped_ptr.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

// Abstract base class for expectation maximization algorithms.
struct OptimizerOptions {
  OptimizerOptions()
      : min_scans(10),
        max_scans(500),
        log_likelihood_change(2e-4f),
        num_blocks(0),
        eta(0.5f) {}

  OptimizerOptions(const OptimizerOptions& opts)
      : min_scans(opts.min_scans),
        max_scans(opts.max_scans),
        log_likelihood_change(opts.log_likelihood_change),
        num_blocks(opts.num_blocks),
        eta(opts.eta) {}

  // Minimal number of data scans.
  int min_scans;
  // Maximum number of data scans.
  int max_scans;
  // Log-likelihood change per column for convergence.
  float log_likelihood_change;
  // Number of blocks into which the training data are divided (def=N^3/8).
  int num_blocks;
  // Learning rate eta
  float eta;
};

// Encapsulation of CRF parameter optimization by stochastic gradient descent.
template< class Alphabet, template<class> class Subject >
class Optimizer {
 public:
  typedef typename std::vector< shared_ptr< Subject<Alphabet> > > DataVec;
  typedef typename std::vector<DataVec> BlocksVec;
  typedef matrix<double> ContextDerivs;
  typedef matrix<double> PcDerivs;
  typedef typename CRF<Alphabet>::ConstTransitionIter ConstTransitionIter;

  // Initializes a new optimizer WITHOUT output.
  Optimizer(const BaumWelchOptions& opts,
                    const DataVec& data,
                    const SubstitutionMatrix<Alphabet>* sm,
                    HMM<Alphabet>& hmm);
  // Initializes a new optimizer WITH output.
  Optimizer(const BaumWelchOptions& opts,
            const DataVec& data,
            const SubstitutionMatrix<Alphabet>* sm,
            HMM<Alphabet>& hmm,
            FILE* fout);

  ~Optimizer() {}

  // Runs stochastic gradient optimization until termination criterion is met.
  void Run();

 private:
  // Runs sum-product algorithm on each subject in provided block and updates
  // derivatives according to posterior probabilities.
  void CalculateDerivatives(const DataVec& block);
  // Updates parameters in CRF based on derivatives.
  void UpdateCRF();
  // Returns number of current scan.
  int scan() const { return scan_; }
  // Returns number of completed iterations.
  int iterations() const { return iterations_; }
  // Returns number of training blocks in current scan.
  int num_blocks() const { return blocks_.size(); }
  // Returns learning rate eta in current scan.
  float eta() const { return eta_; }
  // Returns most recent conditional log-likelihood value.
  float log_likelihood() const { return log_likelihood_; }
  // Calculates change between the current log-likelihood and the
  // log-likelihood in previous scan.
  float log_likelihood_change() const {
    return log_likelihood_ - log_likelihood_prev_;
  }
  // Initializes members for running the optimization.
  void Init();
  // Returns parameter wrapper
  const OptimizerOptions& opts() const { return opts_; }
  // Returns true if any termination condition is fullfilled.
  virtual bool IsDone() const;
  // Fills the blocks vector with training data.
  void SetupBlocks(bool force_batch = false);

  // Training data (either sequences or counts profiles)
  const DataVec& data_;
  // Progress table object for output.
  scoped_ptr<ProgressTable> progress_table_;
  // Effective number of training columns.
  float num_eff_cols_;
  // Blocks of training data.
  BlocksVec blocks_;
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
};  // class Optimizer

}  // namespace cs

#endif  // SRC_OPTIMIZER_H_
