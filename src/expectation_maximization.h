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

#ifndef CS_EXPECTATION_MAXIMIZATION_H_
#define CS_EXPECTATION_MAXIMIZATION_H_

#include <algorithm>
#include <vector>

#include "progress_table.h"
#include "scoped_ptr.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs {

// Abstract base class for expectation maximization algorithms.
struct ExpectationMaximizationOptions {
  ExpectationMaximizationOptions()
      : min_scans(10),
        max_scans(500),
        epsilon(1e-4),
        num_blocks(0),
        eta_null(0.5),
        beta(0.2),
        eta_batch(0.05) {}

  ExpectationMaximizationOptions(const ExpectationMaximizationOptions& opts)
      : min_scans(opts.min_scans),
        max_scans(opts.max_scans),
        epsilon(opts.epsilon),
        num_blocks(opts.num_blocks),
        eta_null(opts.eta_null),
        beta(opts.beta),
        eta_batch(opts.eta_batch) {}

  // Minimal number of data scans.
  int min_scans;
  // Maximum number of data scans.
  int max_scans;
  // Log-likelihood change per column for convergence.
  double epsilon;
  // Number of blocks into which the training data are divided (def=N^3/8).
  int num_blocks;
  // Initial value for learning rate epsilon (1-epsilon is preserved fraction of
  // sufficient statistics).
  double eta_null;
  // Paramter governing exponential decay of epsilon per scan.
  double beta;
  // Learning rate epsilon below which training is switched to batch mode.
  double eta_batch;
};

// Abstract base class for expectation maximization algorithms.
template< class Alphabet, template<class> class Subject >
class ExpectationMaximization {
 public:
  typedef std::vector< shared_ptr< Subject<Alphabet> > > DataVec;
  typedef std::vector<DataVec> BlocksVec;

  // Runs EM algorithm until termination criterion is met.
  void Run();
  // Returns number of current scan.
  int scan() const { return scan_; }
  // Returns number of completed EM iterations.
  int iterations() const { return iterations_; }
  // Returns number of training blocks in current scan.
  int num_blocks() const { return blocks_.size(); }
  // Returns learning rate eta in current scan.
  float eta() const { return eta_; }
  // Returns most recent log-likelihood.
  float logp() const { return logp_; }
  // Calculates the relative difference between the  log-likelihood and the
  // log-likelihood in previous scan.
  float rel_diff() const {
    return (logp_ - logp_prev_) / MAX(fabs(logp_), fabs(logp_prev_));
  }

 protected:
  // Constructs a new EM object.
  ExpectationMaximization(const DataVec& data);

  virtual ~ExpectationMaximization() {}

  // Evaluates the responsibilities using the current parameter values.
  virtual void ExpectationStep(const DataVec& block) = 0;
  // Reestimate teh parameters using the current responsibilities.
  virtual void MaximizationStep() = 0;
  // Initializes members for running the EM algorithm.
  virtual void Init() = 0;
  // Returns parameter wrapper
  virtual const ExpectationMaximizationOptions& opts() const = 0;
  // Returns true if any termination condition is fullfilled.
  virtual bool IsDone() const;
  // Fills the blocks vector with training data.
  void SetupBlocks(bool force_batch = false);
  // Sets all statistics variables to their pseudocount values.
  virtual void ResetAndAddPseudocounts() = 0;
  // Updates global sufficient statistics with sufficient statistics calculated
  // on current block.
  virtual void UpdateSufficientStatistics() = 0;

  // Training data (either sequences or counts profiles)
  const DataVec& data_;
  // Progress table object for output.
  scoped_ptr<ProgressTable> progress_table_;
  // Effective number of training columns.
  float num_eff_cols_;
  // Blocks of training data.
  BlocksVec blocks_;
  // Incremental likelihood in current iteration
  double logp_;
  // Complete likelihood after last complete scan of training data
  double logp_prev_;
  // Number of traning iterations performed so far.
  int iterations_;
  // Number of current data scan.
  int scan_;
  // Current learning rate eta.
  double eta_;
};  // class ExpectationMaximization

}  // namespace cs

#endif  // CS_EXPECTATION_MAXIMIZATION_H_
