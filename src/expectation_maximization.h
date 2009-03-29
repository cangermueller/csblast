#ifndef CS_EXPECTATION_MAXIMIZATION_H
#define CS_EXPECTATION_MAXIMIZATION_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// Abstract base class for expectation maximization algorithms.

#include <iostream>
#include <vector>

#include "progress_table.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs
{

struct ExpectationMaximizationParams
{
    ExpectationMaximizationParams()
            : min_scans(10),
              max_scans(500),
              log_likelihood_change(2e-4f),
              num_blocks(0),
              epsilon_null(0.5f),
              beta(0.2f),
              epsilon_batch(0.05f)
    { }

    ExpectationMaximizationParams(const ExpectationMaximizationParams& params)
            : min_scans(params.min_scans),
              max_scans(params.max_scans),
              log_likelihood_change(params.log_likelihood_change),
              num_blocks(params.num_blocks),
              epsilon_null(params.epsilon_null),
              beta(params.beta),
              epsilon_batch(params.epsilon_batch)
    { }

    // Minimal number of data scans.
    int min_scans;
    // Maximum number of data scans.
    int max_scans;
    // Log-likelihood change per column for convergence.
    float log_likelihood_change;
    // Number of blocks into which the training data are divided (default:  B=N^(3/8)).
    int num_blocks;
    // Initial value for learning rate epsilon (1-epsilon is preserved fraction of sufficient statistics).
    float epsilon_null;
    // Paramter governing exponential decay of epsilon per scan.
    float beta;
    // Learning rate epsilon below which training is switched to batch mode.
    float epsilon_batch;
};


template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
class ExpectationMaximization
{
  public:
    typedef typename std::vector< shared_ptr< Subject_T<Alphabet_T> > > data_vector;
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
    // Calculates the change between current log-likelihood and the log-likelihood in previous scan.
    float log_likelihood_change() const;

  protected:
    // Constructs a new EM object.
    ExpectationMaximization(const data_vector& data);

    virtual ~ExpectationMaximization()
    { }

    // Evaluates the responsibilities using the current parameter values.
    virtual void expectation_step(const data_vector& block) = 0;
    // Reestimate teh parameters using the current responsibilities.
    virtual void maximization_step() = 0;
    // Initializes members for running the EM algorithm.
    virtual void init() = 0;
    // Returns parameter wrapper
    virtual const ExpectationMaximizationParams& params() const = 0;
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



template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
ExpectationMaximization<Alphabet_T, Subject_T>::ExpectationMaximization(const data_vector& data)
        : data_(data),
          progress_table_(NULL),
          num_eff_cols_(0),
          log_likelihood_(0.0f),
          log_likelihood_prev_(0.0f),
          iterations_(0),
          scan_(1),
          epsilon_(1.0f)
{ }

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ExpectationMaximization<Alphabet_T, Subject_T>::run()
{
    setup_blocks(true);
    if (progress_table_) progress_table_->print_header();

    // Calculate log-likelihood baseline by batch EM without blocks
    if (progress_table_) progress_table_->print_row_begin();
    expectation_step(blocks_.front());
    maximization_step();
    ++iterations_;
    if (progress_table_) progress_table_->print_row_end();

    // Perform E-step and M-step on each training block until convergence
    setup_blocks();
    while (!terminate()) {
        if (static_cast<int>(blocks_.size()) > 1)
            epsilon_ = params().epsilon_null * exp(-params().beta * (scan_ - 1));
        log_likelihood_prev_ = log_likelihood_;
        log_likelihood_      = 0.0;
        ++scan_;

        if (progress_table_) progress_table_->print_row_begin();
        for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
            expectation_step(blocks_[b]);
            maximization_step();
            ++iterations_;
        }
        if (progress_table_) progress_table_->print_row_end();
        if (epsilon_ < params().epsilon_batch && static_cast<int>(blocks_.size()) > 1) {
            setup_blocks(true);
            epsilon_ = 1.0f;
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
void ExpectationMaximization<Alphabet_T, Subject_T>::setup_blocks(bool force_batch)
{
    if (force_batch && static_cast<int>(blocks_.size()) != 1) {
        blocks_.clear();
        blocks_.push_back(data_);

    } else {
        const int num_blocks = params().num_blocks == 0 ? iround(pow(data_.size(), 3.0/8.0)) : params().num_blocks;
        const int block_size = iround(data_.size() / num_blocks);

        blocks_.clear();
        for (int b = 0; b < num_blocks; ++b) {
            data_vector block;
            if (b == num_blocks - 1)  // last block may differ in block size
                for (int n = b * block_size; n < static_cast<int>(data_.size()); ++n)
                    block.push_back(data_[n]);
            else
                for (int n = b * block_size; n < (b+1) * block_size; ++n)
                    block.push_back(data_[n]);
            blocks_.push_back(block);
        }
    }
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline float ExpectationMaximization<Alphabet_T, Subject_T>::log_likelihood_change() const
{
    return log_likelihood_ - log_likelihood_prev_;
}

template< class Alphabet_T,
          template<class Alphabet_U> class Subject_T >
inline bool ExpectationMaximization<Alphabet_T, Subject_T>::terminate() const
{
    if (scan_ < params().min_scans)
        return false;
    else if (scan_ >= params().max_scans)
        return true;
    else
        return fabs(log_likelihood_change()) <= params().log_likelihood_change;
}

}  // cs

#endif
