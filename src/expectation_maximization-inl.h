// Copyright 2009, Andreas Biegert

#ifndef SRC_EXPECTATION_MAXIMIZATION_INL_H_
#define SRC_EXPECTATION_MAXIMIZATION_INL_H_

#include "expectation_maximization.h"

#include "utils-inl.h"

namespace cs {

template< class Alphabet, template<class> class Subject >
ExpectationMaximization<Alphabet, Subject>::ExpectationMaximization(
    const DataVec& data)
    : data_(data),
      progress_table_(NULL),
      num_eff_cols_(0.0f),
      log_likelihood_(0.0f),
      log_likelihood_prev_(0.0f),
      iterations_(0),
      scan_(1),
      epsilon_(1.0f) {}

template< class Alphabet, template<class> class Subject >
void ExpectationMaximization<Alphabet, Subject>::Run() {
  SetupBlocks(true);
  if (progress_table_) progress_table_->PrintHeader();

  // Calculate log-likelihood baseline by batch EM without blocks
  if (progress_table_) progress_table_->PrintRowBegin();
  ResetAndAddPseudocounts();
  ExpectationStep(blocks_.front());
  MaximizationStep();
  ++iterations_;
  if (progress_table_) progress_table_->PrintRowEnd();

  // Perform E-step and M-step on each training block until convergence
  SetupBlocks();
  while (!IsDone()) {
    if (static_cast<int>(blocks_.size()) > 1)
      epsilon_ = opts().epsilon_null * exp(-opts().beta * (scan_ - 1));
    log_likelihood_prev_ = log_likelihood_;
    log_likelihood_      = 0.0;
    ++scan_;

    if (progress_table_) progress_table_->PrintRowBegin();
    for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
      ResetAndAddPseudocounts();
      ExpectationStep(blocks_[b]);
      MaximizationStep();
      ++iterations_;
    }
    if (progress_table_) progress_table_->PrintRowEnd();
    if (epsilon_ < opts().epsilon_batch &&
        static_cast<int>(blocks_.size()) > 1) {
      SetupBlocks(true);
      epsilon_ = 1.0f;
    }
  }
}

template< class Alphabet, template<class> class Subject >
void ExpectationMaximization<Alphabet, Subject>::SetupBlocks(bool force_batch) {
  if (force_batch && static_cast<int>(blocks_.size()) != 1) {
    blocks_.clear();
    blocks_.push_back(data_);

  } else {
    const int num_blocks = opts().num_blocks == 0 ?
      iround(pow(data_.size(), 3.0/8.0)) : opts().num_blocks;
    const int block_size = iround(data_.size() / num_blocks);

    blocks_.clear();
    for (int b = 0; b < num_blocks; ++b) {
      DataVec block;
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

template< class Alphabet, template<class> class Subject >
inline bool ExpectationMaximization<Alphabet, Subject>::IsDone() const {
  if (scan_ < opts().min_scans)
    return false;
  else if (scan_ >= opts().max_scans)
    return true;
  else
    return fabs(log_likelihood_change()) <= opts().log_likelihood_change;
}

}  // namespace cs

#endif  // SRC_EXPECTATION_MAXIMIZATION_INL_H_
