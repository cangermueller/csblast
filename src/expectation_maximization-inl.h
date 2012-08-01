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

#ifndef CS_EXPECTATION_MAXIMIZATION_INL_H_
#define CS_EXPECTATION_MAXIMIZATION_INL_H_

#include "expectation_maximization.h"

#include "utils.h"

namespace cs {

template< class Alphabet, template<class> class Subject >
ExpectationMaximization<Alphabet, Subject>::ExpectationMaximization(
    const DataVec& data)
    : data_(data),
      num_eff_cols_(0.0f),
      logp_(0.0f),
      logp_prev_(0.0f),
      iterations_(0),
      scan_(1),
      eta_(1.0f) {}

template< class Alphabet, template<class> class Subject >
void ExpectationMaximization<Alphabet, Subject>::Run() {
  SetupBlocks(true);
  if (progress_table_) progress_table_->PrintHeader();

  // Calculate log-likelihood baseline by batch EM without blocks
  if (progress_table_) progress_table_->PrintRowBegin();
  ResetAndAddPseudocounts();
  ExpectationStep(blocks_.front());
  UpdateSufficientStatistics();
  MaximizationStep();
  ++iterations_;
  if (progress_table_) progress_table_->PrintRowEnd();

  // Perform E-step and M-step on each training block until convergence
  SetupBlocks();
  while (!IsDone()) {
    if (static_cast<int>(blocks_.size()) > 1)
      eta_ = opts().eta_null * exp(-opts().beta * (scan_ - 1));
    logp_prev_ = logp_;
    logp_      = 0.0;
    ++scan_;

    if (progress_table_) progress_table_->PrintRowBegin();
    for (int b = 0; b < static_cast<int>(blocks_.size()); ++b) {
      ResetAndAddPseudocounts();
      ExpectationStep(blocks_[b]);
      UpdateSufficientStatistics();
      MaximizationStep();
      ++iterations_;
    }
    if (progress_table_) progress_table_->PrintRowEnd();
    if (eta_ < opts().eta_batch &&
        static_cast<int>(blocks_.size()) > 1) {
      SetupBlocks(true);
      eta_ = 1.0f;
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
      iround(pow(data_.size(), 0.375)) : opts().num_blocks;
    const int block_size = iround(static_cast<float>(data_.size()) / num_blocks);

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
    return fabs(rel_diff()) <= opts().epsilon;
}

}  // namespace cs

#endif  // CS_EXPECTATION_MAXIMIZATION_INL_H_
