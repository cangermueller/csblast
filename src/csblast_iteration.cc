// Copyright 2009, Andreas Biegert

#include "csblast_iteration.h"

#include <algorithm>
#include <string>
#include <vector>

#include "globals.h"
#include "log.h"

using std::string;
using std::vector;

namespace cs {

CSBlastIteration::CSBlastIteration(int num_iterations)
    : iterations_todo_(num_iterations), iterations_done_(0) {}

bool CSBlastIteration::HasConverged() const {
  if ((previous_data_.empty() && current_data_.empty()) ||
      previous_data_.size() != current_data_.size())
    return false;

  bool retval = true;
  SeqIds::const_iterator end  = previous_data_.end();
  SeqIds::const_iterator prev = previous_data_.begin();
  SeqIds::const_iterator curr = current_data_.begin();
  for (; prev != end; ++prev, ++curr) {
    if (*prev != *curr) {
      retval = false;
      break;
    }
  }

  return retval;
}

void CSBlastIteration::Advance(const BlastHits& hits) {
  previous_data_ = current_data_;
  current_data_.clear();

  for (BlastHits::ConstHitIter it = hits.begin(); it != hits.end(); ++it) {
    current_data_.push_back(it->definition);
  }
  sort(current_data_.begin(), current_data_.end());
  ++iterations_done_;
}

}  // namespace cs
