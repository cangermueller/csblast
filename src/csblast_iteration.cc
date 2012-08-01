/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

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

#include "cs.h"
#include "csblast_iteration.h"

using std::string;
using std::vector;

namespace cs {

CSBlastIteration::CSBlastIteration(int niters)
    : iterations_todo_(niters), iterations_done_(0) {}

bool CSBlastIteration::HasConverged() const {
  // For an object that hasn't been 'advanced' or one that only has performed
  // one iteration it doesn't make sense to have converged
  if (iterations_done_ <= 1)
    return false;
  // If the size differs, we obviously have not converged
  if (previous_data_.size() != current_data_.size())
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
  std::sort(current_data_.begin(), current_data_.end());
  ++iterations_done_;
}

}  // namespace cs
