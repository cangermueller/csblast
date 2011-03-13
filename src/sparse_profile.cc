// Copyright 2009, Andreas Biegert

#include "cs.h"
#include "sparse_profile.h"

namespace cs {

// Needed to mark unassigned elements in map.
const uint8_t SparseProfileCol::kUnassigned = 255;

// Maximal number of assigned elements
const size_t SparseProfileCol::kMaxElements = 255;

SparseProfileCol::SparseProfileCol(const double* pp, size_t n, double cutoff)
    : map_(n, kUnassigned) {
  // Add all posteriors above cutoff to elems table
  size_t k_max = 0;
  double p_max = 0.0;
  for (size_t k = 0; k < n; ++k) {
      if (pp[k] >= cutoff) {
        elements_.push_back(Element(k, pp[k]));
      } else if(pp[k] > p_max) {
        k_max = k;
        p_max = pp[k];
      }
  }
  // If no posterior was above the cutoff, add at least the larget one
  if (elements_.empty())
    elements_.push_back(Element(k_max, p_max));
  // Order elements by descending posterior probs
  std::sort(elements_.begin(), elements_.end());
  // Remove all elements above capacity
  if (elements_.size() > kMaxElements)
    elements_.erase(elements_.begin() + kMaxElements, elements_.end());
  // Shrink to fit to save space
  std::vector<Element>(elements_).swap(elements_);
  // Record offsets of assigned elements
  for (size_t i = 0; i < elements_.size(); ++i)
    map_[elements_[i].index] = i;
}

std::ostream& operator<< (std::ostream& out, const SparseProfileCol& asc) {
  out << StringifyRange(asc.begin(), asc.end());
  return out;
}

std::ostream& operator<< (std::ostream& out,
                          const SparseProfileCol::Element& elem) {
  out << strprintf("[index=%-4d prob=%4.2f]", elem.index, elem.prob);
  return out;
}

// AbstractStateProfile::AbstractStateProfile(int num_cols, int num_states)
//     : num_cols_(num_cols),
//       num_states_(num_states),
//       profile_(NULL) {
//   Init();
// }

// AbstractStateProfile::AbstractStateProfile(const matrix<double>& m,
//                                            double prob_cutoff)
//     : num_cols_(m.num_rows()),
//       num_states_(m.num_cols()),
//       profile_(NULL) {
//   Init();
//   InitWithPosteriors(m, prob_cutoff);
// }

// AbstractStateProfile::~AbstractStateProfile() {
//   for (int i = 0; i < num_cols_; ++i) {
//     delete profile_[i];
//     profile_[i] = NULL;
//   }
//   delete[] profile_;
//   profile_ = NULL;
// }

// void AbstractStateProfile::Init() {
//   profile_ = new ColType*[num_cols_];
//   for (int i = 0; i < num_cols_; ++i)
//     profile_[i] = new ColType(num_states_);
// }

// void AbstractStateProfile::InitWithPosteriors(const matrix<double>& m,
//                                               double prob_cutoff) {
//   for (int i = 0; i < num_cols_; ++i) {
//     int index_max = 0;
//     double prob_max = 0.0;

//     for (int k = 0; k < num_states_; ++k) {
//       if (m[i][k] >= prob_cutoff) {
//         set(i, k, m[i][k]);
//       } else if(m[i][k] > prob_max) {
//         index_max = k;
//         prob_max  = m[i][k];
//       }
//     }
//     if (num_elements(i) == 0)
//       set(i, index_max, prob_max);
//   }
// }


// void AbstractStateProfile::Print(std::ostream& out) const {
//   for (int i = 0; i < num_cols_; ++i) {
//     out << strprintf("% 4i:  ", i+1);
//     for (ColIter it = col_begin(i); it != col_end(i); ++it)
//       out << strprintf("   %4i -> %5.1f%%", it->index, it->prob * 100);
//     out << std::endl;
//   }
// }

// void AbstractStateProfile::Write(FILE* fp) const {
//   for (int i = 0; i < num_cols_; ++i) {
//     fprintf(fp, "% 4i:  ", i+1);
//     for (ColIter it = col_begin(i); it != col_end(i); ++it)
//       fprintf(fp, "   %4i -> %5.1f%%", it->index, it->prob * 100);
//     fputs("\n", fp);
//   }
// }

}  // namespace cs
