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

#ifndef CS_SPARSE_PROFILE_H_
#define CS_SPARSE_PROFILE_H_

namespace cs {

// Implementation of a sparse-column specifically tailored to hold abstract-state
// probabilities. The maximum number of elements per column is limited to 255 in
// order to keep the mapping table reasonably small.
class SparseProfileCol {
  public:
    // Simple POD for storing index and probability sparsely
    struct Element {
        Element() : index(0), prob(0.0) {}
        Element(unsigned int i, double p) : index(i), prob(p) {}

        operator double() const { return prob; }
        bool operator< (const SparseProfileCol::Element& other) const {
            return prob > other.prob;
        }

        unsigned int index;
        double prob;
    };

    // Element iterators are read-only!
    typedef std::vector<Element>::const_iterator ElementIter;

    // Constructs an empty sparse column
    SparseProfileCol() {}

    // Constructs a sparse column from 'n' probabilities in array 'pp'. Note that
    // only those probabilities above the cutoff are retained.
    SparseProfileCol(const double* pp, size_t n, double cutoff = 0.0);

    // Fast accessor to probability at index 'k'. Iff element 'k' is assigned we
    // return its probability, otherwise we return zero.
    double operator[] (size_t k) const {
        return (map_[k] == kUnassigned) ? 0.0 : elements_[map_[k]];
    }

    // Returns the the maximum column index 'k' plus one.
    size_t size() const { return map_.size(); }

    // Returns true iff no elements are nonempty.
    bool empty() const { return elements_.empty(); }

    // Returns the number of assigned elements
    size_t num_nonempty() const { return elements_.size(); }

    // Returns a const iterator over assigned elements.
    ElementIter begin() const { return elements_.begin(); }

    // Returns a const iterator past the end of all assigned elements.
    ElementIter end() const { return elements_.end(); }

  private:
    // Needed to mark unassigned elements in map.
    static const uint8_t kUnassigned;
    // Maximal number of assigned elements
    static const size_t kMaxElements;

    std::vector<Element> elements_;  // assigned elements in sparse column
    Vector<uint8_t> map_;            // maps an index to its offset in 'elements'
};


// Prints an abstract state column in human-readable format for debugging.
std::ostream& operator<< (std::ostream& out, const SparseProfileCol& asc);

// Prints a column element in human-readable format for debugging.
std::ostream& operator<< (std::ostream& out,
                          const SparseProfileCol::Element& elem);

// // A class representing a sparse profile over abstract state probabilities.
// class AbstractStateProfile {
//  public:
//   // Simple struct for profile values.
//   struct Element {
//     Element() : index(0), prob(0.0f) {}
//     Element(int i, float p) : index(i), prob(p) {}
//     ~Element() {}

//     operator float() const { return prob; }

//     const int index;
//     float prob;
//   };  // class Element

//   typedef google::sparsetable<Element> ColType;
//   //typedef const google::sparsetable<const Element> ConstColType;
//   typedef ColType::const_nonempty_iterator ColIter;
//   //typedef ConstColType::const_nonempty_iterator ConstColIter;

//   // Constructs a empty profile with num_cols columns.
//   AbstractStateProfile(int num_cols, int num_states);
//   // Constructs a profile from a posterior probability matrix.
//   AbstractStateProfile(const matrix<double>& m, double prob_cutoff = 0.0);
//   // Deallocates columns
//   ~AbstractStateProfile();

//   // Access methods to get the (i,j) element
//   ColType operator[](int i) { return *profile_[i]; }
//   ColType operator[](int i) const { return *profile_[i]; }
//   // Returns the number columns in the profile
//   int num_cols() const { return num_cols_; }
//   int length() const { return num_cols_; }
//   // Returns maximal number of elements per column
//   int num_states() const { return num_states_; }
//   // Number of assigned elements in column i.
//   int num_elements(int i) const { return profile_[i]->num_nonempty(); }
//   // Returns a non-empty iterator to the first element in column i.
//   ColIter col_begin(int i) { return profile_[i]->nonempty_begin(); }
//   // Returns a non-empty iterator just past the end of column i.
//   ColIter col_end(int i) { return profile_[i]->nonempty_end(); }
//   // Returns a const non-empty iterator to the first element in column i.
//   ColIter col_begin(int i) const { return profile_[i]->nonempty_begin(); }
//   // Returns a const non-empty iterator just past the end of column i.
//   ColIter col_end(int i) const { return profile_[i]->nonempty_end(); }
//   // Sets matrix element at position (i,k) to given probability
//   void set(int i, int k, float p) { profile_[i]->set(k, Element(k,p)); }
//   // Clears all elements from the profile.
//   void clear();
//   // Writes the profile in human-readable format to output stream.
//   void Write(FILE* fp) const;

//   // Prints profile in human-readable format for debugging.
//   friend std::ostream& operator<< (std::ostream& out,
//                                    const AbstractStateProfile& p) {
//     p.Print(out);
//     return out;
//   }

//  private:
//   // Prints the profile in human-readable format to output stream.
//   void Print(std::ostream& out) const;
//   // Sets up profile data structure ready to be filled.
//   void Init();
//   // Initializes profile from posteriors in matrix.
//   void InitWithPosteriors(const matrix<double>& m, double prob_cutoff);

//   // Number of columns
//   const int num_cols_;
//   // Number of states per column
//   const int num_states_;
//   // Rows of profile matrix, each pointing to columns of elements.
//   ColType** profile_;
// };  // AbstractStateProfile

}  // namespace cs

#endif  // CS_SPARSE_PROFILE_H_
