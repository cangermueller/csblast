// Copyright 2009, Andreas Biegert

#ifndef SRC_SUBSTITUTION_MATRIX_H_
#define SRC_SUBSTITUTION_MATRIX_H_

#include <cmath>
#include <cstdio>

#include <iostream>
#include <vector>

#include "matrix.h"

namespace cs {

// Abstract base class for substitution matrix classes.
template<class Alphabet>
class SubstitutionMatrix {
 public:
  SubstitutionMatrix();
  virtual ~SubstitutionMatrix() = 0;

  // Return substitution score S(a,b).
  float s(int a, int b) const { return s_[a][b]; }
  // Returns joint probability P(a,b).
  float p(int a, int b) const { return p_[a][b]; }
  // Returns probability P(a|b) "a given b".
  float r(int a, int b) const { return r_[a][b]; }
  // Returns background frequency of a.
  float f(int a) const { return f_[a]; }
  // Returns the size of the substitution matrix.
  int size() const { return size_; }

  // Prints the substitution matrix in human readable format to stream.
  friend std::ostream& operator<< (std::ostream& out,
                                   const SubstitutionMatrix& m) {
    m.Print(out);
    return out;
  }

 protected:
  // Initializes the other matrix data members from target frequencies.
  void init_from_target_frequencies();
  // Initializes the other matrix data members from score matrix and background
  // frequencies.
  void init_from_scores_and_background_frequencies();

  // Size of the matrix.
  const int size_;
  // Target frequency matrix P(a,b).
  matrix<float> p_;
  // Substitution matrix S(a,b).
  matrix<float> s_;
  // "a given b" probability matrix.
  matrix<float> r_;
  // Background frequencies of alphabet.
  std::vector<float> f_;

 private:
  // Prints the substitution matrix in human-readable format to output stream.
  virtual void Print(std::ostream& out) const;
};

}  // namespace cs

#endif  // SRC_SUBSTITUTION_MATRIX_H_
