// Copyright 2009, Andreas Biegert

#ifndef SRC_PSEUDOCOUNTS_H_
#define SRC_PSEUDOCOUNTS_H_

#include <algorithm>

#include "globals.h"

namespace cs
{

// Forward declarations
template<class Alphabet>
class Sequence;
template<class Alphabet>
class Profile;
template<class Alphabet>
class CountProfile;

// Calculates pseudocount admixture for profile column.
class Admixture {
 public:
  Admixture() {}
  virtual ~Admixture() {};

  virtual float operator() (float neff) const = 0;
};

// An abstract base class for pseudocount factories.
template<class Alphabet>
class Pseudocounts {
 public:
  Pseudocounts() {}
  virtual ~Pseudocounts() {}

  // Adds pseudocounts to sequence and stores resulting frequencies in given
  // profile.
  virtual void add_to_sequence(const Sequence<Alphabet>& seq,
                               const Admixture& pca,
                               Profile<Alphabet>* profile) const = 0;
  // Adds pseudocounts to alignment derived profile.
  virtual void add_to_profile(const Admixture& pca,
                              CountProfile<Alphabet>* profile) const = 0;

 private:
  DISALLOW_COPY_AND_ASSIGN(Pseudocounts);
};  // Pseudocounts



// Calculates constant pseudocount admixture independent of number of effective
// sequences.
class ConstantAdmixture : public Admixture {
 public:
  ConstantAdmixture(float x) : x_(x) {}
  ~ConstantAdmixture() {};

  virtual float operator() (float) const { return x_; }

 private:
  const float x_;
};

// Calculates divergence-dependent pseudocount admixture as
// tau = A * (B + 1) / (B + Neff)
class DivergenceDependentAdmixture : public Admixture {
 public:
  DivergenceDependentAdmixture(float a, float b) : a_(a), b_(b) {}
  ~DivergenceDependentAdmixture() {};

  virtual float operator() (float neff) const {
    return std::min(1.0f, a_ * (b_ + 1.0f) / (b_ + neff));
  }

 private:
  const float a_;
  const float b_;
};

}  // namespace cs

#endif  // SRC_PSEUDOCOUNTS_H_
