#ifndef CS_PSEUDOCOUNTS_H
#define CS_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// An abstract base class for pseudocount methods.

#include <algorithm>

#include "shared_ptr.h"

namespace cs
{

// Forward declarations
template<class Alphabet_T>
class Sequence;
template<class Alphabet_T>
class Profile;
template<class Alphabet_T>
class CountsProfile;

// Calculates pseudocount admixture for profile column.
class Admixture
{
  public:
    Admixture() {}
    virtual ~Admixture() {};

    virtual float calculate(float neff) const = 0;
};

template<class Alphabet_T>
class Pseudocounts
{
  public:
    Pseudocounts(shared_ptr<Admixture> pca) : pca_(pca) {}
    virtual ~Pseudocounts() {}

    // Adds pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence<Alphabet_T>& seq,
                                 Profile<Alphabet_T>& p) = 0;
    // Adds pseudocounts to alignment derived profile.
    virtual void add_to_profile(CountsProfile<Alphabet_T>& p) = 0;

    // Sets pseudocount admixture formula to be used.
    void set_admixture(shared_ptr<Admixture> pca) { pca_ = pca; }

  protected:
    // Calculates pseudocount admixute depending on alignment diversity.
    float admixture(float neff) const;

  private:
    // Disallow copy and assign
    Pseudocounts(const Pseudocounts&);
    void operator=(const Pseudocounts&);

    // The admixture formula to be used.
    shared_ptr<Admixture> pca_;
};  // Pseudocounts



template<class Alphabet_T>
float Pseudocounts<Alphabet_T>::admixture(float neff) const
{
    return pca_->calculate(neff);
}


// Calculates constant pseudocount admixture independent of number of effective sequences.
class ConstantAdmixture : public Admixture
{
  public:
    ConstantAdmixture(float x) : x_(x) {}
    ~ConstantAdmixture() {};

    virtual float calculate(float neff) const { return 0 * neff + std::min(1.0f, x_); }

  private:
    const float x_;
};

// Calculates divergence-dependent pseudocount admixture as tau = A * (B + 1) / (B + Neff)
class DivergenceDependentAdmixture : public Admixture
{
  public:
    DivergenceDependentAdmixture(float a, float b) : a_(a), b_(b) {}
    ~DivergenceDependentAdmixture() {};

    virtual float calculate(float neff) const
    { return std::min(1.0f, a_ * (b_ + 1.0f) / (b_ + neff)); }

  private:
    const float a_;
    const float b_;
};

}  // cs

#endif
