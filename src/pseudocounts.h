#ifndef CS_PSEUDOCOUNTS_H
#define CS_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// An abstract base class for pseudocount methods.

#include <algorithm>

namespace cs
{

// Forward declarations
class Sequence;
class Profile;
class CountsProfile;
class AdmixtureCalculator;

class Pseudocounts
{
  public:
    Pseudocounts() {}
    ~Pseudocounts() {}

    // Adds pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence& seq, const AdmixtureCalculator& pca, Profile& p) = 0;
    // Adds pseudocounts to alignment derived profile.
    virtual void add_to_profile(CountsProfile& p, const AdmixtureCalculator& pca) = 0;

  private:
    // Disallow copy and assign
    Pseudocounts(const Pseudocounts&);
    void operator=(const Pseudocounts&);
};  // Pseudocounts



// Calculates pseudocount admixture for profile column.
class AdmixtureCalculator
{
  public:
    virtual float operator() (float neff) const = 0;
};

// Calculates constant pseudocount admixture independent of number of effective sequences.
class ConstantAdmixture : public AdmixtureCalculator
{
  public:
    ConstantAdmixture(float x) : x_(x) {}
    ~ConstantAdmixture() {};

    virtual float operator() (float neff) const { return 0.0f * neff + std::min(1.0f, x_); }

  private:
    const float x_;
};

// Calculates divergence-dependent pseudocount admixture as tau = A * (B + 1) / (B + Neff)
class DivergenceDependentAdmixture : public AdmixtureCalculator
{
  public:
    DivergenceDependentAdmixture(float a, float b) : a_(a), b_(b) {}
    ~DivergenceDependentAdmixture() {};

    virtual float operator() (float neff) const
    { return std::min(1.0f, a_ * (b_ + 1.0f) / (b_ + neff)); }

  private:
    const float a_;
    const float b_;
};

}  // cs

#endif
