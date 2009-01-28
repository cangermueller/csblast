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

class ConstantAdmixture : public AdmixtureCalculator
{
  public:
    ConstantAdmixture(float x) : x_(x) {}
    virtual float operator() (float neff) const
    { return 0.0f * neff + std::min(1.0f, x_); }
  private:
    float x_;
};

class ProfileSequenceAdmixture : public AdmixtureCalculator
{
  public:
    ProfileSequenceAdmixture(float a, float b) : a_(a), b_(b) {}
    virtual float operator() (float neff) const
    { return std::min(1.0f, a_ * (1.0f + 1.0f / b_) / (1.0f + neff / b_)); }
  private:
    float a_;
    float b_;
};

}  // cs

#endif
