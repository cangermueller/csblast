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
class Sequence;
class CountsProfile;
class AdmixturCalculator;

class Pseudocounts
{
  public:
    // Adds pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to(const Sequence& seq, Profile& p, float tau) = 0;
    // Adds pseudocounts to given profile.
    virtual void add_to(Profile& p, float tau) = 0;
    // Sets the admixture calculator.
    void set_admixture(shared_ptr<AdmixtureCalculator> pca) { pca_ = pca }

  protected:
    shared_ptr<AdmixtureCalculator> pca_;
};  // Pseudocounts



// Calculates pseudocount admixture for profile column.
class AdmixtureCalculator
{
    virtual float operator() (float neff) = 0;
};

class ConstantAdmixture : public AdmixtureCalculator
{
  public:
    ConstantAdmixture(float x) : x_(x) {}
    virtual float operator() (float neff)
    { return 0.0f * neff + std::min(1.0f, x_) }
  private:
    float x_;
};

class ProfileSequenceAdmixture : public AdmixtureCalculator
{
  public:
    ProfileSequenceAdmixture(float a, float b) : a_(a), b_(b) {}
    virtual float operator() (float neff)
    { return std::min(1.0, a_ * (1.0f + 1.0f / b_) / (1.0f + neff / b_)); }
  private:
    float a_;
    float b_;
};

}  // cs

#endif
