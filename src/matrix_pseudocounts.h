#ifndef CS_MATRIX_PSEUDOCOUNTS_H
#define CS_MATRIX_PSEUDOCOUNTS_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// An abstract base class for pseudocount methods.

#include "pseudocounts.h"
#include "substitution_matrix.h"

namespace cs
{

// Forward declarations
class Sequence;
class Profile;
class CountsProfile;
class AdmixturCalculator;

class MatrixPseudocounts : public Pseudocounts
{
  public:
    MatrixPseudocounts(const SubstitutionMatrix& m);
    ~MatrixPseudocounts() {}

    // Adds substitution matrix pseudocounts to sequence and stores resulting frequencies in given profile.
    virtual void add_to_sequence(const Sequence& seq, const AdmixtureCalculator& pca, Profile& p);
    // Adds substitution matrix pseudocounts to alignment derived profile.
    virtual void add_to_profile(CountsProfile& p, const AdmixtureCalculator& pca);

  private:
    // Disallow copy and assign
    MatrixPseudocounts(const MatrixPseudocounts&);
    void operator=(const MatrixPseudocounts&);

    // Substitution matrix with conditional probabilities for pseudocounts.
    const SubstitutionMatrix& m_;
};  // MatrixPseudocounts

}  // cs

#endif
