// Copyright 2009, Andreas Biegert

#ifndef CS_MATRIX_PSEUDOCOUNTS_H_
#define CS_MATRIX_PSEUDOCOUNTS_H_

#include "count_profile.h"
#include "profile-inl.h"
#include "pseudocounts-inl.h"
#include "sequence.h"
#include "substitution_matrix.h"

namespace cs {

// Substitution matrix pseudocounts factory.
template<class Abc>
class MatrixPseudocounts : public Pseudocounts<Abc> {
  public:
    MatrixPseudocounts(const SubstitutionMatrix<Abc>& m) : m_(m) {};

    virtual ~MatrixPseudocounts() {}

  private:
    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const;

    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const;

    virtual void AddToPOHmm(const POHmm<Abc>* hmm, Profile<Abc>& p) const;

    // Substitution matrix with conditional probabilities for pseudocounts.
    const SubstitutionMatrix<Abc>& m_;

    DISALLOW_COPY_AND_ASSIGN(MatrixPseudocounts);
};  // MatrixPseudocounts

}  // namespace cs

#endif  // CS_MATRIX_PSEUDOCOUNTS_H_
