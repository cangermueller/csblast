// Copyright 2009, Andreas Biegert

#ifndef CS_LR_PSEUDOCOUNTS_H_
#define CS_LR_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "lr_func-inl.h"
#include "pseudocounts-inl.h"
#include "blosum_matrix.h"
#include "tamura_nei_matrix.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from lr params
template<class Abc>
class LrPseudocounts : public Pseudocounts<Abc> {
 public:
  LrPseudocounts(const LrParams<Abc>& params) : params_(params) {}

  virtual ~LrPseudocounts() {}

  virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const;

  virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const;

 private:
  Sequence<Abc> ConsensusSeq(const CountProfile<Abc>& cp) const;
  
 private:
  const LrParams<Abc>& params_;

  DISALLOW_COPY_AND_ASSIGN(LrPseudocounts);
};  // LrPseudocounts

}  // namespace cs

#endif  // CS_LR_PSEUDOCOUNTS_H_
