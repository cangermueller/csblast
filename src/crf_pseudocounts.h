// Copyright 2009, Andreas Biegert

#ifndef CS_CRF_PSEUDOCOUNTS_H_
#define CS_CRF_PSEUDOCOUNTS_H_

#include "count_profile-inl.h"
#include "emission.h"
#include "profile-inl.h"
#include "pseudocounts.h"
#include "sequence-inl.h"
#include "crf-inl.h"

namespace cs {

// Encapsulation of context-specific pseudocounts calculated from a context CRF
template<class Abc>
class CrfPseudocounts : public Pseudocounts<Abc> {
 public:
  CrfPseudocounts(const Crf<Abc>& lib);

  virtual ~CrfPseudocounts() {}

  virtual void AddToSequence(const Sequence<Abc>& seq,
                             const Admix& pca,
                             Profile<Abc>& p) const;

  virtual void AddToProfile(const CountProfile<Abc>& cp,
                            const Admix& pca,
                            Profile<Abc>& p) const;

 private:
  // CRF with context weights and pseudocount emission weights.
  const Crf<Abc>& crf_;

  DISALLOW_COPY_AND_ASSIGN(CrfPseudocounts);
};  // CrfPseudocounts

}  // namespace cs

#endif  // CS_CRF_PSEUDOCOUNTS_H_
