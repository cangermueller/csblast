// Copyright 2009, Andreas Biegert

#ifndef CS_MATRIX_PSEUDOCOUNTS_INL_H_
#define CS_MATRIX_PSEUDOCOUNTS_INL_H_

#include "matrix_pseudocounts.h"

#include "count_profile-inl.h"
#include "po_hmm-inl.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

template<class Abc>
void MatrixPseudocounts<Abc>::AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const {
    assert_eq(seq.length(), p.length());
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to sequence ...";

    for(size_t i = 0; i < p.length(); ++i) {
        for(size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = m_.r(a, seq[i]);
    }
}

template<class Abc>
void MatrixPseudocounts<Abc>::AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const {
    assert_eq(cp.counts.length(), p.length());
    LOG(DEBUG) << "Adding substitution matrix pseudocounts to profile ...";

    for(size_t i = 0; i < cp.counts.length(); ++i) {
        for(size_t a = 0; a < Abc::kSize; ++a) {
            double sum = 0.0;
            for(size_t b = 0; b < Abc::kSize; ++b)
                sum += m_.r(a,b) * cp.counts[i][b] / cp.neff[i];
            p[i][a] = sum;
        }
    }
}

template<class Abc>
void MatrixPseudocounts<Abc>::AddToPOHmm(const POHmm<Abc>* hmm, Profile<Abc>& p) const {
    typedef typename POHmm<Abc>::Graph Graph;
    const Graph& g = hmm->g;

    for (size_t i = 1; i <= hmm->size(); ++i) {
        for(size_t a = 0; a < Abc::kSize; ++a) {
            double sum = 0.0;
            for(size_t b = 0; b < Abc::kSize; ++b)
                sum += m_.r(a,b) * g[i].counts[b] / g[i].neff;
            p[i - 1][a] = sum;
        }
    }
}

}  // namespace cs

#endif  // CS_MATRIX_PSEUDOCOUNTS_INL_H_
