// Copyright 2009, Andreas Biegert

#ifndef CS_PSEUDOCOUNTS_INL_H_
#define CS_PSEUDOCOUNTS_INL_H_
#include "pseudocounts.h"

namespace cs {


// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, const Admix& admix) const {
    Profile<Abc> p(seq.length());
    AddToSequence(seq, p);
    Mix(p, seq, admix);
    for(size_t i = 0; i < seq.length(); ++i) p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to sequence using target Neff and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, CSBlastAdmix& admix,
        double neff, double delta) const {
    Profile<Abc> p(seq.length());
    AddToSequence(seq, p);
    AdjustNeff(p, seq, admix, neff, delta);
    for(size_t i = 0; i < seq.length(); ++i) p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, const Admix& admix) const {
    Profile<Abc> p(cp.counts.length());
    AddToProfile(cp, p);
    Mix(p, cp, admix);
    for(size_t i = 0; i < cp.counts.length(); ++i)
        p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to sequence using target Neff and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, CSBlastAdmix& admix,
        double neff, double delta) const {
    Profile<Abc> p(cp.counts.length());
    AddToProfile(cp, p);
    AdjustNeff(p, cp, admix, neff, delta);
    for(size_t i = 0; i < cp.counts.length(); ++i)
        p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to counts in PO-HMM vertices using admixture and stores results in 'probs' vector.
template<class Abc>
void Pseudocounts<Abc>::AddTo(POHmm<Abc>* hmm, const Admix& admix) const {
    typename POHmm<Abc>::Graph& g = hmm->g;
    size_t size = hmm->size();

    Profile<Abc> p(size);
    AddToPOHmm(hmm, p);
    Mix(p, *hmm, admix);

    for (size_t i = 1; i <= size; ++i) {
        memcpy(&g[i].probs[0], p[i - 1], Abc::kSize * sizeof(double));
        g[i].probs[Abc::kAny] = 1.0;
        Normalize(g[i].probs, 1.0);
    }
}

// Adds pseudocounts to counts in PO-HMM vertices and stores results in 'probs' vector.
template<class Abc>
void Pseudocounts<Abc>::AddTo(POHmm<Abc>* hmm, CSBlastAdmix& admix,
        double neff, double delta) const {
    typename POHmm<Abc>::Graph& g = hmm->g;
    size_t size = hmm->size();

    Profile<Abc> p(size);
    AddToPOHmm(hmm, p);
    AdjustNeff(p, *hmm, admix, neff, delta);
    for (size_t i = 1; i <= size; ++i) {
        memcpy(&g[i].probs[0], p[i - 1], Abc::kSize * sizeof(double));
        g[i].probs[Abc::kAny] = 1.0;
        Normalize(hmm->g[i].probs, 1.0);
    }
}

template<class Abc>
// Mixes Profile p and Sequence q
void Pseudocounts<Abc>::Mix(Profile<Abc>& p, const Sequence<Abc>& q, const Admix& admix) const {
    double tau = admix(1.0);
    double t = 1 - tau;
    for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] *= tau;
        p[i][q[i]] += t;
    }
}

template<class Abc>
// Mixes Profile p and CountProfile q
void Pseudocounts<Abc>::Mix(Profile<Abc>& p, const CountProfile<Abc>& q, const Admix& admix) const {
    for (size_t i = 0; i < p.length(); ++i) {
        double tau = admix(q.neff[i]);
        double t = 1 - tau;
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = tau * p[i][a] + t * q.counts[i][a] / q.neff[i];
    }
}

template<class Abc>
// Mixes Profile p and POHmm q
void Pseudocounts<Abc>::Mix(Profile<Abc>& p, const POHmm<Abc>& q, const Admix& admix) const {
    const typename POHmm<Abc>::Graph& g = q.g;

    for (size_t i = 1; i <= q.size(); ++i) {
        double tau = admix(g[i].neff);
        double t = 1 - tau;
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i - 1][a] = tau * p[i - 1][a] + t * g[i].counts[a] / g[i].neff;
    }
}

// Adjusts the Neff in 'p' to 'neff' by admixing q and returns tau
template<class Abc>
template<class T>
double Pseudocounts<Abc>::AdjustNeff(Profile<Abc>& p, const T& q, CSBlastAdmix& admix, 
    double neff, double delta) const {

    double l  = kAdjustMin;
    double r  = kAdjustMax;
    admix.pca = kAdjustInit;
    Profile<Abc> pp;
    while (l < kAdjustMax - kAdjustEps && r > kAdjustMin + kAdjustEps) {
        pp = p;
        Mix(pp, q, admix);
        double ne = Neff(pp);
        if (fabs(ne - neff) <= delta) {
            break;
        } else {
            if (ne < neff) l = admix.pca;
            else r = admix.pca;
        }
        admix.pca = 0.5 * (l + r);
    }
    if (l > kAdjustMax - kAdjustEps) {
        admix.pca = kAdjustMax;
        Mix(p, q, admix);
    } else if (r < kAdjustMin + kAdjustEps) {
        admix.pca = kAdjustMin;
        Mix(p, q, admix);
    } else {
        p = pp;
    }
    return admix.pca;
}


}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_INL_H_
