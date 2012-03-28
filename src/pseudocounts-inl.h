// Copyright 2009, Andreas Biegert

#ifndef CS_PSEUDOCOUNTS_INL_H_
#define CS_PSEUDOCOUNTS_INL_H_
#include "pseudocounts.h"

namespace cs {


// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, Admix& admix) const {
    Profile<Abc> p(seq.length());
    AddToSequence(seq, p);
    if (target_neff_ >= 1.0) {
      AdmixToTargetNeff(seq, p, admix);
    } else {
      AdmixTo(seq, p, admix);
    }
    for(size_t i = 0; i < seq.length(); ++i) p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, Admix& admix) const {
    Profile<Abc> p(cp.counts.length());
    AddToProfile(cp, p);
    if (target_neff_ >= 1.0) {
      AdmixToTargetNeff(cp, p, admix);
    } else {
      AdmixTo(cp, p, admix);
    }
    for(size_t i = 0; i < cp.counts.length(); ++i)
        p[i][Abc::kAny] = 0.0;
    Normalize(p, 1.0);
    return p;
}

// Adds pseudocounts to counts in PO-HMM vertices using admixture and stores results in 'probs' vector.
template<class Abc>
void Pseudocounts<Abc>::AddTo(POHmm<Abc>* hmm, Admix& admix) const {
    typename POHmm<Abc>::Graph& g = hmm->g;
    size_t size = hmm->size();

    Profile<Abc> p(size);
    AddToPOHmm(hmm, p);
    if (target_neff_ >= 1.0) {
      AdmixToTargetNeff(*hmm, p, admix);
    } else {
      AdmixTo(*hmm, p, admix);
    }
    for (size_t i = 1; i <= size; ++i) {
        memcpy(&g[i].probs[0], p[i - 1], Abc::kSize * sizeof(double));
        g[i].probs[Abc::kAny] = 1.0;
        Normalize(g[i].probs, 1.0);
    }
}

template<class Abc>
void Pseudocounts<Abc>::AdmixTo(const Sequence<Abc>& q, Profile<Abc>& p, const Admix& admix) const {
    double tau = admix(1.0);
    double t = 1 - tau;
    for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] *= tau;
        p[i][q[i]] += t;
    }
}

template<class Abc>
void Pseudocounts<Abc>::AdmixTo(const CountProfile<Abc>& q, Profile<Abc>& p, const Admix& admix) const {
    for (size_t i = 0; i < p.length(); ++i) {
        double tau = admix(q.neff[i]);
        double t = 1 - tau;
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = tau * p[i][a] + t * q.counts[i][a] / q.neff[i];
    }
}

template<class Abc>
void Pseudocounts<Abc>::AdmixTo(const POHmm<Abc>& q, Profile<Abc>& p, const Admix& admix) const {
    const typename POHmm<Abc>::Graph& g = q.g;

    for (size_t i = 1; i <= q.size(); ++i) {
        double tau = admix(g[i].neff);
        double t = 1 - tau;
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i - 1][a] = tau * p[i - 1][a] + t * g[i].counts[a] / g[i].neff;
    }
}

// Adjusts the Neff in 'p' to 'neff' by admixing q and returns tau.
template<class Abc>
template<class T>
double Pseudocounts<Abc>::AdmixToTargetNeff(const T& q, Profile<Abc>& p, Admix& admix) const {

    double l = kTargetNeffParamMin;
    double r = kTargetNeffParamMax;
    admix.SetTargetNeffParam(kTargetNeffParamInit);
    Profile<Abc> pp;
    while (l < kTargetNeffParamMax - kTargetNeffEps && r > kTargetNeffParamMin + kTargetNeffEps) {
        pp = p;
        AdmixTo(q, pp, admix);
        double ne = Neff(pp);
        if (fabs(ne - target_neff_) <= target_neff_delta_) {
            break;
        } else {
            if (ne < target_neff_) l = admix.GetTargetNeffParam();
            else r = admix.GetTargetNeffParam();
        }
        admix.SetTargetNeffParam(0.5 * (l + r));
    }
    if (l > kTargetNeffParamMax - kTargetNeffEps) {
        admix.SetTargetNeffParam(kTargetNeffParamMax);
        AdmixTo(q, p, admix);
    } else if (r < kTargetNeffParamMin + kTargetNeffEps) {
        admix.SetTargetNeffParam(kTargetNeffParamMin);
        AdmixTo(q, p, admix);
    } else {
        p = pp;
    }
    return admix.GetTargetNeffParam();
}


}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_INL_H_
