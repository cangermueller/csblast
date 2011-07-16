// Copyright 2009, Andreas Biegert

#ifndef CS_PSEUDOCOUNTS_INL_H_
#define CS_PSEUDOCOUNTS_INL_H_
#include "pseudocounts.h"

namespace cs {


// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, const Admix& admix) const {
    Profile<Abc> rv(seq.length());
    AddToSequence(seq, rv);
    Mix(rv, seq, admix(1.0));
    for(size_t i = 0; i < seq.length(); ++i) rv[i][Abc::kAny] = 0.0;
    Normalize(rv, 1.0);
    return rv;
}

// Adds pseudocounts to sequence using target Neff and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const Sequence<Abc>& seq, double neff, 
                                      double delta, double &tau) const {
    Profile<Abc> rv(seq.length());
    AddToSequence(seq, rv);
    tau = AdjustNeff(rv, seq, neff, delta);
    for(size_t i = 0; i < seq.length(); ++i) rv[i][Abc::kAny] = 0.0;
    Normalize(rv, 1.0);
    return rv;
}

// Adds pseudocounts to sequence using admixture and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, const Admix& admix) const {
    Profile<Abc> rv(cp.counts.length());
    AddToProfile(cp, rv);
    for (size_t i = 0; i < rv.length(); ++i) {
        double tau = admix(cp.neff[i]);
        for (size_t a = 0; a < Abc::kSize; ++a) {
            rv[i][a] = (1 - tau) * cp.counts[i][a] / cp.neff[i] +
                tau * rv[i][a];
        }
    }
    for(size_t i = 0; i < cp.counts.length(); ++i)
        rv[i][Abc::kAny] = 0.0;
    Normalize(rv, 1.0);
    return rv;
}

// Adds pseudocounts to sequence using target Neff and returns normalized profile.
template<class Abc>
Profile<Abc> Pseudocounts<Abc>::AddTo(const CountProfile<Abc>& cp, double neff, 
                                      double delta, double &tau) const {
    Profile<Abc> rv(cp.counts.length());
    AddToProfile(cp, rv);
    tau = AdjustNeff(rv, Profile<Abc>(cp), neff, delta);
    for(size_t i = 0; i < cp.counts.length(); ++i)
        rv[i][Abc::kAny] = 0.0;
    Normalize(rv, 1.0);
    return rv;
}

// Adds pseudocounts to counts in PO-HMM vertices using admixture and stores results in 'probs' vector.
template<class Abc>
void Pseudocounts<Abc>::AddTo(POHmm<Abc>* hmm, const Admix& admix) const {
    typedef typename POHmm<Abc>::Graph Graph;
    Graph& g = hmm->g;
    size_t size = hmm->size();

    Profile<Abc> rv(size);
    AddToPOHmm(hmm, rv);
    for (size_t i = 1; i <= size; ++i) {
        double tau = admix(g[i].neff);
        for (size_t a = 0; a < Abc::kSize; ++a) {
            g[i].probs[a] = (1 - tau) * g[i].counts[a] / g[i].neff +
                tau * rv[i - 1][a];
        }
        g[i].probs[Abc::kAny] = 1.0;
        Normalize(g[i].probs, 1.0);
    }
}

// Adds pseudocounts to counts in PO-HMM vertices and stores results in 'probs' vector.
template<class Abc>
void Pseudocounts<Abc>::AddTo(POHmm<Abc>* hmm, double neff, double delta, double &tau) const {
    typedef typename POHmm<Abc>::Graph Graph;
    Graph& g = hmm->g;
    size_t size = hmm->size();

    Profile<Abc> rv(size);
    AddToPOHmm(hmm, rv);
    Profile<Abc> p_hmm(size);
    for (size_t i = 1; i <= size; ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p_hmm[i - 1][a] = g[i].counts[a] / g[i].neff;
    }
    tau = AdjustNeff(rv, p_hmm, neff, delta);
    for (size_t i = 1; i <= size; ++i) {
        memcpy(&g[i].probs[0], rv[i - 1], Abc::kSize * sizeof(double));
        g[i].probs[Abc::kAny] = 1.0;
        Normalize(hmm->g[i].probs, 1.0);
    }
}

// Mixes profile 'p' and sequence 'q': tau * p + (1 - tau) * q
template<class Abc>
void Pseudocounts<Abc>::Mix(Profile<Abc>& p, const Sequence<Abc>& q, double tau) const {
    double t = 1.0 - tau;
    for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] *= tau;
        p[i][q[i]] += t;
    }
}

// Mixes profile 'p' and profile 'q': tau * p + (1 - tau) * q
template<class Abc>
void Pseudocounts<Abc>::Mix(Profile<Abc>& p, const Profile<Abc>& q, double tau) const {
    double t = 1.0 - tau;
    for (size_t i = 0; i < p.length(); ++i) {
        for (size_t a = 0; a < Abc::kSize; ++a)
            p[i][a] = tau * p[i][a] + t * q[i][a];
    }
}

// Adjusts the Neff in 'p' to 'neff' by admixing q and returns tau
template<class Abc>
template<class T>
double Pseudocounts<Abc>::AdjustNeff(Profile<Abc>& p, const T& q, double neff, double delta) const {
    double l         = kNeffTauMin;
    double r         = kNeffTauMax;
    double tau       = kNeffTauInit;
    const double EPS = 0.01;
    Profile<Abc> pp;
    while (l < kNeffTauMax - EPS && r > kNeffTauMin + EPS) {
        pp = p;
        Mix(pp, q, tau);
        double ne = Neff(pp);
        if (fabs(ne - neff) <= delta) {
            break;
        } else {
            if (ne < neff) l = tau;
            else r = tau;
        }
        tau = 0.5 * (l + r);
    }
    if (l > kNeffTauMax - EPS) {
        tau = 1.0; 
    } else if (r < kNeffTauMin + EPS) {
        tau = 0.0;
        Mix(p, q, tau);
    } else {
        p = pp;
    }
    return tau;
}



// Admixture methods



// Calculates constant pseudocount admixture independent of number of effective
// sequences.
struct ConstantAdmix : public Admix {
    ConstantAdmix(double a) : pca(a) {}
    virtual ~ConstantAdmix() {}

    virtual double operator() (double) const { return pca; }

    const double pca;
};

// Calculates divergence-dependent pseudocount admixture as in CS-BLAST
// tau = A * (B + 1) / (B + Neff)
struct CSBlastAdmix : public Admix {
    CSBlastAdmix(double a, double b) : pca(a), pcb(b) {}
    virtual ~CSBlastAdmix() {}

    virtual double operator() (double neff) const {
        return MIN(1.0, pca * (pcb + 1.0) / (pcb + neff));
    }

    double pca, pcb;
};

// Calculates divergence-dependent pseudocount admixture as in HHsearch
struct HHsearchAdmix : public Admix {
    HHsearchAdmix(double a, double b, double c = 1.0) : pca(a), pcb(b), pcc(c) {}
    virtual ~HHsearchAdmix() {}

    virtual double operator() (double neff) const {
        double rv = 0.0;
        if (pcc == 1.0)
            rv = MIN(1.0, pca / (1.0 + neff / pcb));
        else
            rv = MIN(1.0, pca / (1.0 + pow(neff / pcb, pcc)));
        return rv;
    }

    double pca, pcb, pcc;
};

template<class Abc>
double Pseudocounts<Abc>::NeffTauDump = 0.0;

}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_INL_H_
