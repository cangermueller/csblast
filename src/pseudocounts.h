// Copyright 2009, Andreas Biegert

#ifndef CS_PSEUDOCOUNTS_H_
#define CS_PSEUDOCOUNTS_H_

namespace cs {

// Forward declarations
template<class Abc>
class Sequence;
template<class Abc>
class Profile;
template<class Abc>
class CountProfile;
template<class Abc>
class POHmm;

// Calculates pseudocount admixture for profile column.
struct Admix {
  public:
    Admix() {}
    virtual ~Admix() {}

    virtual double operator() (double neff) const = 0;
};

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



// An abstract base class for pseudocount factories.
template<class Abc>
class Pseudocounts {
  public:
    Pseudocounts() {}
    virtual ~Pseudocounts() {}

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, const Admix& admix) const;

    // Adds pseudocounts to sequence using target Neff and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, CSBlastAdmix& admix,
        double neff, double delta = kAdjustDelta) const;

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, const Admix& admix) const;

    // Adds pseudocounts to sequence using target Neff and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, CSBlastAdmix& admix,
        double neff, double delta = kAdjustDelta) const;

    // Adds pseudocounts to counts in PO-HMM vertices using admixture and stores results in 'probs' vector.
    void AddTo(POHmm<Abc>* hmm, const Admix& admix) const;

    // Adds pseudocounts to counts in PO-HMM vertices using target Neff and stores results in 'probs' vector.
    void AddTo(POHmm<Abc>* hmm, CSBlastAdmix& admix,
        double neff, double delta = kAdjustDelta) const;

  private:
    // Adds pseudocounts to sequence and stores resulting frequencies in given
    // profile.
    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToPOHmm(const POHmm<Abc>* hhm, Profile<Abc>& p) const {};

    // Admixes Sequence q to Profile p
    void AdmixTo(const Sequence<Abc>& q, Profile<Abc>& p, const Admix& admix) const;

    // Admixes CountProfile q to Profile q
    void AdmixTo(const CountProfile<Abc>& q, Profile<Abc>& p, const Admix& admix) const;

    // Admixes POHmm q to profile q
    void AdmixTo(const POHmm<Abc>& q, Profile<Abc>& p, const Admix& admix) const;

    // Admixes q to Profile p such that the Neff in p is equal to neff
    template<class T>
    double AdmixToNeff(const T& q, Profile<Abc>& p, CSBlastAdmix& admix, 
        double neff, double delta = kAdjustDelta) const;

  private:
    static const double kNormalize  = 1e-5; // Normalization threshold
    static const double kAdjustMin  = 0.0;  // Minimal paramater value for adjusting the Neff
    static const double kAdjustMax  = 1.0;  // Maximal parameter value for adjusting the Neff
    static const double kAdjustInit = 0.5;  // Initial parameter value for adjusting the Neff
    static const double kAdjustEps  = 0.01; // Convergence threshold for adjusting the Neff

  public:
    static const double kAdjustDelta   = 0.025;  // Tolerance for adjusting the Neff

    DISALLOW_COPY_AND_ASSIGN(Pseudocounts);
};  // Pseudocounts

}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_H_
