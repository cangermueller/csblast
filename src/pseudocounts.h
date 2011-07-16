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

// An abstract base class for pseudocount factories.
template<class Abc>
class Pseudocounts {
  public:
    Pseudocounts() {}
    virtual ~Pseudocounts() {}

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, const Admix& admix) const;

    // Adds pseudocounts to sequence using target Neff and returns normalized profile.
    Profile<Abc> AddTo(const Sequence<Abc>& seq, double neff, 
        double delta = kNeffDelta, double &tau = NeffTauDump) const;

    // Adds pseudocounts to sequence using admixture and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, const Admix& admix) const;

    // Adds pseudocounts to sequence using target Neff and returns normalized profile.
    Profile<Abc> AddTo(const CountProfile<Abc>& cp, double neff, 
        double delta = kNeffDelta, double &tau = NeffTauDump) const;

    // Adds pseudocounts to counts in PO-HMM vertices using admixture and stores results in 'probs' vector.
    void AddTo(POHmm<Abc>* hmm, const Admix& admix) const;

    // Adds pseudocounts to counts in PO-HMM vertices using target Neff and stores results in 'probs' vector.
    void AddTo(POHmm<Abc>* hmm, double neff, double delta = kNeffDelta, double &tau = NeffTauDump) const;

  private:
    // Adds pseudocounts to sequence and stores resulting frequencies in given
    // profile.
    virtual void AddToSequence(const Sequence<Abc>& seq, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToProfile(const CountProfile<Abc>& cp, Profile<Abc>& p) const = 0;

    // Adds pseudocounts to alignment derived profile.
    virtual void AddToPOHmm(const POHmm<Abc>* hhm, Profile<Abc>& p) const {};

    // Mixes profile 'p' and sequence 'q': tau * p + (1 - tau) * q
    void Mix(Profile<Abc>& p, const Sequence<Abc>& q, double tau) const;

    // Mixes profile 'p' and profile 'q': tau * p + (1 - tau) * q
    void Mix(Profile<Abc>& p, const Profile<Abc>& q, double tau) const;

    // Adjusts the Neff in 'p' to 'neff_' by admixing q and returns tau
    template<class T>
    double AdjustNeff(Profile<Abc>& p, const T& q, double neff, double delta) const;


  private:
    static const double kNormalize   = 1e-5; // Normalization threshold
    static const double kNeffTauInit = 0.8;  // Initial pseudocounts admixture for adjusting the Neff
    static const double kNeffTauMin  = 0.0;  // Minimal pseudocounts admixture for adjusting the Neff
    static const double kNeffTauMax  = 1.0;  // Maximal pseudocounts admixture for adjusting the Neff
    static double NeffTauDump;

  public:
    static const double kNeffDelta   = 0.1;  // Tolerance for adjusting the Neff

    DISALLOW_COPY_AND_ASSIGN(Pseudocounts);
};  // Pseudocounts

}  // namespace cs

#endif  // CS_PSEUDOCOUNTS_H_
