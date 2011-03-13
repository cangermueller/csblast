// Copyright 2009, Andreas Biegert

#ifndef CS_CRF_H_
#define CS_CRF_H_

#include "count_profile.h"
#include "context_profile-inl.h"
#include "context_library-inl.h"
#include "crf_state.h"
#include "pseudocounts.h"
#include "sequence.h"

namespace cs {

// Forward declarations
template<class Abc>
class Crf;

// Strategy class for initializing a CRF
template<class Abc>
class CrfInit {
  public:
    CrfInit() {}
    virtual ~CrfInit() {}
    virtual void operator() (Crf<Abc>& crf) const = 0;
};


// A container class for CRF states to represent the most common
// sequence motifs in a training database of proteins/DNA sequences.
template<class Abc>
class Crf {
  public:
    typedef CrfState<Abc>* StateIter;
    typedef const CrfState<Abc>* ConstStateIter;

    // Constructs an empty CRF of given dimenions.
    Crf(size_t size, size_t wlen);

    // Constructs a CRF from serialized data read from input stream.
    explicit Crf(FILE* fin);

    // Constructs CRF with a specific init-strategy encapsulated by an
    // initializer.
    Crf(size_t size, size_t wlen, const CrfInit<Abc>& init);

    // Deallocates states in profile vector
    virtual ~Crf() {}

    // Returns the number of profiles in the fully assembled CRF
    size_t size() const { return states_.size(); }

    // Returns the number of columns in each context profile.
    size_t wlen() const { return wlen_; }

    // Returns total number of weights in this CRf. Note that context weights
    // and pseudocount weights of letter ANY are not accounted for since these are
    // held fix at zero anyway.
    size_t nweights() const { return size() * (1 + (wlen() + 1) * Abc::kSize); }

    // Returns index of central profile column.
    size_t center() const { return (wlen_ - 1) / 2; }

    // Accessor methods for state i, where i is from interval [0,size].
    CrfState<Abc>& operator[](size_t i) { return states_[i]; }
    const CrfState<Abc>& operator[](size_t i) const { return states_[i]; }

    // Initializes profile at index 'idx' with given profile.
    void SetState(size_t idx, const CrfState<Abc>& s);

    // Returns an iterator to a list of pointers to profiles.
    StateIter begin() { return &states_[0]; }

    // Returns an iterator pointing past the end of pointers to profiles.
    StateIter end() { return &states_[0] + states_.size(); }

    // Returns a const iterator over pointers of profiles.
    ConstStateIter begin() const { return &states_[0]; }

    // Returns a const iterator pointing past the end of pointers to profiles.
    ConstStateIter end() const { return &states_[0] + states_.size(); }

    // Writes the CRF in serialization format to output stream.
    void Write(FILE* fout) const;

  private:
    // Initializes the library from serialized data read from stream.
    void Read(FILE* fin);

    size_t wlen_;                           // size of context window.
    Vector< CrfState<Abc> > states_;  // states ordered by index.
};  // Crf


// Prints the library in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const Crf<Abc>& crf) {
    out << "CRF" << std::endl;
    out << "size:\t" << crf.size() << std::endl;
    out << "wlen:\t" << crf.wlen() << std::endl;
    for (size_t k = 0; k < crf.size(); ++k) out << crf[k];
    return out;
}


// Strategy for initializing CRF by sampling from training set, optionally
// adding pseudocounts.
template<class Abc, class TrainingPair>
class SamplingCrfInit : public CrfInit<Abc> {
  public:
    typedef std::vector<TrainingPair> TrainingSet;

    SamplingCrfInit(const TrainingSet& trainset,
                    const Pseudocounts<Abc>& pc,
                    const Admix& admix,
                    const SubstitutionMatrix<Abc>& sm,
                    unsigned int seed = 0)
            : trainset_(trainset),
              pc_(pc),
              admix_(admix),
              sm_(sm),
              seed_(seed) {}

    virtual ~SamplingCrfInit() {}

    virtual void operator() (Crf<Abc>& crf) const;

  private:
    const TrainingSet& trainset_;
    const Pseudocounts<Abc>& pc_;
    const Admix& admix_;
    const SubstitutionMatrix<Abc>& sm_;
    const unsigned int seed_;
};  // SamplingCrfInit


// FIXME: why doesn't this compile?!
// Comparison function for sorting context profiles by prior in descending order
// template<class Abc>
// struct ContextProfilePriorCompare :
//       public std::binary_function<ContextProfile<Abc>, ContextProfile<Abc>, bool> {
//   bool operator() (const ContextProfile<Abc>& lh,
//                    const ContextProfile<Abc>& rh) const {
//     return rh.prior < lh.prior;
//   }
// };

// template<class Abc>
// bool ContextProfilePriorCompare(const ContextProfile<Abc>& lhs,
//                                 const ContextProfile<Abc>& rhs) {
//   return rhs.prior < lhs.prior;
// }

// Strategy that uses context profiles from a profile library to initialize
// CRF states.
template<class Abc>
class LibraryBasedCrfInit : public CrfInit<Abc> {
  public:
    LibraryBasedCrfInit(const ContextLibrary<Abc>& lib,
                        const SubstitutionMatrix<Abc>& sm)
            : profiles_(lib.begin(), lib.end()), sm_(sm) {
        // std::sort(profiles_.begin(), profiles_.end(),
        //           ContextProfilePriorCompare<Abc>());
    }

    virtual ~LibraryBasedCrfInit() {}

    virtual void operator() (Crf<Abc>& crf) const {
        if (profiles_.size() < crf.size())
            throw Exception("Too few profiles in context library for CRF initialization!");
        for (size_t k = 0; k < crf.size(); ++k)
            crf.SetState(k, CrfState<Abc>(profiles_[k], sm_));
    }

  private:
    std::vector<ContextProfile<Abc> > profiles_;
    const SubstitutionMatrix<Abc>& sm_;
};  // class LibraryBasedCrfInit


// Strategy that initializes CRF weights by sammpling from gaussian distribution
template<class Abc>
class GaussianCrfInit : public CrfInit<Abc> {
  public:
    GaussianCrfInit(double sigma,
                    const SubstitutionMatrix<Abc>& sm,
                    unsigned int seed = 0)
            : sigma_(sigma), sm_(sm), seed_(seed) {}

    virtual ~GaussianCrfInit() {}

    virtual void operator() (Crf<Abc>& crf) const;

  protected:
    double sigma_;
    const SubstitutionMatrix<Abc>& sm_;
    unsigned int seed_;
};  // class GaussianCrfInit

}  // namespace cs

#endif  // CS_CRF_H_
