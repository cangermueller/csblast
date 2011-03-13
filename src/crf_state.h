// Copyright 2009, Andreas Biegert

#ifndef CS_CRF_STATE_H_
#define CS_CRF_STATE_H_

#include "profile_column.h"
#include "profile-inl.h"
#include "context_profile-inl.h"
#include "substitution_matrix.h"

namespace cs {

template<class Abc>
struct CrfState {
    // Default construction
    CrfState() : bias_weight(0.0) {};

    // Constructs a CRF state with 'len' columns
    explicit CrfState(size_t len)
            : bias_weight(0.0), context_weights(len) {
        assert(len & 1);
    }

    // Construction from serialized profile read from input stream.
    explicit CrfState(FILE* fin) { Read(fin); }

    // Constructs CRF with context weights initialized based on values in profile 'p'
    // and pseudocount weights initialized based on counts in column profile. We
    // assume that 'prof' and 'col' are both in lin-space.
    CrfState(const Profile<Abc>& prof,
             const ProfileColumn<Abc>& col,
             const SubstitutionMatrix<Abc>& sm)
            : bias_weight(0.0), context_weights(prof.length()) {
        assert(prof.length() & 1);

        Profile<Abc> p(prof);  // local copy for normalization
        Normalize(p, 1.0);
        ProfileColumn<Abc> c(col);  // local copy for normalization
        Normalize(c, 1.0);

        for (size_t j = 0; j < p.length(); ++j) {
            for (size_t a = 0; a < Abc::kSize; ++a)
                context_weights[j][a] = log(p[j][a] / sm.p(a));
            context_weights[j][Abc::kAny] = 0.0;
        }

        for (size_t a = 0; a < Abc::kSize; ++a)
            pc_weights[a] = log(c[a]);
        pc_weights[Abc::kAny] = 0.0;

        UpdatePseudocounts(*this);
    }

    // Constructs a CRF from probabilities in profile 'p' and background freqs.
    CrfState(const ContextProfile<Abc>& p, const SubstitutionMatrix<Abc>& sm)
            : bias_weight(p.is_log ? p.prior : log(p.prior)),
              context_weights(p.probs.length()) {
        assert(p.probs.length() & 1);

        for (size_t j = 0; j < p.probs.length(); ++j) {
            for (size_t a = 0; a < Abc::kSize; ++a) {
                double prob = p.is_log ? exp(p.probs[j][a]) : p.probs[j][a];
                context_weights[j][a] = log(prob / sm.p(a));
            }
            context_weights[j][Abc::kAny] = 0.0;
        }

        const int c = (p.probs.length() - 1) / 2;
        for (size_t a = 0; a < Abc::kSize; ++a) {
            double prob = p.is_log ? exp(p.probs[c][a]) : p.probs[c][a];
            pc_weights[a] = log(prob);
        }
        pc_weights[Abc::kAny] = 0.0;

        UpdatePseudocounts(*this);
    }

    // Initializes count profile with a serialized profile read from stream.
    void Read(FILE* fin);

    // Initializes count profile with a serialized profile read from stream.
    void Write(FILE* fin) const;

    std::string name;               // name of this state
    double bias_weight;             // bias weight lamda_k of this state
    Profile<Abc> context_weights;   // context weights lamda_k(j,a)
    ProfileColumn<Abc> pc_weights;  // unnormalized logs of pseudocounts
    ProfileColumn<Abc> pc;          // predicted pseudocounts at central column
};

// Prints CRF state weights in human-readable format for debugging.
template<class Abc>
std::ostream& operator<< (std::ostream& out, const CrfState<Abc>& crf);

// Updates pseudocount emission probs in given CRF state based on pc_weights.
template<class Abc>
void UpdatePseudocounts(CrfState<Abc>& state);

// Calculates context score between a CRF state and a sequence window
template<class Abc>
double ContextScore(const Profile<Abc>& context_weights,
                    const Sequence<Abc>& seq,
                    size_t idx,
                    size_t center);

// Calculates context score between a CRF state and a count profile window
template<class Abc>
double ContextScore(const Profile<Abc>& context_weights,
                    const CountProfile<Abc>& cp,
                    size_t idx,
                    size_t center);

}  // namespace cs

#endif  // CS_CRF_STATE_H_
