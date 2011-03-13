// Copyright 2009, Andreas Biegert

#ifndef CS_FORWARD_BACKWARD_ALGORITHM_INL_H_
#define CS_FORWARD_BACKWARD_ALGORITHM_INL_H_

#include <valarray>

#include "context_profile-inl.h"
#include "mult_emission-inl.h"
#include "hmm-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils.h"

namespace cs {

// POD with forward backward matrices, emission probs, and log-likelihood.
struct ForwardBackwardMatrices {
  ForwardBackwardMatrices(int slen, int nstates)
      : subject_length(slen),
        num_states(nstates),
        f(slen, nstates, 0.0),
        b(slen, nstates, 0.0),
        pp(slen, nstates, 0.0),
        e(slen, nstates, 0.0),
        s(0.0, slen),
        logp(0.0) {}

  friend std::ostream& operator<< (std::ostream& out,
                                   const ForwardBackwardMatrices& m) {
    out << "Forward matrix f[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%7.5f  ", m.f[i][k]);
      }
      out << std::endl;
    }
    out << "Backward matrix b[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%7.5f  ", m.b[i][k]);
      }
      out << std::endl;
    }
    out << "Forward row-sums before scaling s[i]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      out << strprintf("%9.3g  ", m.s[i]);
    }
    out << "\nMultEmission probabilities e[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%8.3g  ", m.e[i][k]);
      }
      out << std::endl;
    }
    out << "Posterior probabilities pp[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        double p = 100.0 * m.pp[i][k];
        out << strprintf("%6.2f  ", p);
      }
      out << std::endl;
    }
    out << strprintf("Log-likelihood = %-10.5f\n", m.logp);
    return out;
  }

  // length of the subject whose Fwd-Bwd information is stored
  const int subject_length;
  // Number of states in the HMM against which the subject is aligned
  const int num_states;
  // Forward matrix
  matrix<double> f;
  // Backward matrix
  matrix<double> b;
  // Posterior probability matrix
  matrix<double> pp;
  // Emission probability matrix
  matrix<double> e;
  // Forward matrix column sum s(n) before normalization
  std::valarray<double> s;
  // Log-likelihood calculated with forward algorithm.
  double logp;
};


// Forward-Backward algorithm encapsulation.
template< class Alphabet, template<class> class Subject >
void ForwardBackwardAlgorithm(const Hmm<Alphabet>& hmm,
                              const Subject<Alphabet>& subject,
                              const MultEmission<Alphabet>& emission,
                              ForwardBackwardMatrices* fbm) {
  LOG(DEBUG) << "Running forward-backward algorithm ...";
  LOG(DEBUG1) << hmm;
  LOG(DEBUG1) << subject;

  ForwardAlgorithm(hmm, subject, emission, fbm);
  BackwardAlgorithm(hmm, subject, fbm);
}

template< class Alphabet, template<class> class Subject >
void ForwardAlgorithm(const Hmm<Alphabet>& hmm,
                      const Subject<Alphabet>& subject,
                      const MultEmission<Alphabet>& emission,
                      ForwardBackwardMatrices* fbm) {
  typedef ContextProfileState<Alphabet> State;
  typedef typename Hmm<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Forward algorithm ...";
  const int length     = subject.length();
  const int num_states = hmm.num_states();
  ForwardBackwardMatrices& m = *fbm;
  m.logp = 0.0;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", 0);
  for (int k = 0; k < num_states; ++k) {
    m.e[0][k] = pow(2.0, emission(hmm[k], subject, 0));
    m.f[0][k] = hmm[k].prior() * m.e[0][k];
    m.s[0] += m.f[0][k];

    LOG(DEBUG2) << strprintf("f[%i][%i] = %-7.2e", 0, k, m.f[0][k]);
  }
  // Rescale forward values
  double scale_fac = 1.0 / m.s[0];
  for (int k = 0; k < num_states; ++k)
    m.f[0][k] *= scale_fac;
  m.logp += log2(m.s[0]);

  // Recursion
  for (int i = 1; i < length; ++i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    m.s[i] = 0.0;

    for (int l = 0; l < num_states; ++l) {
      double f_il = 0.0;
      LOG(DEBUG2) << strprintf("f[%i][%i] = 0", i, l);

      for (ConstTransitionIter t_kl = hmm[l].in_transitions_begin();
           t_kl != hmm[l].in_transitions_end(); ++t_kl) {
        f_il += m.f[i-1][t_kl->state] * t_kl->weight;

        LOG(DEBUG3) << strprintf("f[%i][%i] += f[%i][%i]=%-7.2e * tr[%i][%i]=%-7.5f",
                                 i, l, i-1, t_kl->state, m.f[i-1][t_kl->state],
                                 t_kl->state, l, t_kl->weight);
      }

      m.e[i][l] = pow(2.0, emission(hmm[l], subject, i));
      f_il *= m.e[i][l];

      LOG(DEBUG3) << strprintf("f[%i][%i] *= e[%i][%i]=%-7.2e",
                               i, l, i, l, m.e[i][l]);

      m.f[i][l] = f_il;
      m.s[i] += f_il;
    }

    // Rescale forward values
    scale_fac = 1.0 / m.s[i];
    for (int l = 0; l < num_states; ++l)
      m.f[i][l] *= scale_fac;
    m.logp += log2(m.s[i]);
  }

  LOG(DEBUG) << strprintf("log(L) = %-7.2g", m.logp);
}

template< class Alphabet, template<class> class Subject >
void BackwardAlgorithm(const Hmm<Alphabet>& hmm,
                       const Subject<Alphabet>& subject,
                       ForwardBackwardMatrices* fbm) {
  typedef ContextProfileState<Alphabet> State;
  typedef typename Hmm<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Backward algorithm ...";
  const int length     = subject.length();
  const int num_states = hmm.num_states();
  ForwardBackwardMatrices& m = *fbm;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", length-1);
  for (int l = 0; l < num_states; ++l) {
    m.b[length-1][l] = 1.0;
    m.pp[length-1][l] = m.f[length-1][l] * m.b[length-1][l];
    LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", length-1, l, m.b[length-1][l]);
  }

  // Recursion
  for (int i = length-2; i >= 0; --i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    for (int k = 0; k < num_states; ++k) {
      double b_ik = 0.0;
      LOG(DEBUG2) << strprintf("b[%i][%i] = 0", i, k);

      for (ConstTransitionIter t_kl = hmm[k].out_transitions_begin();
           t_kl != hmm[k].out_transitions_end(); ++t_kl) {
        b_ik += t_kl->weight * m.e[i+1][t_kl->state] * m.b[i+1][t_kl->state];

        LOG(DEBUG3) << strprintf("b[%i][%i] += tr[%i][%i]=%-7.5f *"
                                 " e[%i][%i]=%-7.5f * f[%i][%i]=%-7.2e",
                                 i, k, t_kl->state, k, t_kl->weight,
                                 i+1, t_kl->state, m.e[i+1][t_kl->state],
                                 i+1, t_kl->state, m.b[i+1][t_kl->state]);
      }
      m.b[i][k] = b_ik / m.s[i+1];
      m.pp[i][k] = m.f[i][k] * m.b[i][k];
      LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", i, k, m.b[i][k]);
    }
  }
}

}  // namespace cs

#endif  // CS_FORWARD_BACKWARD_ALGORITHM_INL_H_
