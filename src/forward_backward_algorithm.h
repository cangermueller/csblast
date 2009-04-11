// Copyright 2009, Andreas Biegert

#ifndef SRC_FORWARD_BACKWARD_ALGORITHM_H_
#define SRC_FORWARD_BACKWARD_ALGORITHM_H_

#include <valarray>

#include "context_profile-inl.h"
#include "emitter-inl.h"
#include "hmm-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

// POD with forward backward matrices, emission probs, and log-likelihood.
struct ForwardBackwardMatrices {
  typedef double value_type;

  ForwardBackwardMatrices()
      : log_likelihood(0.0) { }

  ForwardBackwardMatrices(int slen, int nstates)
      : f(slen, nstates, 0.0),
        b(slen, nstates, 0.0),
        e(slen, nstates, 1.0),
        s(0.0, slen),
        log_likelihood(0.0) { }

  friend std::ostream& operator<< (std::ostream& out,
                                   const ForwardBackwardMatrices& m) {
    out << "Forward matrix f[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%6.4f  ", m.f[i][k]);
      }
      out << std::endl;
    }
    out << "Backward matrix b[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%6.4f  ", m.b[i][k]);
      }
      out << std::endl;
    }
    out << "Forward row-sums before scaling s[i]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      out << strprintf("%6.2f  ", m.s[i]);
    }
    out << "\nEmitter probabilities e[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%6.4f  ", m.e[i][k]);
      }
      out << std::endl;
    }
    out << "Posterior probability matrix p[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        value_type p = 100.0 * m.f[i][k] * m.b[i][k];
        out << strprintf("%6.2f  ", p);
      }
      out << std::endl;
    }
    out << strprintf("Log-likelihood = %-7.2g\n", m.log_likelihood);
    return out;
  }

  // forward matrix
  matrix<value_type> f;
  // backward matrix
  matrix<value_type> b;
  // emission probability matrix
  matrix<value_type> e;
  // forward matrix column sum s(n) before normalization
  std::valarray<value_type> s;
  // likelihood P(x)
  value_type log_likelihood;
};


// Forward-Backward algorithm encapsulation.
template< class Alphabet,
          template<class A> class Subject >
void forward_backward_algorithm(const HMM<Alphabet>& hmm,
                                const Subject<Alphabet>& subject,
                                const Emitter<Alphabet>& emitter,
                                ForwardBackwardMatrices* fbm) {
  LOG(DEBUG) << "Running forward-backward algorithm ...";
  LOG(DEBUG1) << hmm;
  LOG(DEBUG1) << subject;

  forward_algorithm(hmm, subject, emitter, fbm);
  backward_algorithm(hmm, subject, fbm);
}

template< class Alphabet,
          template<class A> class Subject >
void forward_algorithm(const HMM<Alphabet>& hmm,
                       const Subject<Alphabet>& subject,
                       const Emitter<Alphabet>& emitter,
                       ForwardBackwardMatrices* fbm) {
  typedef typename ForwardBackwardMatrices::value_type value_type;
  typedef typename HMM<Alphabet>::const_state_iterator const_state_iterator;
  typedef typename State<Alphabet>::const_transition_iterator const_transition_iterator;

  LOG(DEBUG1) << "Forward algorithm ...";
  const int length     = subject.length();
  const int num_states = hmm.num_states();
  ForwardBackwardMatrices& m = *fbm;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", 0);
  for (int k = 0; k < num_states; ++k) {
    m.e[0][k] = fast_pow2(emitter(hmm[k], subject, 0));
    // m.e[0][k] = pow(2.0, emitter(hmm[k], subject, 0));
    LOG(DEBUG1) << strprintf("%-7.2e  %i  %i  %7.3f",
                             emitter(hmm[k], subject, 0),
                             hmm.states_logspace() ? 1 : 0,
                             hmm[k].logspace() ? 1 : 0,
                             hmm[k][0][0]);
    m.f[0][k] = hmm[k].prior() * m.e[0][k];
    LOG(DEBUG2) << strprintf("f[%i][%i] = %-7.2e", 0, k, m.f[0][k]);
  }
  m.log_likelihood = 0.0;

  // Recursion
  for (int i = 1; i < length; ++i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    m.s[i] = 0.0;

    for (int l = 0; l < num_states; ++l) {
      value_type f_il = 0.0;
      LOG(DEBUG2) << strprintf("f[%i][%i] = 0", i, l);

      for (const_transition_iterator t_kl = hmm[l].in_transitions_begin();
           t_kl != hmm[l].in_transitions_end(); ++t_kl) {
        f_il += m.f[i-1][t_kl->state] * t_kl->probability;

        LOG(DEBUG3) <<
          strprintf("f[%i][%i] += f[%i][%i]=%-7.2e * tr[%i][%i]=%-7.5f",
                    i, l, i-1, t_kl->state, m.f[i-1][t_kl->state],
                    t_kl->state, l, t_kl->probability);
      }

      m.e[i][l] = fast_pow2(emitter(hmm[l], subject, i));
      //m.e[i][l] = pow(2.0, emitter(hmm[l], subject, i));
      f_il *= m.e[i][l];
      LOG(DEBUG3) << strprintf("f[%i][%i] *= e[%i][%i]=%-7.5f",
                               i, l, i, l, m.e[i][l]);

      m.f[i][l] = f_il;
      m.s[i] += f_il;
      LOG(DEBUG2) << strprintf("f[%i][%i] = %-7.2e", i, l, f_il);
    }

    // Rescale forward values
    value_type scale_fac = 1.0 / m.s[i];
    for (int l = 0; l < num_states; ++l)
      m.f[i][l] *= scale_fac;
    m.log_likelihood += fast_log2(m.s[i]);
    //m.log_likelihood += log2(m.s[i]);
  }

  LOG(DEBUG) << strprintf("log(L) = %-7.2g", m.log_likelihood);
}

template< class Alphabet,
          template<class A> class Subject >
void backward_algorithm(const HMM<Alphabet>& hmm,
                        const Subject<Alphabet>& subject,
                        ForwardBackwardMatrices* fbm) {
  typedef typename ForwardBackwardMatrices::value_type value_type;
  typedef typename HMM<Alphabet>::const_state_iterator const_state_iterator;
  typedef typename State<Alphabet>::const_transition_iterator const_transition_iterator;

  LOG(DEBUG1) << "Backward algorithm ...";
  const int length     = subject.length();
  const int num_states = hmm.num_states();
  ForwardBackwardMatrices& m = *fbm;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", length-1);
  for (int l = 0; l < num_states; ++l) {
    m.b[length-1][l] = 1.0;
    LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", length-1, l, m.b[length-1][l]);
  }

  // Recursion
  for (int i = length-2; i >= 0; --i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    for (int k = 0; k < num_states; ++k) {
      value_type b_ik = 0.0;
      LOG(DEBUG2) << strprintf("b[%i][%i] = 0", i, k);

      for (const_transition_iterator t_kl = hmm[k].out_transitions_begin();
           t_kl != hmm[k].out_transitions_end(); ++t_kl) {
        b_ik += t_kl->probability * m.e[i+1][t_kl->state] * m.b[i+1][t_kl->state];

        LOG(DEBUG3) << strprintf("b[%i][%i] += tr[%i][%i]=%-7.5f *"
                                 " e[%i][%i]=%-7.5f * f[%i][%i]=%-7.2e",
                                 i, k, t_kl->state, k, t_kl->probability,
                                 i+1, t_kl->state, m.e[i+1][t_kl->state],
                                 i+1, t_kl->state, m.b[i+1][t_kl->state]);
      }
      m.b[i][k] = b_ik / m.s[i+1];
      LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", i, k, m.b[i][k]);
    }
  }
}

}  // namespace cs

#endif  // SRC_FORWARD_BACKWARD_ALGORITHM_H_
