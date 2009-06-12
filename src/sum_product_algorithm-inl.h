// Copyright 2009, Andreas Biegert

#ifndef SRC_SUM_PRODUCT_ALGORITHM_INL_H_
#define SRC_SUM_PRODUCT_ALGORITHM_INL_H_

#include <valarray>

#include "context_weight_state-inl.h"
#include "crf-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

// POD needed for output of sum product algorithm.
struct SumProductMatrices {
  SumProductMatrices(int slen, int nstates)
      : subject_length(slen),
        num_states(nstates),
        alpha(slen, nstates, 0.0),
        beta(slen, nstates, 0.0),
        alpha_pc(slen, nstates, 0.0),
        beta_pc(slen, nstates, 0.0),
        context_prob(slen, nstates, 0.0),
        pc_prob(slen, nstates, 0.0),
        alpha_sum(0.0, slen),
        alpha_sum_pc(0.0, slen),
        log_likelihood(0.0),
        log_likelihood_pc(0.0) {}

  friend std::ostream& operator<< (std::ostream& out,
                                   const SumProductMatrices& m) {
    out << "Alpha matrix alpha[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%7.5f  ", m.alpha[i][k]);
      }
      out << std::endl;
    }
    out << "Beta matrix b[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%7.5f  ", m.beta[i][k]);
      }
      out << std::endl;
    }
    out << "Alpha row-sums before scaling alpha_sum[i]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      out << strprintf("%9.3g  ", m.alpha_sum[i]);
    }
    out << "\nContext probability matrix context_prob[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%8.3g  ", m.context_prob[i][k]);
      }
      out << std::endl;
    }
    out << "\nPseudocount probability matrix pc_prob[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        out << strprintf("%8.3g  ", m.pc_prob[i][k]);
      }
      out << std::endl;
    }
    out << "Posterior probabilities pp[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        double p = 100.0 * m.alpha[i][k] * m.beta[i][k];
        out << strprintf("%6.2f  ", p);
      }
      out << std::endl;
    }
    out << "Posterior probabilities WITH pseudocounts pp[i][k]:" << std::endl;
    for (int i = 0; i < m.subject_length; ++i) {
      for (int k = 0; k < m.num_states; ++k) {
        double p = 100.0 * m.alpha_pc[i][k] * m.beta_pc[i][k];
        out << strprintf("%6.2f  ", p);
      }
      out << std::endl;
    }
    out << strprintf("Conditional log-likelihood = %-10.5f\n", m.log_likelihood);
    out << strprintf("Conditional log-likelihood WITH pseudocounts = %-10.5f\n",
                     m.log_likelihood_pc);
    return out;
  }

  // length of the subject whose Fwd-Bwd information is stored
  const int subject_length;
  // Number of states in the CRF against which the subject is aligned
  const int num_states;
  // Matrix with alpha messages from forward pass
  matrix<double> alpha;
  // Matrix with beta messages from backward pass
  matrix<double> beta;
  // Matrix with alpha messages from forward pass WITH pseudocount emissions
  matrix<double> alpha_pc;
  // Matrix with beta messages from backward pass WITH pseudocount emissions
  matrix<double> beta_pc;
  // Context-probability matrix
  matrix<double> context_prob;
  // Matrix with pseudocount emission probabilites
  matrix<double> pc_prob;
  // Row sums of alpha matrix before re-scaling
  std::valarray<double> alpha_sum;
  // Row sums of alpha_pc matrix before re-scaling
  std::valarray<double> alpha_sum_pc;
  // Conditional log-likelihood computed from forward pass
  double log_likelihood;
  // Conditional log-likelihood computed from forward pass WITH pseudocount
  // emissions
  double log_likelihood_pc;
};


// Calculates the probability of the context window centered around position i
// in the count profile given the context weights in the CRF state.
template<class Alphabet>
inline double ContextProb(const ContextWeightState<Alphabet>& state,
                          const CountProfile<Alphabet>& count_profile,
                          int idx) {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();
  const int center = state.center();
  const int beg = std::max(0, idx - center);
  const int end = std::min(count_profile.num_cols() - 1, idx + center);
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - idx + center;
    for (int a = 0; a < alphabet_size; ++a)
      rv += count_profile.counts(i, a) * state[j][a];
  }

  return pow(2.0, rv);
}

// Calculates the probability of counts in column i of given count profile to
// be emitted by the central pseudocount column in the provided CRF state.
template<class Alphabet>
inline double PseudocountProb(const ContextWeightState<Alphabet>& state,
                              const CountProfile<Alphabet>& count_profile,
                              int i) {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();
  double pc_sum = 0.0;
  for (int a = 0; a < alphabet_size; ++a)
    pc_sum += pow(2.0, state(a));
  pc_sum = log2(pc_sum);

  double rv = 0.0;
  for (int a = 0; a < alphabet_size; ++a)
    rv += count_profile.counts(i, a) * (state(a) - pc_sum);

  return pow(2.0, rv);
}


// Sum-product algorithm encapsulation.
template< class Alphabet, template<class> class Subject >
void SumProductAlgorithm(const CRF<Alphabet>& crf,
                         const Subject<Alphabet>& subject,
                         SumProductMatrices* spm) {
  LOG(DEBUG) << "Running sum-product algorithm ...";
  LOG(DEBUG1) << crf;
  LOG(DEBUG1) << subject;

  SumProductForwardPass(crf, subject, spm);
  SumProductBackwardPass(crf, subject, spm);
}

template< class Alphabet, template<class> class Subject >
void SumProductForwardPass(const CRF<Alphabet>& crf,
                           const Subject<Alphabet>& subject,
                           SumProductMatrices* spm) {
  typedef ContextWeightState<Alphabet> State;
  typedef typename CRF<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product forward pass ...";
  const int length     = subject.length();
  const int num_states = crf.num_states();
  double scale_fac = 1.0;
  SumProductMatrices& m = *spm;
  m.log_likelihood    = 0.0;
  m.log_likelihood_pc = 0.0;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", 0);
  for (int k = 0; k < num_states; ++k) {
    m.context_prob[0][k] = ContextProb(crf[k], subject, 0);
    m.pc_prob[0][k]      = PseudocountProb(crf[k], subject, 0);

    m.alpha[0][k] = m.context_prob[0][k];
    m.alpha_sum[0] += m.alpha[0][k];

    m.alpha_pc[0][k] = m.context_prob[0][k] * m.pc_prob[0][k];;
    m.alpha_sum_pc[0] += m.alpha_pc[0][k];

    LOG(DEBUG2) << strprintf("alpha[%i][%i] = %-7.2e", 0, k, m.alpha[0][k]);
  }
  m.log_likelihood    += log2(m.alpha_sum[0]);
  m.log_likelihood_pc += log2(m.alpha_sum_pc[0]);

  // Rescale alpha values
  scale_fac = 1.0 / m.alpha_sum[0];
  for (int k = 0; k < num_states; ++k)
    m.alpha[0][k] *= scale_fac;
  scale_fac = 1.0 / m.alpha_sum_pc[0];
  for (int k = 0; k < num_states; ++k)
    m.alpha_pc[0][k] *= scale_fac;

  // Recursion
  for (int i = 1; i < length; ++i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    m.alpha_sum[i]    = 0.0;
    m.alpha_sum_pc[i] = 0.0;

    for (int l = 0; l < num_states; ++l) {
      double alpha_il    = 0.0;
      double alpha_il_pc = 0.0;

      LOG(DEBUG2) << strprintf("alpha[%i][%i] = 0", i, l);

      for (ConstTransitionIter t_kl = crf[l].in_transitions_begin();
           t_kl != crf[l].in_transitions_end(); ++t_kl) {
        alpha_il    += m.alpha[i-1][t_kl->state]    * t_kl->weight;
        alpha_il_pc += m.alpha_pc[i-1][t_kl->state] * t_kl->weight;

        LOG(DEBUG3) << strprintf("alpha[%i][%i] += alpha[%i][%i]=%-7.2e * "
                                 "tr[%i][%i]=%-7.5f", i, l, i-1, t_kl->state,
                                 m.alpha[i-1][t_kl->state], t_kl->state, l,
                                 t_kl->weight);
      }

      m.context_prob[i][l] = ContextProb(crf[l], subject, i);
      m.pc_prob[i][l]      = PseudocountProb(crf[l], subject, i);

      alpha_il    *= m.context_prob[i][l];
      alpha_il_pc *= m.context_prob[i][l] * m.pc_prob[i][l];

      LOG(DEBUG3) << strprintf("alpha[%i][%i] *= context_prob[%i][%i]=%-7.2e",
                               i, l, i, l, m.context_prob[i][l]);

      m.alpha[i][l]    = alpha_il;
      m.alpha_pc[i][l] = alpha_il_pc;

      m.alpha_sum[i]    += alpha_il;
      m.alpha_sum_pc[i] += alpha_il_pc;
    }
    m.log_likelihood    += log2(m.alpha_sum[i]);
    m.log_likelihood_pc += log2(m.alpha_sum_pc[i]);

    // Rescale alpha values
    scale_fac = 1.0 / m.alpha_sum[i];
    for (int l = 0; l < num_states; ++l)
      m.alpha[i][l] *= scale_fac;
    scale_fac = 1.0 / m.alpha_sum_pc[i];
    for (int l = 0; l < num_states; ++l)
      m.alpha_pc[i][l] *= scale_fac;
  }
  LOG(DEBUG) << strprintf("log(L) = %-10.5f", m.log_likelihood);
}

template< class Alphabet, template<class> class Subject >
void SumProductBackwardPass(const CRF<Alphabet>& crf,
                            const Subject<Alphabet>& subject,
                            SumProductMatrices* spm) {
  typedef ContextWeightState<Alphabet> State;
  typedef typename CRF<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product backward pass ...";
  const int length     = subject.length();
  const int num_states = crf.num_states();
  SumProductMatrices& m = *spm;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", length-1);
  for (int l = 0; l < num_states; ++l) {
    m.beta[length-1][l]    = 1.0;
    m.beta_pc[length-1][l] = 1.0;

    LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", length-1, l,
                             m.beta[length-1][l]);
  }

  // Recursion
  for (int i = length-2; i >= 0; --i) {
    LOG(DEBUG1) << strprintf("i=%i", i);

    for (int k = 0; k < num_states; ++k) {
      double beta_ik    = 0.0;
      double beta_ik_pc = 0.0;
      LOG(DEBUG2) << strprintf("beta[%i][%i] = 0", i, k);

      for (ConstTransitionIter t_kl = crf[k].out_transitions_begin();
           t_kl != crf[k].out_transitions_end(); ++t_kl) {
        beta_ik += t_kl->weight * m.context_prob[i+1][t_kl->state] *
          m.beta[i+1][t_kl->state];
        beta_ik_pc += t_kl->weight * m.context_prob[i+1][t_kl->state] *
          m.pc_prob[i+1][t_kl->state] * m.beta_pc[i+1][t_kl->state];

        LOG(DEBUG3) << strprintf("beta[%i][%i] += tr[%i][%i]=%-7.5f * "
                                 "context_prob[%i][%i]=%-7.5f * "
                                 "beta[%i][%i]=%-7.2e", i, k, t_kl->state,
                                 k, t_kl->weight, i+1, t_kl->state,
                                 m.context_prob[i+1][t_kl->state], i+1,
                                 t_kl->state, m.beta[i+1][t_kl->state]);
      }
      m.beta[i][k]    = beta_ik / m.alpha_sum[i+1];
      m.beta_pc[i][k] = beta_ik_pc / m.alpha_sum_pc[i+1];

      LOG(DEBUG2) << strprintf("beta[%i][%i] = %-7.2e", i, k, m.beta[i][k]);
    }
  }
}

}  // namespace cs

#endif  // SRC_SUM_PRODUCT_ALGORITHM_H_
