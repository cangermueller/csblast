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
      : alpha(slen, nstates, 0.0),
        beta(slen, nstates, 0.0),
        alpha_pc(slen, nstates, 0.0),
        beta_pc(slen, nstates, 0.0),
        context_score(slen, nstates, 0.0),
        pc_prob(slen, nstates, 0.0),
        alpha_sum(0.0, slen),
        alpha_sum_pc(0.0, slen),
        log_likelihood(0.0),
        log_likelihood_pc(0.0) {}

  friend std::ostream& operator<< (std::ostream& out,
                                   const SumProductMatrices& m) {
    out << "Alpha matrix alpha[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%7.5f  ", m.alpha[i][k]);
      }
      out << std::endl;
    }
    out << "Beta matrix b[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%7.5f  ", m.beta[i][k]);
      }
      out << std::endl;
    }
    out << "Alpha row-sums before scaling alpha_sum[i]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      out << strprintf("%7.5f  ", m.alpha_sum[i]);
    }
    out << "\nContext-score matrix context_score[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        out << strprintf("%7.5f  ", m.context_score[i][k]);
      }
      out << std::endl;
    }
    out << "Posterior probabilities pp[i][k]:" << std::endl;
    for (int i = 0; i < m.f.num_rows(); ++i) {
      for (int k = 0; k < m.f.num_cols(); ++k) {
        double p = 100.0 * m.alpha[i][k] * m.beta[i][k];
        out << strprintf("%7.5f  ", p);
      }
      out << std::endl;
    }
    out << strprintf("Conditional log-likelihood = %-10.5f\n", m.log_likelihood);
    return out;
  }

  // Matrix with alpha messages from forward pass
  matrix<double> alpha;
  // Matrix with beta messages from backward pass
  matrix<double> beta;
  // Matrix with alpha messages from forward pass WITH pseudocount emissions
  matrix<double> alpha_pc;
  // Matrix with beta messages from backward pass WITH pseudocount emissions
  matrix<double> beta_pc;
  // Context-score matrix
  matrix<double> context_score;
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

template<class Alphabet>
inline double CalculateContextScore<Alphabet>(
    const ContextWeightState<Alphabet>& state,
    const CountProfile<Alphabet>& count_profile,
    int index) const {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();
  const int center = state.center();
  const int beg = std::max(0, index - center);
  const int end = std::min(count_profile.num_cols() - 1, index + center);
  double rv = 0.0;

  for(int i = beg; i <= end; ++i) {
    const int j = i - index + center;
    for (int a = 0; a < alphabet_size; ++a)
      rv += count_profile.counts(i, a) * state[j][a];
  }

  return rv;
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
  SumProductBackwardPass(hmm, subject, spm);
}

template< class Alphabet, template<class> class Subject >
void SumProductForwardPass(const CRF<Alphabet>& crf,
                           const Subject<Alphabet>& subject,
                           SumProductMatrices* spm) {
  typedef CRF<Alphabet>::State State;
  typedef typename CRF<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product forward pass ...";
  const int length     = subject.length();
  const int num_states = crf.num_states();
  SumProductMatrices& m = *spm;
  m.log_likelihood    = 0.0;
  m.log_likelihood_pc = 0.0;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", 0);
  for (int k = 0; k < num_states; ++k) {
    m.context_score[0][k] = pow(2.0, CalculateContextScore(crf[k], subject, 0));
    m.alpha[0][k] = m.context_score[0][k];
    m.alpha_sum[0] += m.alpha[0][k];

    LOG(DEBUG2) << strprintf("alpha[%i][%i] = %-7.2e", 0, k, m.alpha[0][k]);
  }

  // Rescale alpha values
  double scale_fac = 1.0 / m.alpha_sum[0];
  for (int k = 0; k < num_states; ++k)
    m.alpha[0][k] *= scale_fac;
  m.log_likelihood += log2(m.alpha_sum[0]);

  // Recursion
  for (int i = 1; i < length; ++i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    m.alpha_sum[i] = 0.0;

    for (int l = 0; l < num_states; ++l) {
      double alpha_il = 0.0;
      LOG(DEBUG2) << strprintf("alpha[%i][%i] = 0", i, l);

      for (ConstTransitionIter t_kl = crf[l].in_transitions_begin();
           t_kl != crf[l].in_transitions_end(); ++t_kl) {
        alpha_il += m.alpha[i-1][t_kl->state] * t_kl->weight;

        LOG(DEBUG3) << strprintf("alpha[%i][%i] += alpha[%i][%i]=%-7.2e * "
                                 "tr[%i][%i]=%-7.5f", i, l, i-1, t_kl->state,
                                 m.f[i-1][t_kl->state], t_kl->state, l,
                                 t_kl->weight);
      }

      m.context_score[i][l] = pow(2.0, CalculateContextScore(crf[l], subject, i));
      alpha_il *= m.context_score[i][l];

      LOG(DEBUG3) << strprintf("alpha[%i][%i] *= context_score[%i][%i]=%-7.2e",
                               i, l, i, l, m.context_score[i][l]);

      m.alpha[i][l] = alpha_il;
      m.alpha_sum[i] += alpha_il;
    }

    // Rescale forward values
    scale_fac = 1.0 / m.alpha_sum[i];
    for (int l = 0; l < num_states; ++l)
      m.alpha[i][l] *= scale_fac;
    m.log_likelihood += log2(m.alpha_sum[i]);
  }

  LOG(DEBUG) << strprintf("log(L) = %-10.5f", m.log_likelihood);
}

template< class Alphabet, template<class> class Subject >
void SumProductBackwardPass(const CRF<Alphabet>& crf,
                            const Subject<Alphabet>& subject,
                            SumProductMatrices* spm) {
  typedef CRF<Alphabet>::State State;
  typedef typename CRF<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product backward pass ...";
  const int length     = subject.length();
  const int num_states = crf.num_states();
  SumProductMatrices& m = *spm;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", length-1);
  for (int l = 0; l < num_states; ++l) {
    m.beta[length-1][l] = 1.0;
    LOG(DEBUG2) << strprintf("b[%i][%i] = %-7.2e", length-1, l,
                             m.beta[length-1][l]);
  }

  // Recursion
  for (int i = length-2; i >= 0; --i) {
    LOG(DEBUG1) << strprintf("i=%i", i);
    for (int k = 0; k < num_states; ++k) {
      double beta_ik = 0.0;
      LOG(DEBUG2) << strprintf("beta[%i][%i] = 0", i, k);

      for (ConstTransitionIter t_kl = crf[k].out_transitions_begin();
           t_kl != crf[k].out_transitions_end(); ++t_kl) {
        beta_ik +=
          t_kl->weight * m.context_score[i+1][t_kl->state] *
          m.beta[i+1][t_kl->state];

        LOG(DEBUG3) << strprintf("beta[%i][%i] += tr[%i][%i]=%-7.5f * "
                                 "context_score[%i][%i]=%-7.5f * "
                                 "beta[%i][%i]=%-7.2e", i, k, t_kl->state,
                                 k, t_kl->weight, i+1, t_kl->state,
                                 m.context_score[i+1][t_kl->state], i+1,
                                 t_kl->state, m.beta[i+1][t_kl->state]);
      }
      m.beta[i][k] = beta_ik / m.alpha_sum[i+1];
      LOG(DEBUG2) << strprintf("beta[%i][%i] = %-7.2e", i, k, m.beta[i][k]);
    }
  }
}

}  // namespace cs

#endif  // SRC_SUM_PRODUCT_ALGORITHM_H_
