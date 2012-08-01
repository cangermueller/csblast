/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_SUM_PRODUCT_ALGORITHM_INL_H_
#define CS_SUM_PRODUCT_ALGORITHM_INL_H_

#include <algorithm>
#include <valarray>

#include "crf_state-inl.h"
#include "count_profile-inl.h"
#include "crf-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "utils.h"

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
        logp(0.0),
        logp_pc(0.0) {}

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
    out << strprintf("Conditional log-likelihood = %-10.5f\n", m.logp);
    out << strprintf("Conditional log-likelihood WITH pseudocounts = %-10.5f\n",
                     m.logp_pc);
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
  double logp;
  // Conditional log-likelihood computed from forward pass WITH pseudocount
  // emissions
  double logp_pc;
};


// Calculates the probability of the context window centered around position i
// in the count profile given the context weights in the CRF state.
template<class Alphabet>
inline double ContextProb(const CrfState<Alphabet>& state,
                          const CountProfile<Alphabet>& count_profile,
                          int idx) {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();
  const int center = state.center();
  const int beg = MAX(0, idx - center);
  const int end = MIN(count_profile.num_cols() - 1, idx + center);
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
inline double PseudocountProb(const CrfState<Alphabet>& state,
                              const CountProfile<Alphabet>& count_profile,
                              int i) {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();

  double rv = 0.0;
  for (int a = 0; a < alphabet_size; ++a)
    // FIXME: is the following line really correct -> check formula!!!
    rv += count_profile.counts(i, a) * log2(state.pc(a));

  return pow(2.0, rv);
}


// Sum-product algorithm encapsulation.
template< class Alphabet, template<class> class Subject >
void SumProductAlgorithm(const Crf<Alphabet>& crf,
                         const Subject<Alphabet>& subject,
                         const CountProfile<Alphabet>& pred,
                         SumProductMatrices* spm) {
  LOG(DEBUG) << "Running sum-product algorithm ...";
  LOG(DEBUG1) << crf;
  LOG(DEBUG1) << subject;

  SumProductForwardPass(crf, subject, pred, spm);
  SumProductBackwardPass(crf, spm);
}

template< class Alphabet, template<class> class Subject >
void SumProductForwardPass(const Crf<Alphabet>& crf,
                           const Subject<Alphabet>& subject,
                           const CountProfile<Alphabet>& pred,
                           SumProductMatrices* spm) {
  typedef CrfState<Alphabet> State;
  typedef typename Crf<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product forward pass ...";
  const int length     = spm->subject_length;
  const int num_states = crf.num_states();
  double scale_fac = 1.0;
  SumProductMatrices& m = *spm;
  m.logp    = 0.0;
  m.logp_pc = 0.0;

  // Initialization
  LOG(DEBUG1) << strprintf("i=%i", 0);
  for (int k = 0; k < num_states; ++k) {
    m.context_prob[0][k] = ContextProb(crf[k], subject, 0);
    m.pc_prob[0][k]      = PseudocountProb(crf[k], pred, 0);

    m.alpha[0][k] = m.context_prob[0][k];
    m.alpha_sum[0] += m.alpha[0][k];

    m.alpha_pc[0][k] = m.context_prob[0][k] * m.pc_prob[0][k];;
    m.alpha_sum_pc[0] += m.alpha_pc[0][k];

    LOG(DEBUG2) << strprintf("alpha[%i][%i] = %-7.2e", 0, k, m.alpha[0][k]);
  }
  m.logp    += log2(m.alpha_sum[0]);
  m.logp_pc += log2(m.alpha_sum_pc[0]);

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
      m.pc_prob[i][l]      = PseudocountProb(crf[l], pred, i);

      alpha_il    *= m.context_prob[i][l];
      alpha_il_pc *= m.context_prob[i][l] * m.pc_prob[i][l];

      LOG(DEBUG3) << strprintf("alpha[%i][%i] *= context_prob[%i][%i]=%-7.2e",
                               i, l, i, l, m.context_prob[i][l]);

      m.alpha[i][l]    = alpha_il;
      m.alpha_pc[i][l] = alpha_il_pc;

      m.alpha_sum[i]    += alpha_il;
      m.alpha_sum_pc[i] += alpha_il_pc;
    }
    m.logp    += log2(m.alpha_sum[i]);
    m.logp_pc += log2(m.alpha_sum_pc[i]);

    // Rescale alpha values
    scale_fac = 1.0 / m.alpha_sum[i];
    for (int l = 0; l < num_states; ++l)
      m.alpha[i][l] *= scale_fac;
    scale_fac = 1.0 / m.alpha_sum_pc[i];
    for (int l = 0; l < num_states; ++l)
      m.alpha_pc[i][l] *= scale_fac;
  }
  LOG(DEBUG) << strprintf("log(L) = %-10.5f", m.logp);
}

template<class Alphabet>
void SumProductBackwardPass(const Crf<Alphabet>& crf,
                            SumProductMatrices* spm) {
  typedef CrfState<Alphabet> State;
  typedef typename Crf<Alphabet>::ConstStateIter ConstStateIter;
  typedef typename State::ConstTransitionIter ConstTransitionIter;

  LOG(DEBUG1) << "Sum-product backward pass ...";
  const int length     = spm->subject_length;
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

#endif  // CS_SUM_PRODUCT_ALGORITHM_H_
