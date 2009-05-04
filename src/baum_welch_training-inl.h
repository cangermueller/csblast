// Copyright 2009, Andreas Biegert

#ifndef SRC_BAUM_WELCH_TRAINING_INL_H_
#define SRC_BAUM_WELCH_TRAINING_INL_H_

#include "baum_welch_training.h"

#include <cmath>
#include <cstdio>

#include <string>

#include "context_profile-inl.h"
#include "forward_backward_algorithm.h"
#include "hmm-inl.h"
#include "log.h"
#include "progress_table.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template< class Alphabet,
          template<class A> class Subject >
BaumWelchTraining<Alphabet, Subject>::BaumWelchTraining(
    const BaumWelchOptions& opts,
    const data_vector& data,
    HMM<Alphabet>& hmm)
    : ExpectationMaximization<Alphabet, Subject>(opts, data),
      opts_(opts),
      hmm_(hmm),
      emission_(hmm.num_cols(), opts.weight_center, opts.weight_decay),
      transition_stats_(hmm.num_states(), hmm.num_states()),
      profile_stats_(),
      transition_stats_block_(hmm.num_states(), hmm.num_states()),
      profile_stats_block_() {
  init();
}

template< class Alphabet,
          template<class A> class Subject >
BaumWelchTraining<Alphabet, Subject>::BaumWelchTraining(
    const BaumWelchOptions& opts,
    const data_vector& data,
    HMM<Alphabet>& hmm,
    FILE* fout)
    : ExpectationMaximization<Alphabet, Subject>(data),
      opts_(opts),
      hmm_(hmm),
      emission_(hmm.num_cols(), opts.weight_center, opts.weight_decay),
      transition_stats_(hmm.num_states(), hmm.num_states()),
      profile_stats_(),
      transition_stats_block_(hmm.num_states(), hmm.num_states()),
      profile_stats_block_() {
  progress_table_ = new BaumWelchProgressTable<Alphabet, Subject>(this, fout);
  init();
}

template< class Alphabet,
          template<class A> class Subject >
BaumWelchTraining<Alphabet, Subject>::~BaumWelchTraining() {
  if (progress_table_) delete progress_table_;
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::expectation_step(
    const data_vector& block) {
  // Run forward and backward algorithm on each subject in current block
  for (typename data_vector::const_iterator bi = block.begin();
       bi != block.end(); ++bi) {
    ForwardBackwardMatrices fbm((*bi)->length(), hmm_.num_states());
    forward_backward_algorithm(hmm_, **bi, emission_, &fbm);
    add_contribution_to_priors(fbm);
    add_contribution_to_transitions(fbm);
    add_contribution_to_emissions(fbm, **bi);
    log_likelihood_ += fbm.log_likelihood / num_eff_cols_;

    if (progress_table_)
      progress_table_->print_progress(hmm_.num_states() * (**bi).length());
  }

  update_sufficient_statistics();
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::maximization_step() {
  const int num_states    = hmm_.num_states();
  const int num_cols      = hmm_.num_cols();
  const int alphabet_size = hmm_.alphabet_size();

  // Calculate normalization factor for priors
  float sum = 0.0f;
  for (int k = 0; k < num_states; ++k) sum += profile_stats_[k]->prior();
  float fac = 1.0f / sum;

  // Assign new priors and emission probabilities
  for (int k = 0; k < num_states; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_[k];

    hmm_[k].set_prior(p_k.prior() * fac);
    ContextProfile<Alphabet> tmp(p_k);
    if (normalize(&tmp)) {  // don't update profiles that did'n get evidence
      tmp.transform_to_logspace();
      for (int i = 0; i < num_cols; ++i)
        for (int a = 0; a < alphabet_size; ++a)
          hmm_[k][i][a] = tmp[i][a];
    }
  }

  // Calculate and assign new transition probabilities
  hmm_.clear_transitions();
  for (int k = 0; k < num_states; ++k) {
    sum = 0.0f;
    for (int l = 0; l < num_states; ++l) {
      if (transition_stats_.test(k,l)) {
        const float a_kl =
          transition_stats_[k][l] + opts_.transition_pc - 1.0f;
        if (a_kl > 0.0f) {
          transition_stats_[k][l] = a_kl;
          sum += a_kl;
        } else {
          transition_stats_.erase(k,l);
        }
      }
    }
    if (sum != 0.0f) {
      fac = 1.0f / sum;
      for (int l = 0; l < num_states; ++l) {
        if (transition_stats_.test(k,l)) {
          hmm_(k,l) = transition_stats_[k][l] * fac;
          LOG(DEBUG2) << strprintf("tr[%i][%i]=%-8.5f",
                                   k, l, static_cast<float>(hmm_(k,l)));
        }
      }
    }
  }

  // Increment iteration counter in HMM
  ++hmm_;
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::add_contribution_to_priors(
    const ForwardBackwardMatrices& m) {
  const int num_states = hmm_.num_states();

  for (int k = 0; k < num_states; ++k)
    profile_stats_block_[k]->set_prior(profile_stats_block_[k]->prior()
                                       + m.f[0][k] * m.b[0][k]);
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::add_contribution_to_transitions(
    const ForwardBackwardMatrices& m) {
  const int slen = m.f.num_rows();

  for (const_transition_iterator ti = hmm_.transitions_begin();
       ti != hmm_.transitions_end(); ++ti) {
    double w_kl = 0.0;
    for (int i = 0; i < slen-1; ++i) {
      w_kl += m.f[i][ti->from] * m.b[i+1][ti->to] * ti->probability
        * m.e[i+1][ti->to] / m.s[i+1];
    }

    if (w_kl != 0.0) {
      if (!transition_stats_block_.test(ti->from, ti->to))
        transition_stats_block_[ti->from][ti->to] = 0.0f;
      transition_stats_block_[ti->from][ti->to] =
        transition_stats_block_[ti->from][ti->to] + w_kl;
    }
  }
}

template< class Alphabet,
          template<class A> class Subject >
inline void BaumWelchTraining<Alphabet, Subject>::add_contribution_to_emissions(
    const ForwardBackwardMatrices& m,
    const CountProfile<Alphabet>& c) {
  const int slen       = c.length();
  const int num_states = hmm_.num_states();

  for (int k = 0; k < num_states; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_block_[k];
    const int ci = p_k.center();

    for (int i = 0; i < slen; ++i) {
      const int beg = std::max(0, i - ci);
      const int end = std::min(c.num_cols() - 1, i + ci);

      for(int h = beg; h <= end; ++h) {
        const int j = h - i + ci;
        const int alphabet_size = p_k.alphabet_size();
        for (int a = 0; a < alphabet_size; ++a) {
          p_k[j][a] += c[h][a] * m.f[i][k] * m.b[i][k];
        }
      }
    }
  }
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::add_contribution_to_emissions(
    const ForwardBackwardMatrices& m,
    const Sequence<Alphabet>& s) {
  const int slen       = s.length();
  const int num_states = hmm_.num_states();

  for (int k = 0; k < num_states; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_block_[k];
    const int ci = p_k.center();

    for (int i = 0; i < slen; ++i) {
      const int beg = std::max(0, i - ci);
      const int end = std::min(s.length() - 1, i + ci);

      for(int h = beg; h <= end; ++h) {
        const int j = h - i + ci;
        p_k[j][s[h]] += m.f[i][k] * m.b[i][k];
      }
    }
  }
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::update_sufficient_statistics() {
  const float gamma       = 1.0f - epsilon_;
  const int num_states    = hmm_.num_states();
  const int num_cols      = hmm_.num_cols();
  const int alphabet_size = hmm_.alphabet_size();

  // Update transition statistics
  for (const_transition_iterator ti = hmm_.transitions_begin();
       ti != hmm_.transitions_end(); ++ti) {
    if (transition_stats_block_.test(ti->from, ti->to)) {
      if (!transition_stats_.test(ti->from, ti->to))
        transition_stats_[ti->from][ti->to] = 0.0f;
      transition_stats_[ti->from][ti->to] =
        gamma * transition_stats_[ti->from][ti->to] +
        transition_stats_block_[ti->from][ti->to];
    }
  }
  transition_stats_block_.clear();

  // Update priors and emissions statistics
  for (int k = 0; k < num_states; ++k) {
    ContextProfile<Alphabet>& p_block = *profile_stats_block_[k];
    ContextProfile<Alphabet>& p       = *profile_stats_[k];

    p.set_prior(p.prior() * gamma + p_block.prior());
    for (int j = 0; j < num_cols; ++j) {
      for (int a = 0; a < alphabet_size; ++a) {
        p[j][a] = gamma * p[j][a]  + p_block[j][a];
      }
    }
    reset(&p_block);
  }
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchTraining<Alphabet, Subject>::init() {
  // Create profiles for global and block-level sufficient statistics
  for (int k = 0; k < hmm_.num_states(); ++k) {
    profile_stats_.push_back(
        shared_ptr< ContextProfile<Alphabet> >(
            new ContextProfile<Alphabet>(k, hmm_.num_cols())));
    profile_stats_block_.push_back(
        shared_ptr< ContextProfile<Alphabet> >(
            new ContextProfile<Alphabet>(k, hmm_.num_cols())));
  }

  // Compute total number of data columns
  int num_cols = 0;
  for (typename data_vector::const_iterator di = data_.begin();
       di != data_.end(); ++di)
    num_cols += (**di).length();
  if (progress_table_)
    progress_table_->set_total_work(hmm_.num_states() * num_cols);

  // Set number of effective columsn for log-likelihood calculation
  num_eff_cols_ = emission_.SumWeights() * num_cols;
}

template< class Alphabet,
          template<class A> class Subject >
inline bool BaumWelchTraining<Alphabet, Subject>::terminate() const {
  if (scan_ < opts_.min_scans)
    return false;
  else if (scan_ >= opts_.max_scans)
    return true;
  else if (opts_.max_connectivity == 0)
    return fabs(log_likelihood_change()) <= opts_.log_likelihood_change;
  else
    return
      fabs(log_likelihood_change()) <= opts_.log_likelihood_change
      && hmm_.connectivity() <= opts_.max_connectivity;
}


template< class Alphabet,
          template<class A> class Subject >
BaumWelchProgressTable<Alphabet, Subject>::BaumWelchProgressTable(
    const BaumWelchTraining<Alphabet, Subject>* training,
    FILE* fout,
    int width)
    : ProgressTable(fout, width),
      training_(training) {}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::print_header() {
  fprintf(fout_, "%-4s %4s %6s %4s %7s  %-30s  %9s  %8s\n",
          "Scan", "Itrs", "Conn", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
  fputs(std::string(82, '-').c_str(), fout_);
  fputc('\n', fout_);
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::print_row_begin() {
  reset();
  fprintf(fout_, "%-4i %4i %6.1f %4i %7.4f  ", training_->scan(),
          training_->iterations(), training_->hmm_.connectivity(),
          training_->num_blocks(), training_->epsilon());
  fflush(fout_);
}

template< class Alphabet,
          template<class A> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::print_row_end() {
  if (training_->scan() == 1)
    fprintf(fout_, "  %9.5f\n", training_->log_likelihood());
  else
    fprintf(fout_, "  %9.5f  %+8.5f\n", training_->log_likelihood(),
            training_->log_likelihood_change());
  fflush(fout_);
}

}  // namespace cs

#endif  // SRC_BAUM_WELCH_TRAINING_INL_H_
