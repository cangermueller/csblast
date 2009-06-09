// Copyright 2009, Andreas Biegert

#ifndef SRC_BAUM_WELCH_TRAINING_INL_H_
#define SRC_BAUM_WELCH_TRAINING_INL_H_

#include "baum_welch_training.h"

#include <cmath>
#include <cstdio>

#include <numeric>
#include <string>

#include "context_profile-inl.h"
#include "forward_backward_algorithm.h"
#include "hmm-inl.h"
#include "log.h"
#include "progress_table.h"
#include "shared_ptr.h"
#include "substitution_matrix-inl.h"
#include "utils-inl.h"

namespace cs {

template< class Alphabet, template<class> class Subject >
BaumWelchTraining<Alphabet, Subject>::BaumWelchTraining(
    const BaumWelchOptions& opts,
    const DataVec& data,
    const SubstitutionMatrix<Alphabet>* sm,
    HMM<Alphabet>& hmm)
    : ExpectationMaximization<Alphabet, Subject>(data),
      opts_(opts),
      sm_(sm),
      hmm_(hmm),
      emission_(hmm.num_cols(), opts.weight_center, opts.weight_decay),
      transition_stats_(hmm.num_states(), hmm.num_states()),
      prior_stats_(0.0, hmm.num_states()),
      transition_stats_block_(hmm.num_states(), hmm.num_states()),
      prior_stats_block_(0.0, hmm.num_states()) {
  Init();
}

template< class Alphabet, template<class> class Subject >
BaumWelchTraining<Alphabet, Subject>::BaumWelchTraining(
    const BaumWelchOptions& opts,
    const DataVec& data,
    const SubstitutionMatrix<Alphabet>* sm,
    HMM<Alphabet>& hmm,
    FILE* fout)
    : ExpectationMaximization<Alphabet, Subject>(data),
      opts_(opts),
      sm_(sm),
      hmm_(hmm),
      emission_(hmm.num_cols(), opts.weight_center, opts.weight_decay),
      transition_stats_(hmm.num_states(), hmm.num_states()),
      prior_stats_(0.0, hmm.num_states()),
      transition_stats_block_(hmm.num_states(), hmm.num_states()),
      prior_stats_block_(0.0, hmm.num_states()) {
  progress_table_.reset(new BaumWelchProgressTable<Alphabet, Subject>(this, fout));
  Init();
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::ExpectationStep(
    const DataVec& block) {
  LOG(INFO) << "Starting E-step ...";

  // Run forward and backward algorithm on each subject in current block
  const int block_size = block.size();
#pragma omp parallel for schedule(static)
  for (int n = 0; n < block_size; ++n) {
    ForwardBackwardMatrices fbm(block[n]->length(), hmm_.num_states());
    ForwardBackwardAlgorithm(hmm_, *block[n], emission_, &fbm);

    AddContributionToTransitions(fbm);
    AddContributionToStates(fbm, *block[n]);

#pragma omp atomic
    log_likelihood_ += fbm.log_likelihood / num_eff_cols_;
    if (progress_table_) {
#pragma omp critical (print_progress)
      progress_table_->print_progress(hmm_.num_states() * block[n]->length());
    }
  }

  UpdateSufficientStatistics();
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::AddContributionToTransitions(
    const ForwardBackwardMatrices& m) {
  const int slen = m.f.num_rows();

  for (ConstTransitionIter ti = hmm_.transitions_begin();
       ti != hmm_.transitions_end(); ++ti) {
    double w_kl = 0.0;
    for (int i = 0; i < slen-1; ++i) {
      w_kl += m.f[i][ti->source] * m.b[i+1][ti->target] * ti->weight *
        m.e[i+1][ti->target] / m.s[i+1];
    }

#pragma omp critical (add_contribution_to_transition)
    if (!transition_stats_block_.test(ti->source, ti->target)) {
      transition_stats_block_[ti->source][ti->target] = w_kl;
    } else {
      transition_stats_block_[ti->source][ti->target] = w_kl +
        transition_stats_block_[ti->source][ti->target];
    }
  }
}

template< class Alphabet, template<class> class Subject >
inline void BaumWelchTraining<Alphabet, Subject>::AddContributionToStates(
    const ForwardBackwardMatrices& m,
    const CountProfile<Alphabet>& c) {
  const int slen          = c.length();
  const int num_states    = hmm_.num_states();
  const int alphabet_size = hmm_.alphabet_size();
  const int ci            = (hmm_.num_cols() - 1) / 2;

  for (int k = 0; k < num_states; ++k) {
#pragma omp atomic
    prior_stats_block_[k] += m.f[0][k] * m.b[0][k];

    ProfileStats& p_k = *profile_stats_block_[k];
    for (int i = 0; i < slen; ++i) {
      const int beg = std::max(0, i - ci);
      const int end = std::min(c.num_cols() - 1, i + ci);

      for(int h = beg; h <= end; ++h) {
        const int j = h - i + ci;
        for (int a = 0; a < alphabet_size; ++a) {
#pragma omp atomic
          p_k[j][a] += c.counts(h,a) * m.f[i][k] * m.b[i][k];
        }
      }
    }
  }
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::AddContributionToStates(
    const ForwardBackwardMatrices& m,
    const Sequence<Alphabet>& s) {
  const int slen       = s.length();
  const int num_states    = hmm_.num_states();
  const int ci            = (hmm_.num_cols() - 1) / 2;

  for (int k = 0; k < num_states; ++k) {
#pragma omp atomic
    prior_stats_block_[k] += m.f[0][k] * m.b[0][k];

    ProfileStats& p_k = *profile_stats_block_[k];
    for (int i = 0; i < slen; ++i) {
      const int beg = std::max(0, i - ci);
      const int end = std::min(s.length() - 1, i + ci);

      for(int h = beg; h <= end; ++h) {
        const int j = h - i + ci;
#pragma omp atomic
        p_k[j][s[h]] += m.f[i][k] * m.b[i][k];
      }
    }
  }
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::UpdateSufficientStatistics() {
  const double gamma      = 1.0 - epsilon_;
  const int num_states    = hmm_.num_states();
  const int num_cols      = hmm_.num_cols();
  const int alphabet_size = hmm_.alphabet_size();

  // Update transition statistics
  for (ConstTransitionIter ti = hmm_.transitions_begin();
       ti != hmm_.transitions_end(); ++ti) {
    if (transition_stats_block_.test(ti->source, ti->target)) {
      if (!transition_stats_.test(ti->source, ti->target))
        transition_stats_[ti->source][ti->target] = 0.0f;
      transition_stats_[ti->source][ti->target] =
        gamma * transition_stats_[ti->source][ti->target] +
        transition_stats_block_[ti->source][ti->target];
    }
  }

  // Update priors and emissions statistics
  #pragma omp parallel for schedule(static)
  for (int k = 0; k < num_states; ++k) {
    prior_stats_[k] = prior_stats_[k] * gamma + prior_stats_block_[k];

    ProfileStats& p_block = *profile_stats_block_[k];
    ProfileStats& p       = *profile_stats_[k];
    for (int j = 0; j < num_cols; ++j) {
      for (int a = 0; a < alphabet_size; ++a) {
        p[j][a] = gamma * p[j][a] + p_block[j][a];
          }
    }
  }
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::MaximizationStep() {
  const int num_states    = hmm_.num_states();
  const int num_cols      = hmm_.num_cols();
  const int alphabet_size = hmm_.alphabet_size();

  // Normalize prior state probabilities
  prior_stats_ /= prior_stats_.sum();

  // Assign new priors and emission probabilities
#pragma omp parallel for schedule(static)
  for (int k = 0; k < num_states; ++k) {
    // Update prior probability
    hmm_[k].set_prior(prior_stats_[k]);

    // Update emissions
    ProfileStats& p_k = *profile_stats_[k];
    for (int j = 0; j < num_cols; ++j) {
      double sum = std::accumulate(&p_k[j][0], &p_k[j][0] + alphabet_size, 0.0);

      if (sum != 0.0) {
        double norm_fac = 1.0 / sum;
        for (int a = 0; a < alphabet_size; ++a)
          hmm_[k][j][a] = log2(p_k[j][a] * norm_fac);
      }
    }
  }

  // Calculate and assign new transition probabilities
  hmm_.ClearTransitions();
  for (int k = 0; k < num_states; ++k) {
    double tr_sum = 0.0;

    for (int l = 0; l < num_states; ++l) {
      if (transition_stats_.test(k,l)) {
        double a_kl = transition_stats_[k][l];

        if (opts_.transition_pc != 1.0f &&
            hmm_.connectivity() > opts_.max_connectivity)
          a_kl += opts_.transition_pc - 1.0f;

        if (a_kl > 0.0) {
          transition_stats_[k][l] = a_kl;
          tr_sum += a_kl;
        } else {
          transition_stats_.erase(k,l);
        }
      }
    }

    if (tr_sum != 0.0) {
      double tr_fac = 1.0 / tr_sum;

      for (int l = 0; l < num_states; ++l) {
        if (transition_stats_.test(k,l)) {
          hmm_(k,l) = static_cast<float>(transition_stats_[k][l] * tr_fac);
          LOG(DEBUG2) << strprintf("tr[%i][%i]=%-15.10f",
                                   k, l, static_cast<float>(hmm_(k,l)));
        }
      }
    }
  }

  hmm_.increment_iterations();
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::Init() {
  const int alphabet_size = hmm_.alphabet_size();

  // Create profiles for global and block-level sufficient statistics
  for (int k = 0; k < hmm_.num_states(); ++k) {
    ProfileStatsPtr p(new ProfileStats(hmm_.num_cols(), alphabet_size, 0.0));
    profile_stats_.push_back(p);

    ProfileStatsPtr pb(new ProfileStats(hmm_.num_cols(), alphabet_size, 0.0));
    profile_stats_block_.push_back(pb);
  }

  // Compute total number of data columns
  long num_cols = 0;
  for (typename DataVec::const_iterator di = data_.begin();
       di != data_.end(); ++di)
    num_cols += (**di).length();

  if (progress_table_)
    progress_table_->set_total_work(num_cols * hmm_.num_states());

  // Set number of effective columsn for log-likelihood calculation
  num_eff_cols_ = emission_.SumWeights() * num_cols;
}

template< class Alphabet, template<class> class Subject >
void BaumWelchTraining<Alphabet, Subject>::ResetAndAddPseudocounts() {
  const int num_cols      = hmm_.num_cols();
  const int alphabet_size = hmm_.alphabet_size();

  // Erase all transitions. Pseudocounts are taken care of in M-step.
  transition_stats_block_.clear();

  // Reset profile and prior evidences to their pseudocount values for numeric
  // stability.
  for (int k = 0; k < hmm_.num_states(); ++k) {
    prior_stats_block_[k] = opts_.prior_pc;

    ProfileStats& p = *profile_stats_block_[k];
    for (int j = 0; j < num_cols; ++j) {
      for (int a = 0; a < alphabet_size; ++a) {
        p[j][a] = opts_.profile_pc * sm_->f(a);
      }
    }
  }
}

template< class Alphabet, template<class> class Subject >
inline bool BaumWelchTraining<Alphabet, Subject>::IsDone() const {
  if (scan_ < opts_.min_scans)
    return false;
  else if (scan_ >= opts_.max_scans)
    return true;
  else if (opts_.max_connectivity == 0)
    return fabs(log_likelihood_change()) <= opts_.log_likelihood_change;
  else
    return (fabs(log_likelihood_change()) <= opts_.log_likelihood_change &&
            hmm_.connectivity() <= opts_.max_connectivity);
}


template< class Alphabet, template<class> class Subject >
BaumWelchProgressTable<Alphabet, Subject>::BaumWelchProgressTable(
    const BaumWelchTraining<Alphabet, Subject>* training,
    FILE* fout,
    int width)
    : ProgressTable(fout, width),
      training_(training) {}

template< class Alphabet, template<class> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::PrintHeader() {
  fprintf(fout_, "%-4s %4s %6s %4s %7s  %-30s  %9s  %8s\n",
          "Scan", "Itrs", "Conn", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
  fputs(std::string(82, '-').c_str(), fout_);
  fputc('\n', fout_);
}

template< class Alphabet, template<class> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::PrintRowBegin() {
  Reset();
  fprintf(fout_, "%-4i %4i %6.1f %4i %7.4f  ", training_->scan(),
          training_->iterations(), training_->hmm_.connectivity(),
          training_->num_blocks(), training_->epsilon());
  fflush(fout_);
}

template< class Alphabet, template<class> class Subject >
void BaumWelchProgressTable<Alphabet, Subject>::PrintRowEnd() {
  if (training_->scan() == 1)
    fprintf(fout_, "  %9.5f\n", training_->log_likelihood());
  else
    fprintf(fout_, "  %9.5f  %+8.5f\n", training_->log_likelihood(),
            training_->log_likelihood_change());
  fflush(fout_);
}

}  // namespace cs

#endif  // SRC_BAUM_WELCH_TRAINING_INL_H_
