// Copyright 2009, Andreas Biegert

#ifndef SRC_CLUSTERING_INL_H_
#define SRC_CLUSTERING_INL_H_

#include "clustering.h"

#include <cmath>
#include <cstdio>

#include <string>

#include "context_profile-inl.h"
#include "mult_emission-inl.h"
#include "profile_library-inl.h"
#include "log.h"
#include "progress_table.h"
#include "shared_ptr.h"
#include "utils-inl.h"

namespace cs {

template< class Alphabet,
          template<class A> class Subject >
Clustering<Alphabet, Subject>::Clustering(const ClusteringOptions& opts,
                                          const data_vector& data,
                                          ProfileLibrary<Alphabet>& lib)
    : ExpectationMaximization<Alphabet, Subject>(data),
      opts_(opts),
      lib_(lib),
      emission_(lib.num_cols(), opts.weight_center, opts.weight_decay),
      profile_stats_(),
      profile_stats_block_() {
  init();
}

template< class Alphabet,
          template<class A> class Subject >
Clustering<Alphabet, Subject>::Clustering(const ClusteringOptions& opts,
                                          const data_vector& data,
                                          ProfileLibrary<Alphabet>& lib,
                                          FILE* fout)
    : ExpectationMaximization<Alphabet, Subject>(data),
      opts_(opts),
      lib_(lib),
      emission_(lib.num_cols(), opts.weight_center, opts.weight_decay),
      profile_stats_(),
      profile_stats_block_() {
  progress_table_ = new ClusteringProgressTable<Alphabet, Subject>(this, fout);
  init();
}

template< class Alphabet,
          template<class A> class Subject >
Clustering<Alphabet, Subject>::~Clustering() {
  if (progress_table_) delete progress_table_;
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::expectation_step(const data_vector& block) {
  const int num_profiles = lib_.num_profiles();
  std::valarray<double> p_zn(0.0f, lib_.num_profiles());

  // Compute posterior probabilities p_zn[k] of profile k for counts n
  for (typename data_vector::const_iterator bi = block.begin();
       bi != block.end(); ++bi) {
    double sum = 0.0f;
    for (int k = 0; k < num_profiles; ++k) {
      p_zn[k] = lib_[k].prior() *
        fast_pow2(emission_(lib_[k], **bi, lib_[k].center()));
      sum += p_zn[k];

      LOG(DEBUG2) << strprintf("a(%i)=%-8.5g P(c_n|p_%i)=%-8.5g P(z_n=%-4i)=%-8.5g",
                               k, lib_[k].prior(), k,
                               emission_(lib_[k], **bi, lib_[k].center()), k, p_zn[k]);
    }
    p_zn /= sum;
    add_contribution_to_priors(p_zn);
    add_contribution_to_emissions(p_zn, **bi);
    log_likelihood_ += log2(sum) / num_eff_cols_;
    LOG(DEBUG1) << strprintf("log(L)=%-8.5g", log_likelihood_);

    if (progress_table_)
      progress_table_->print_progress(lib_.num_profiles());
  }

  update_sufficient_statistics();
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::maximization_step() {
  const int num_profiles  = lib_.num_profiles();
  const int num_cols      = lib_.num_cols();
  const int alphabet_size = lib_.alphabet_size();

  float sum = 0.0f;
  for (int k = 0; k < num_profiles; ++k) sum += profile_stats_[k]->prior();
  float fac = 1.0f / sum;

  // Assign new priors and emission probabilities
  for (int k = 0; k < num_profiles; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_[k];
    LOG(DEBUG1) << p_k;

    lib_[k].set_prior(p_k.prior() * fac);
    ContextProfile<Alphabet> tmp(p_k);
    if (normalize(&tmp)) {  // don't update profiles that did'n get evidence
      tmp.transform_to_logspace();
      for (int i = 0; i < num_cols; ++i)
        for (int a = 0; a < alphabet_size; ++a)
          lib_[k][i][a] = tmp[i][a];
    }
  }

  ++lib_;  // increment iteration counter
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::add_contribution_to_priors(
    const std::valarray<double>& p_zn) {
  const int num_profiles = lib_.num_profiles();

  for (int k = 0; k < num_profiles; ++k)
    profile_stats_block_[k]->set_prior(profile_stats_block_[k]->prior() + p_zn[k]);
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::add_contribution_to_emissions(
    const std::valarray<double>& p_zn,
    const CountProfile<Alphabet>& c) {
  const int num_profiles  = lib_.num_profiles();
  const int num_cols      = lib_.num_cols();
  const int alphabet_size = lib_.alphabet_size();

  for (int k = 0; k < num_profiles; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_block_[k];

    for (int j = 0; j < num_cols; ++j) {
      for (int a = 0; a < alphabet_size; ++a) {
        p_k[j][a] += c[j][a] * p_zn[k];
      }
    }
  }
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::add_contribution_to_emissions(
    const std::valarray<double>& p_zn,
    const Sequence<Alphabet>& s) {
  const int num_profiles  = lib_.num_profiles();
  const int num_cols      = lib_.num_cols();
  const int alphabet_size = lib_.alphabet_size();

  for (int k = 0; k < num_profiles; ++k) {
    ContextProfile<Alphabet>& p_k = *profile_stats_block_[k];

    for (int j = 0; j < num_cols; ++j) {
      p_k[j][s[j]] += p_zn[k];
    }
  }
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::update_sufficient_statistics() {
  const float gamma       = 1.0f - epsilon_;
  const int num_profiles  = lib_.num_profiles();
  const int num_cols      = lib_.num_cols();
  const int alphabet_size = lib_.alphabet_size();

  // Update priors and emissions statistics
  for (int k = 0; k < num_profiles; ++k) {
    ContextProfile<Alphabet>& p_block = *profile_stats_block_[k];
    ContextProfile<Alphabet>& p       = *profile_stats_[k];

    p.set_prior(p.prior() * gamma + p_block.prior());
    for (int j = 0; j < num_cols; ++j) {
      for (int a = 0; a < alphabet_size; ++a) {
        p[j][a] = gamma * p[j][a] + p_block[j][a];
      }
    }
    reset(&p_block);
  }
}

template< class Alphabet,
          template<class A> class Subject >
void Clustering<Alphabet, Subject>::init() {
  // Create profiles for global and block-level sufficient statistics
  for (int k = 0; k < lib_.num_profiles(); ++k) {
    profile_stats_.push_back(
        shared_ptr< ContextProfile<Alphabet> >(
            new ContextProfile<Alphabet>(k, lib_.num_cols())));
    profile_stats_block_.push_back(
        shared_ptr< ContextProfile<Alphabet> >(
            new ContextProfile<Alphabet>(k, lib_.num_cols())));
  }

  // Initialize total amount of work per scan
  if (progress_table_)
    progress_table_->set_total_work(static_cast<long>(lib_.num_profiles()) *
                                    data_.size());

  // Compute effective number of training data columns
  num_eff_cols_ = emission_.SumWeights() * data_.size();
}



template< class Alphabet,
          template<class A> class Subject >
ClusteringProgressTable<Alphabet, Subject>::ClusteringProgressTable(
    const Clustering<Alphabet, Subject>* clustering,
    FILE* fout,
    int width)
    : ProgressTable(fout, width),
      clustering_(clustering) {}

template< class Alphabet,
          template<class A> class Subject >
void ClusteringProgressTable<Alphabet, Subject>::print_header() {
  fprintf(fout_, "%-4s %4s %4s %7s  %-30s  %9s  %8s\n",
          "Scan", "Itrs", "Blks", "Epsilon", "E-Step", "log(L)", "+/-");
  fputs(std::string(75, '-').c_str(), fout_);
  fputc('\n', fout_);
}

template< class Alphabet,
          template<class A> class Subject >
void ClusteringProgressTable<Alphabet, Subject>::print_row_begin() {
  reset();
  fprintf(fout_, "%-4i %4i %4i %7.4f  ", clustering_->scan(),
          clustering_->iterations(), clustering_->num_blocks(),
          clustering_->epsilon());
  fflush(fout_);
}

template< class Alphabet,
          template<class A> class Subject >
void ClusteringProgressTable<Alphabet, Subject>::print_row_end() {
  if (clustering_->scan() == 1)
    fprintf(fout_, "  %9.5f\n", clustering_->log_likelihood());
  else
    fprintf(fout_, "  %9.5f  %+8.5f\n", clustering_->log_likelihood(),
            clustering_->log_likelihood_change());
  fflush(fout_);
}

}  // namespace cs

#endif  // SRC_CLUSTERING_INL_H_
