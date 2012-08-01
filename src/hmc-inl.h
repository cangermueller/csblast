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

#ifndef CS_HMC_INL_H_
#define CS_HMC_INL_H_

#include "hmc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "globals.h"
#include "crf-inl.h"
#include "log.h"
#include "sequence-inl.h"
#include "shared_ptr.h"
#include "substitution_matrix-inl.h"
#include "utils.h"

namespace cs {

template< class Alphabet, template<class> class Subject >
Hmc<Alphabet, Subject>::Hmc(const std::vector<DataPtr>& data,
                            const std::vector<PredPtr>& pred,
                            const HmcOptions* opts,
                            const SubstitutionMatrix<Alphabet>* sm,
                            Crf<Alphabet>* crf,
                            FILE* fout)
    : alphabet_size_(Alphabet::instance().size()),
      num_cols_(crf->num_cols()),
      num_states_(crf->num_states()),
      num_data_(data.size()),
      num_weights_(num_states_ * (1 + (num_cols_ + 1) * alphabet_size_)),
      center_((crf->num_cols() - 1) / 2),
      data_(data),
      pred_(pred),
      opts_(opts),
      sm_(sm),
      crf_(*crf),
      progress_table_(new HmcProgressTable<Alphabet, Subject>(this, fout)),
      grad1_(0.0, num_weights_),
      grad2_(0.0, num_weights_),
      mom_(0.0, num_weights_),
      best_(NULL),
      memento_(*crf, num_weights_),
      logp_(0.0),
      logr_(0.0),
      hamiltonian_end_(0.0),
      hamiltonian_beg_(0.0),
      iter_(0),
      num_scans_(0),
      integration_time_(0.0f),
      num_accepted_(0),
      epsilon_(opts->epsilon),
      temperature_(opts->temperature),
      num_replicas_(1),
      myrank_(0),
      rng_local_(NULL),
      rng_replica_exch_(NULL),
      rng_participants_(NULL),
      rng_swap_(NULL) {
  Init();
}

template< class Alphabet, template<class> class Subject >
Hmc<Alphabet, Subject>::~Hmc() {
  gsl_rng_free(rng_local_);
  gsl_rng_free(rng_replica_exch_);
  gsl_rng_free(rng_participants_);
  gsl_rng_free(rng_swap_);
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::Run() {
  assert(opts_->leapfrog_steps % (2 * opts_->num_batches) == 0);

  // Calculate initial likelihood and gradient of first batch
  for (RevBatchIter bi = batches_.rbegin(); bi != batches_.rend(); ++bi) {
    CalculateLikelihoodAndGradient(*bi, false);
  }

  // Store initial solution
  SamplePtr sample(new Sample<Alphabet>(crf_, iter_, logp_, true));
  samples_.push_back(sample);
  best_ = sample;
  if (progress_table_)
    progress_table_->PrintHeader();

  // Perform HMC sampling iterations
  for (iter_ = 1; iter_ <= opts_->sampling_iters; ++iter_) {
    // Store current state in case we don't accept the next sample
    SaveToMemento();
    // Print table row up to progress bar
    if (progress_table_)
      progress_table_->PrintRowBegin();
    // Sample momentum terms from isotropic gaussian with stddev 1
    for (int i = 0; i < num_weights_; ++i)
      mom_[i] = gsl_ran_gaussian(rng_local_, 1.0);
    // Save hamiltonian at beginning of leapfrog
    hamiltonian_beg_ = GetHamiltonian();
    // Generate next parameter sample by running the leapfrog integration
    RunLeapfrogIntegration();
    // Save hamiltonian at the end of leapfrog
    hamiltonian_end_ = GetHamiltonian();
    // Evaluate Metropolis criterion based on total energies
    double p_acc = MetropolisProb(hamiltonian_beg_, hamiltonian_end_);
    bool accepted = (gsl_ran_flat(rng_local_, 0.0, 1.0) < p_acc);
    if (accepted) ++num_accepted_;
    // Store sample in vector and update pointer to best solution
    sample = SamplePtr(new Sample<Alphabet>(crf_, iter_, logp_, accepted));
    samples_.push_back(sample);
    // Report this sampling step
    if (progress_table_) progress_table_->PrintRowEnd();
    // Keep track of the best sample found so far
    if(logp_ > best_->logp) best_ = sample;
    // Accept or discard sampled solution based on Metropolis decision
    if (accepted) {
      epsilon_ *= opts_->rescale_accept;
    } else {
      epsilon_ *= opts_->rescale_reject;
      RestoreFromMemento();
    }
#ifdef PARALLEL
    // Should we try a replica exchange?
    double r = gsl_ran_flat(rng_replica_exch_, 0.0, 1.0);
    MPI_Barrier(MPI_COMM_WORLD);
    if (r < opts_->theta) {
      ReplicaExchange();
    }
#endif
  }

  // Replace CRF with best sample
  crf_ = best_->crf;
  //crf_ = samples_.back()->crf;
  //crf_ = samples_.front()->crf;
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::RunLeapfrogIntegration() {
  // We keep a running average of CRF weights and momentas to reduce oscillation.
  // const int kAverLen   = 50;
  // const float kUpscale = 1.5f;
  // const float eps_bak  = epsilon_;
  // std::valarray<double> wgt_sum(0.0, crf_.num_weights());  // sum of CRF weights
  // std::valarray<double> mom_sum(0.0, crf_.num_weights());  // sum of momenta
  // int num_aver = 0;

#ifndef NDEBUG
  const float mom_length  = sqrt((mom_ * mom_).sum());
  const float grad_length = sqrt(((grad1_ + grad2_) * (grad1_ + grad2_)).sum());
  const std::valarray<double> mom(mom_);
#endif

  // Start with a half step
  mom_ += 0.5 * epsilon_ * ((grad1_ / temperature_) + grad2_);
  // mom_ += 0.5 * epsilon_ * grad2_;
  // Run remaining leapfrog steps
  for (int n = 0; n < opts_->leapfrog_steps; ++n) {
    int m = n % (2 * opts_->num_batches);
    int b = m < opts_->num_batches ? m : 2 * opts_->num_batches - m - 1;

    // Update CRF weights based on momenta
    for (int k = 0, i = 0; k < num_states_; ++k) {
      CrfState<Alphabet>& state = crf_[k];
      // The order in which we update the weight of each state is important:
      // first bias weight, then column weights, and finally pseudocount weights
      state.bias_weight() += epsilon_ * mom_[i++];
      for (int j = 0; j < num_cols_; ++j)
        for (int a = 0; a < alphabet_size_; ++a)
          state[j][a] += epsilon_ * mom_[i++];
      for (int a = 0; a < alphabet_size_; ++a)
        state(a) += epsilon_ * mom_[i++];
      // Precompute pseudocount emissions and columns weights
      state.UpdatePseudocounts();
    }

    // Update running sum of CRF weights and momentum terms
    // std::valarray<double> wgt(0.0, crf_.num_weights());
    // Pack(crf_, &wgt[0]);
    // wgt_sum += wgt;
    // mom_sum += mom_;

    // Replace CRF weights and momenta with average over last 20 iterations.
    // This should allow us to increase epsilon by a good factor.
    // if ((n + 1) % kAverLen == 0 && num_aver < 6) {
    //   wgt_sum /= static_cast<double>(kAverLen);
    //   Unpack(&wgt_sum[0], &crf_);
    //   mom_ = mom_sum / static_cast<double>(kAverLen);

    //   epsilon_ *= kUpscale;
    //   wgt_sum = 0.0;
    //   mom_sum = 0.0;
    //   ++num_aver;
    // }

    // Reset likelihood at start of new scan
    if (n % opts_->num_batches == 0) {
      logp_ = 0.0;
      logr_ = 0.0;
      ++num_scans_;
    }
    // Calculate new gradient based on batch b
    CalculateLikelihoodAndGradient(batches_[b]);
    // Update momenta based on new
    const double step_len = (n + 1 == opts_->leapfrog_steps) ? 0.5 : 1.0;
    mom_ += step_len * epsilon_ * ((grad1_ / temperature_) + grad2_);
    //mom_ += step_len * epsilon_ * grad2_;

    integration_time_ += epsilon_;

#ifndef NDEBUG
    // Calculate and report decoherence and gradient norm
    float decoh = (mom * mom_).sum() / (mom_length * sqrt((mom_ * mom_).sum()));
    float grad_length_curr = sqrt(((grad1_ + grad2_) * (grad1_ + grad2_)).sum());
    //LOG(ERROR) << strprintf("iter= %-4d n= %-4d grad= %7.5f  decoh= %7.5f", iter_,
    //                        n+1, grad_length_curr / grad_length, decoh);
#endif
  }

  //epsilon_ = eps_bak;
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::CalculateLikelihoodAndGradient(const Batch& batch,
                                                            bool verbose) {
  const int batch_size   = batch.size();
  const float batch_frac = static_cast<float>(batch_size) / num_data_;

  // Allocate helper matrices
  double** mpp = new double*[batch_size];
  double** mpa = new double*[batch_size];
  for (int b = 0; b < batch_size; ++b) {
    mpp[b] = new double[num_states_];
    mpa[b] = new double[alphabet_size_];
    for (int a = 0; a < alphabet_size_; ++a) mpa[b][a] = 0.0;
  }

  double logp_batch = 0.0;  // log likelihood contribution of the whole batch
#pragma omp parallel for schedule(static)
  for (int b = 0; b < batch_size; ++b) {
    const int n = batch[b];
    double* pp = &mpp[b][0];
    double* pa = &mpa[b][0];

    // Calculate posterior probability pp[k] of state k given count profile n
    double pp_sum = 0.0;
    for (int k = 0; k < num_states_; ++k) {
      pp[k] = MatchProb(crf_[k], *data_[n]);
      pp_sum += pp[k];
    }
    // Normalize posteriors and calculate pseudocounts given count profile n
    double norm_fac = 1.0 / pp_sum;
    for (int k = 0; k < num_states_; ++k) {
      pp[k] *= norm_fac;
      for (int a = 0; a < alphabet_size_; ++a)
        pa[a] += crf_[k].pc(a) * pp[k];
    }
    // Update log-likelihood based on training profile n
    double logp_n = 0.0;  // log-likelihood contribution of training profile n
    for (int a = 0; a < alphabet_size_; ++a)
      logp_n += (*pred_[n])[center_][a] * (log(pa[a]) - log(sm_->f(a)));
#pragma omp atomic
    logp_batch += logp_n;
  }

  // Add derivative terms for log likelihood to gradient
  CalculateGradient(batch, mpp, mpa, verbose);

  // Deallocate helper matrices
  for (int b = 0; b < batch_size; ++b) {
    delete[] mpp[b];
    delete[] mpa[b];
  }
  delete[] mpp;
  delete[] mpa;

  // Update log likelihood and regularization member variables
  logp_ += logp_batch;
  logr_ += batch_frac * GetRegularizer();
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::CalculateGradient(const Batch& batch,
                                               double** mpp,
                                               double** mpa,
                                               bool verbose) {
  // Parallelization over the training patterns instead of CRF states allows for
  // better speedup because we don't need to synchronize write access to gradient
  // vector
  const int batch_size = batch.size();
  const int letter_any = Alphabet::instance().any();
  double* g = &grad1_[0];
  grad1_ = 0.0;

#pragma omp parallel for schedule(static)
  for (int k = 0; k < num_states_; ++k) {
    const CrfState<Alphabet>& state = crf_[k];
    const value_type* pc = &state.pc(0);
    double fit, sum = 0.0;

    for (int b = 0; b < batch_size; ++b) {
      const Sequence<Alphabet>& seq = *data_[batch[b]];
      int i = k * (1 + (num_cols_ + 1) * alphabet_size_);
      double* pcp = &mpa[b][0];                       // predicted pseudocounts
      float*  pcl = &pred_[batch[b]]->at(center_,0);  // label pseudocounts

      // Precumpute fit
      fit = 0.0;
      for (int a = 0; a < alphabet_size_; ++a)
        fit += pcl[a] * (pc[a] / pcp[a] - 1.0);
      // Update gradient of bias weight
      g[i++] += mpp[b][k] * fit;
      // Update gradient terms of context weights
      for(int j = 0; j < num_cols_; ++j)
        if (seq[j] != letter_any)
          g[i + j*alphabet_size_ + seq[j]] += mpp[b][k] * fit;
      i += num_cols_ * alphabet_size_;
      // Precompute sum needed for gradient update of pseudocounts weights
      sum = 0.0;
      for (int a = 0; a < alphabet_size_; ++a)
        sum += pc[a] * pcl[a] / pcp[a];
      // Update gradient terms of pseudocount weights
      for (int a = 0; a < alphabet_size_; ++a)
        g[i+a] += mpp[b][k] * pc[a] * (pcl[a] / pcp[a] - sum);
    }
    // Increment progress bar
    if (verbose && progress_table_) {
#pragma omp critical (print_progress)
      progress_table_->print_progress(1);
    }
  }

  // Now that likelihood terms have been included we add derivate terms of
  // regularizers to the gradient. This is fast, so no parallelization is needed.
  double batch_frac  = static_cast<double>(batch_size) / num_data_;
  double fac_bias = -batch_frac / (opts_->sigma_bias * opts_->sigma_bias);
  double fac_pc   = -batch_frac / (opts_->sigma_pc * opts_->sigma_pc);
  grad2_ = 0.0;
  for (int k = 0, i = 0; k < num_states_; ++k) {
    CrfState<Alphabet>& state = crf_[k];
    grad2_[i++] += fac_bias * state.bias_weight();
    for(int j = 0; j < num_cols_; ++j) {
      double s = opts_->sigma_context * pow(opts_->sigma_decay, fabs(j - center_));
      double fac_context = -batch_frac / (s * s);
      for (int a = 0; a < alphabet_size_; ++a)
        grad2_[i++] += fac_context * state[j][a];
    }
    for (int a = 0; a < alphabet_size_; ++a)
      grad2_[i++] += fac_pc * state.pc_weight(a);
  }

  int beg = 1 + center_ * alphabet_size_;
  int end = 1 + (center_ + 1) * alphabet_size_;
  LOG(ERROR) << "------------------------------";
  for (int i = beg; i < beg + 1; ++i){
    LOG(ERROR) << strprintf("grad1_cw[%d]=%+8.4f", i, grad1_[i]);
    LOG(ERROR) << strprintf("grad2_cw[%d]=%+8.4f", i, grad2_[i]);
    LOG(ERROR) << strprintf("  mom_cw[%d]=%+8.4f", i, mom_[i]);
  }

  beg = 1 + num_cols_ * alphabet_size_;
  end = 1 + (num_cols_ + 1) * alphabet_size_;
  LOG(ERROR) << "------------------------------";
  for (int i = beg; i < beg + 1; ++i){
    LOG(ERROR) << strprintf("grad1_pc[%d]=%+8.4f", i, grad1_[i]);
    LOG(ERROR) << strprintf("grad2_pc[%d]=%+8.4f", i, grad2_[i]);
    LOG(ERROR) << strprintf("  mom_pc[%d]=%+8.4f", i, mom_[i]);
  }
}

template< class Alphabet, template<class> class Subject >
double Hmc<Alphabet, Subject>::GetRegularizer() {
  double fac_bias = -0.5 / (opts_->sigma_bias * opts_->sigma_bias);
  double fac_pc   = -0.5 / (opts_->sigma_pc * opts_->sigma_pc);
  double rv = 0.0;

  for (int k = 0; k < num_states_; ++k) {
    CrfState<Alphabet>& state = crf_[k];
    rv += fac_bias * state.bias_weight() * state.bias_weight();
    for(int j = 0; j < num_cols_; ++j) {
      double s = opts_->sigma_context * pow(opts_->sigma_decay, fabs(j - center_));
      double fac_context = -0.5 / (s * s);
      for (int a = 0; a < alphabet_size_; ++a)
        rv += fac_context * state[j][a] * state[j][a];
    }
    for (int a = 0; a < alphabet_size_; ++a)
      rv += fac_pc * state(a) * state(a);
  }
  return rv;
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::ReplicaExchange() {
  if (num_replicas_ <= 1) return;

#ifdef PARALLEL
  double r = gsl_ran_flat(rng_participants_, 0.0, num_replicas_ - 2);
  const int l = static_cast<int>(floor(r));  // index of lower temperature level
  const int lower = replicas_[l];            // rank of replica with lower temp
  const int upper = replicas_[l+1];          // rank of replica with upper temp
  double temperature = 0;
  double logp        = 0;
  MPI_Status status;
  MPI_Comm comm_world_swap, comm_world_bcast;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_swap);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm_world_bcast);

  // Draw random number for swapping already now to kepp all rng's in sync
  r = gsl_ran_flat(rng_swap_, 0.0, 1.0);

  if (r < opts_->theta) {
    if (myrank_ == lower) {
      // Swap temperature and likelihood with upper replica
      MPI_Sendrecv(&temperature_, 1, MPI_DOUBLE, upper, kTagTemperatureSwap,
                   &temperature , 1, MPI_DOUBLE, upper, kTagTemperatureSwap,
                   comm_world_swap, &status);
      MPI_Sendrecv(&logp_, 1, MPI_DOUBLE, upper, kTagLikelihoodSwap,
                   &logp , 1, MPI_DOUBLE, upper, kTagLikelihoodSwap,
                   comm_world_swap, &status);

      // Evaluate replicae exchange criterion
      double p_swap = ReplicaExchangeProb(temperature_, logp_, temperature, logp);
      if (r < p_swap) {
        LOG(INFO) << strprintf("%2d: i=%d  LL=%10.2f T=%5.2f => LL=%10.2f T=%5.2f",
                               myrank_, iter_, logp_, temperature_, logp,
                               temperature);
        // Swap temperatures
        temperature_ = temperature;
        // Switch replica positions in replica table
        replicas_[l]   = upper;
        replicas_[l+1] = lower;
      }

    } else if (myrank_ == upper) {
      // Swap temperature and likelihood with lower replica
      MPI_Sendrecv(&temperature_, 1, MPI_DOUBLE, lower, kTagTemperatureSwap,
                   &temperature , 1, MPI_DOUBLE, lower, kTagTemperatureSwap,
                   comm_world_swap, &status);
      MPI_Sendrecv(&logp_, 1, MPI_DOUBLE, lower, kTagLikelihoodSwap,
                   &logp , 1, MPI_DOUBLE, lower, kTagLikelihoodSwap,
                   comm_world_swap, &status);

      // Evaluate replicae exchange criterion
      double p_swap = ReplicaExchangeProb(temperature, logp, temperature_, logp_);
      if (r < p_swap) {
        LOG(INFO) << strprintf("%2d: i=%d  LL=%10.2f T=%5.2f <= LL=%10.2f T=%5.2f",
                               myrank_, iter_, logp, temperature, logp_,
                               temperature_);
        // Swap temperatures
        temperature_ = temperature;
      }
    }

    // Lower replica broadcasts new replica order to all other replicas
    MPI_Bcast (&replicas_[0], num_replicas_, MPI_INT, lower, comm_world_bcast);
  }
#endif
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::Init() {
  // Set the total amount of work to do per scan
  if (progress_table_)
    progress_table_->set_total_work(opts_->leapfrog_steps * num_states_);

#ifdef PARALLEL
  // Init members needed for parallel tempering
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
  MPI_Comm_size(MPI_COMM_WORLD, &num_replicas_);
#endif

  // Set geometrically increasing temperature for this replica
  temperature_ = opts_->temperature * pow(opts_->alpha, myrank_);

  // Initialize replicas vector with initial increasing order
  for (int r = 0; r < num_replicas_; ++r)
    replicas_.push_back(r);

  // Setup random number generators by generating seeds with the build-in rand()
  // method which in turn is seeded with the user provided seed.
  const gsl_rng_type* T = gsl_rng_default;
  srand(opts_->seed + myrank_);
  rng_local_ = gsl_rng_alloc(T);
  gsl_rng_set(rng_local_, rand());
  srand(opts_->seed);
  rng_replica_exch_ = gsl_rng_alloc(T);
  gsl_rng_set(rng_replica_exch_, rand());
  rng_participants_ = gsl_rng_alloc(T);
  gsl_rng_set(rng_participants_, rand());
  rng_swap_ = gsl_rng_alloc(T);
  gsl_rng_set(rng_swap_, rand());

  // Generate vector with data indices in random order
  Batch idx;
  for (int n = 0; n < num_data_; ++n) idx.push_back(n);
  random_shuffle(idx.begin(), idx.end());
  // Setup batches for online training by distributing indices randomly into batches
  batches_.clear();
  const int batch_size = iround(static_cast<float>(num_data_) / opts_->num_batches);
  for (int b = 0; b < opts_->num_batches; ++b) {
    Batch batch;
    // Last batch may differ in length
    const int end = (b == opts_->num_batches - 1) ? num_data_ : (b+1) * batch_size;
    for (int n = b * batch_size; n < end; ++n)
      batch.push_back(idx[n]);
    batches_.push_back(batch);
  }
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::SaveToMemento() {
  memento_.crf   = crf_;
  memento_.grad1 = grad1_;
  memento_.grad2 = grad2_;
  memento_.logp  = logp_;
  memento_.logr  = logr_;
}

template< class Alphabet, template<class> class Subject >
void Hmc<Alphabet, Subject>::RestoreFromMemento() {
  crf_   = memento_.crf;
  grad1_ = memento_.grad1;
  grad2_ = memento_.grad2;
  logp_  = memento_.logp;
  logr_  = memento_.logr;
}



template< class Alphabet, template<class> class Subject >
HmcProgressTable<Alphabet, Subject>::HmcProgressTable(
    const Hmc<Alphabet, Subject>* hmc,
    FILE* fout,
    int width)
    : ProgressTable(fout, width),
      hmc_(hmc) {}

template< class Alphabet, template<class> class Subject >
void HmcProgressTable<Alphabet, Subject>::PrintHeader() {
  fprintf(fout_, "%-3s %5s  %7s %8s  %-30s%9s %10s  %5s %4s  %5s %5s %7s %7s\n",
          "No", "Temp", "Epsilon", "Best LL", "Leapfrog", "Sample LL",
          "Prior", "Prob", "Acc", "Hist", "Scans", "Time", "Dist");
  fputs(std::string(120, '-').c_str(), fout_);
  fputc('\n', fout_);
}

template< class Alphabet, template<class> class Subject >
void HmcProgressTable<Alphabet, Subject>::PrintRowBegin() {
  Reset();
  fprintf(fout_, "%-3i %5.4g  %7.5f %8.4f  ", hmc_->iter_, hmc_->temperature_,
          hmc_->epsilon_, hmc_->best_->logp / hmc_->num_data_);
  fflush(fout_);
}

template< class Alphabet, template<class> class Subject >
void HmcProgressTable<Alphabet, Subject>::PrintRowEnd() {
  fprintf(fout_, " %8.4f %10.4f  %5.1f %4s  %5.1f %5d %7.3f %7.2f\n",
          hmc_->likelihood_normalized(),
          hmc_->regularizer_normalized(),
          100 * MetropolisProb(hmc_->hamiltonian_beg_, hmc_->hamiltonian_end_),
          hmc_->samples_.back()->accepted ? "yes" : "no",
          100.0 * hmc_->num_accepted_ / hmc_->iter_,
          hmc_->num_scans_, hmc_->integration_time_,
          EuclideanDistance(hmc_->memento_.crf, hmc_->crf_));
  fflush(fout_);
}



template<class Alphabet>
inline double MatchProb(const CrfState<Alphabet>& state,
                        const CountProfile<Alphabet>& count_profile) {
  assert(!count_profile.logspace());

  const int alphabet_size = state.alphabet_size();
  const int num_cols      = state.num_cols();
  double rv = state.bias_weight();

  for(int j = 0; j < num_cols; ++j) {
    for (int a = 0; a < alphabet_size; ++a)
      rv += count_profile[j][a] * state[j][a];
  }
  return exp(rv);
}

template<class Alphabet>
inline double MatchProb(const CrfState<Alphabet>& state,
                        const Sequence<Alphabet>& seq) {
  const int num_cols = state.num_cols();
  double rv = state.bias_weight();
  for(int j = 0; j < num_cols; ++j) rv += state[j][seq[j]];
  return exp(rv);
}

inline double MetropolisProb(double hamiltonian_beg, double hamiltonian_end) {
  return MIN(1., exp(hamiltonian_beg - hamiltonian_end));
}

inline double ReplicaExchangeProb(double temp_i, double logp_i,
                                  double temp_j, double logp_j) {
  double p = exp((logp_j - logp_i) * ((1.0 / temp_i) - (1.0 / temp_j)));
  return MIN(1.0, p);
}

template<class Alphabet>
double EuclideanDistance(const Crf<Alphabet>& crf_p, const Crf<Alphabet>& crf_q) {
  const int num_states = crf_p.num_states();
  const int num_cols   = crf_p.num_cols();
  const int alph_size  = Alphabet::instance().size();

  double sum = 0.0;
  for (int k = 0; k < num_states; ++k) {
    const CrfState<Alphabet>& p = crf_p[k];
    const CrfState<Alphabet>& q = crf_q[k];

    sum += (p.bias_weight() - q.bias_weight()) *
      (p.bias_weight() - q.bias_weight());
    for (int a = 0; a < alph_size; ++a)
      sum += (p(a) - q(a)) * (p(a) - q(a));
    for (int j = 0; j < num_cols; ++j)
      for (int a = 0; a < alph_size; ++a)
        sum += (p[j][a] - q[j][a]) * (p[j][a] - q[j][a]);
  }

  return sqrt(sum);
}

template<class Alphabet>
double EuclideanNorm(const Crf<Alphabet>& crf) {
  const int num_states = crf.num_states();
  const int num_cols   = crf.num_cols();
  const int alph_size  = Alphabet::instance().size();

  double sum = 0.0;
  for (int k = 0; k < num_states; ++k) {
    const CrfState<Alphabet>& state = crf[k];

    sum += state.bias_weight() * state.bias_weight();
    for (int a = 0; a < alph_size; ++a)
      sum += state(a) * state(a);
    for (int j = 0; j < num_cols; ++j)
      for (int a = 0; a < alph_size; ++a)
        sum += state[j][a] * state[j][a];
  }

  return sqrt(sum);
}


}  // namespace cs

#endif  // CS_HMC_INL_H_
