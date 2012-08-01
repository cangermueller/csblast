/*
  Copyright 2009-2012 Andreas Biegert, Christof Angermueller

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

#ifndef CS_SGD_INL_H_
#define CS_SGD_INL_H_

#include "sgd.h"

#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "globals.h"
#include "sum_product_algorithm-inl.h"
#include "crf-inl.h"
#include "log.h"
#include "shared_ptr.h"
#include "substitution_matrix-inl.h"
#include "utils.h"

namespace cs {

template< class Alphabet, template<class> class Subject >
Sgd<Alphabet, Subject>::Sgd(const DataVec& data,
                            const PredVec& pred,
                            const SgdOptions* opts,
                            const SubstitutionMatrix<Alphabet>* sm,
                            Crf<Alphabet>* crf)
        : alphabet_size_(Alphabet::instance().size()),
          num_cols_(crf->num_cols()),
          num_states_(crf->num_states()),
          num_data_(data.size()),
          data_(data),
          pred_(pred),
          opts_(opts),
          sm_(sm),
          crf_(*crf),
          d_pc_(crf->num_states(), Alphabet::instance().size(), 0.0),
          d_tr_(crf->num_states(), crf->num_states(), 0.0),
          d_pc_prev_(crf->num_states(), Alphabet::instance().size(), 0.0),
          d_tr_prev_(crf->num_states(), crf->num_states(), 0.0),
          eta_pc_(crf->num_states(), Alphabet::instance().size(), opts->eta),
          eta_tr_(crf->num_states(), crf->num_states(), opts->eta),
          v_pc_(crf->num_states(), Alphabet::instance().size(), 0.0),
          v_tr_(crf->num_states(), crf->num_states(), 0.0),
          num_data_cols_(0),
          logp_(0.0f),
          logp_prev_(0.0f),
          iters_(0),
          scan_(1),
          eta_(opts->eta) {
    Init();
}

template< class Alphabet, template<class> class Subject >
Sgd<Alphabet, Subject>::Sgd(const DataVec& data,
                            const PredVec& pred,
                            const SgdOptions* opts,
                            const SubstitutionMatrix<Alphabet>* sm,
                            Crf<Alphabet>* crf,
                            FILE* fout)
        : alphabet_size_(Alphabet::instance().size()),
          num_cols_(crf->num_cols()),
          num_states_(crf->num_states()),
          num_data_(data.size()),
          data_(data),
          pred_(pred),
          opts_(opts),
          sm_(sm),
          crf_(*crf),
          d_pc_(crf->num_states(), Alphabet::instance().size(), 0.0),
          d_tr_(crf->num_states(), crf->num_states(), 0.0),
          d_pc_prev_(crf->num_states(), Alphabet::instance().size(), 0.0),
          d_tr_prev_(crf->num_states(), crf->num_states(), 0.0),
          eta_pc_(crf->num_states(), Alphabet::instance().size(), opts->eta),
          eta_tr_(crf->num_states(), crf->num_states(), opts->eta_tr),
          v_pc_(crf->num_states(), Alphabet::instance().size(), 0.0),
          v_tr_(crf->num_states(), crf->num_states(), 0.0),
          num_data_cols_(0),
          logp_(0.0f),
          logp_prev_(0.0f),
          iters_(0),
          scan_(1),
          eta_(opts->eta) {
    progress_table_.reset(new SgdProgressTable<Alphabet, Subject>(this, fout));
    Init();
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::Run() {
    crf_.TransformTransitionsToLogSpace();
    if (progress_table_) progress_table_->PrintHeader();

    LOG(ERROR) << crf_;

    while (!IsDone()) {
        SetupBatches(opts_->num_batches);
        logp_prev_ = logp_;
        logp_      = 0.0;

        if (progress_table_) progress_table_->PrintRowBegin();

        for (BatchVec::iterator it = batches_.begin(); it != batches_.end(); ++it) {
            AddRegularizersToLikelihood(1.0f / opts_->num_batches);
            ResetDerivatives(1.0f / opts_->num_batches);

            CalculateGradient(*it);
            if (iters_ == 0)
                InitGradientAverage();
            else
                UpdateLearningRates();

            UpdateCRF();

            ++iters_;
        }
        if (progress_table_) progress_table_->PrintRowEnd();
        ++scan_;
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::CalculateGradient(const Batch& batch) {
    LOG(ERROR) << "Calculating gradient ...";
    crf_.TransformTransitionsToLinSpace();

    const int batch_size = batch.size();
#pragma omp parallel for schedule(static)
    for (int n = 0; n < batch_size; ++n) {
        SumProductMatrices spm(data_[batch[n]]->length(), num_states_);
        SumProductAlgorithm(crf_, *data_[batch[n]], *pred_[batch[n]], &spm);

        UpdateContextWeightDerivatives(spm, *data_[batch[n]]);
        UpdatePseudocountDerivatives(spm, *pred_[batch[n]]);
        UpdateTransitionDerivatives(spm);

#pragma omp atomic
        logp_ += (spm.logp_pc - spm.logp) / num_data_cols_;

        if (progress_table_) {
#pragma omp critical (print_progress)
            progress_table_->print_progress(num_states_ * data_[batch[n]]->length());
        }
    }

    crf_.TransformTransitionsToLogSpace();
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::UpdateContextWeightDerivatives(
        const SumProductMatrices& m,
        const CountProfile<Alphabet>& c) {
    const int slen   = c.length();
    const int center = (num_cols_ - 1) / 2;

    for (int k = 0; k < num_states_; ++k) {
        Matrix& w_k = *d_cw_[k];

        for (int i = 0; i < slen; ++i) {
            const int beg = MAX(0, i - center);
            const int end = MIN(c.num_cols() - 1, i + center);
            const double delta = (m.alpha_pc[i][k] * m.beta_pc[i][k]) -
                (m.alpha[i][k] * m.beta[i][k]);

            for(int h = beg; h <= end; ++h) {
                const int j = h - i + center;
                for (int a = 0; a < alphabet_size_; ++a) {
#pragma omp atomic
                    w_k[j][a] += delta * c.counts(h,a);
                }
            }
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::UpdatePseudocountDerivatives(
        const SumProductMatrices& m,
        const CountProfile<Alphabet>& c) {
    const int slen   = c.length();

    for (int k = 0; k < num_states_; ++k) {
        CrfState<Alphabet>& s_k = crf_[k];

        for (int i = 0; i < slen; ++i) {
            const double pp = m.alpha_pc[i][k] * m.beta_pc[i][k];

            for (int a = 0; a < alphabet_size_; ++a) {
#pragma omp atomic
                d_pc_[k][a] += pp * (c.counts(i,a) - s_k.pc(a) * c.neff(i));
            }
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::UpdateTransitionDerivatives(
        const SumProductMatrices& m) {

    const int slen = m.subject_length;  // length of subject
    int k;         // index of source state
    int l;         // index of target state
    double pp;     // posterior probability of going from k to l
    double pp_pc;  // posterior probability of going from k to l WITH pseudocounts

    for (ConstTransitionIter ti = crf_.transitions_begin();
         ti != crf_.transitions_end(); ++ti) {
        k = ti->source;
        l = ti->target;

        for (int i = 0; i < slen - 1; ++i) {
            pp = m.alpha[i][k] * m.beta[i+1][l] * ti->weight * m.context_prob[i+1][l] /
                m.alpha_sum[i+1];
            pp_pc = m.alpha_pc[i][k] * m.beta_pc[i+1][l] * ti->weight *
                m.context_prob[i+1][l] * m.pc_prob[i+1][l] / m.alpha_sum_pc[i+1];

#pragma omp atomic
            d_tr_[k][l] += pp_pc - pp;
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::UpdateCRF() {
    LOG(ERROR) << "Updating CRF ...";

#pragma omp parallel for schedule(static)
    for (int k = 0; k < num_states_; ++k) {
        CrfState<Alphabet>& s_k = crf_[k];
        Matrix& d_cw      = *d_cw_[k];
        Matrix& d_cw_prev = *d_cw_prev_[k];
        Matrix& eta_cw    = *eta_cw_[k];

        for (int a = 0; a < alphabet_size_; ++a) {
            s_k(a) += eta_pc_[k][a] * d_pc_[k][a];

            d_pc_prev_[k][a] = d_pc_[k][a];
            d_pc_[k][a] = 0.0;
        }

        for (int j = 0; j < num_cols_; ++j) {
            for (int a = 0; a < alphabet_size_; ++a) {
                s_k[j][a] += eta_cw[j][a] * d_cw[j][a];

                d_cw_prev[j][a] = d_cw[j][a];
                d_cw[j][a] = 0.0;
            }
        }

        s_k.UpdatePseudocounts();
    }

    // assert(crf_.transitions_logspace());
    // for (int k = 0; k < num_states_; ++k) {
    //   for (int l = 0; l < num_states_; ++l) {
    //     if (crf_.test_transition(k,l)) {
    //       crf_(k,l) = crf_(k,l) + eta_tr_[k][l] * d_tr_[k][l];

    //       d_tr_prev_[k][l] = d_tr_[k][l];
    //       d_tr_[k][l] = 0.0;
    //     }
    //   }
    // }

    crf_.increment_iterations();
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::SetupBatches(int num_batches) {
    batches_.clear();

    Batch idx;
    for (int n = 0; n < num_data_; ++n)
        idx.push_back(n);
    random_shuffle(idx.begin(), idx.end());

    // Sort into indexes batches
    const int batch_size = iround(static_cast<float>(num_data_) / num_batches);
    for (int b = 0; b < num_batches; ++b) {
        Batch batch;
        const int end = (b == num_batches - 1) ? num_data_ : (b+1) * batch_size;
        for (int n = b * batch_size; n < end; ++n)
            batch.push_back(idx[n]);
        batches_.push_back(batch);
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::ResetDerivatives(float w) {
    const float sigma_squared = opts_->sigma * opts_->sigma;

    assert(crf_.transitions_logspace());
    Reset(&d_tr_);  // FIXME: is this really needed?
    for (ConstTransitionIter ti = crf_.transitions_begin();
         ti != crf_.transitions_end(); ++ti) {
        d_tr_[ti->source][ti->target] = -w * ti->weight / sigma_squared;
    }

    for (int k = 0; k < num_states_; ++k) {
        CrfState<Alphabet>& s_k = crf_[k];
        Matrix& w_k = *d_cw_[k];

        for (int a = 0; a < alphabet_size_; ++a) {
            const float log_fa = fast_log2(sm_->f(a));

            d_pc_[k][a] = -w * (s_k(a) - log_fa) / sigma_squared;

            for (int j = 0; j < num_cols_; ++j) {
                w_k[j][a] = -w * (s_k[j][a] - log_fa) / sigma_squared;
            }
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::UpdateLearningRates() {
    LOG(ERROR) << "Updating learning rates ...";

    const double r = opts_->rho;
    const double m = opts_->mu;
    const double g = opts_->gamma;
    const double e = opts_->eta_max;
    int k,l = 0;
    double f = 0.0;

    for (ConstTransitionIter ti = crf_.transitions_begin();
         ti != crf_.transitions_end(); ++ti) {
        k = ti->source;
        l = ti->target;

        v_tr_[k][l] = g * v_tr_[k][l] + (1.0 - g) * d_tr_[k][l] * d_tr_[k][l];
        f = 1.0 + m * d_tr_[k][l] * d_tr_prev_[k][l] / v_tr_[k][l];
        eta_tr_[k][l] = MIN(e, eta_tr_[k][l] * MAX(r, f));

        //LOG(ERROR) << strprintf("eta_tr[%i][%i]=%7.2g", k, l, eta_tr_[k][l]);
    }

    for (int k = 0; k < num_states_; ++k) {
        for (int a = 0; a < alphabet_size_; ++a) {
            v_pc_[k][a] = g * v_pc_[k][a] + (1.0 - g) * d_pc_[k][a] * d_pc_[k][a];
            f = 1.0 + m * d_pc_[k][a] * d_pc_prev_[k][a] / v_pc_[k][a];
            eta_pc_[k][a] = MIN(e, eta_pc_[k][a] * MAX(r, f));

            // LOG(ERROR) << strprintf("eta_pc[%i][%i]=%7.2g", k, a, eta_pc_[k][a]);
        }
    }

    for (int k = 0; k < num_states_; ++k) {
        Matrix& d_cw      = *d_cw_[k];
        Matrix& d_cw_prev = *d_cw_prev_[k];
        Matrix& eta_cw    = *eta_cw_[k];
        Matrix& v_cw      = *v_cw_[k];

        for (int j = 0; j < num_cols_; ++j) {
            for (int a = 0; a < alphabet_size_; ++a) {
                v_cw[j][a] = g * v_cw[j][a] + (1.0 - g) * d_cw[j][a] * d_cw[j][a];
                f = 1.0 + m * d_cw[j][a] * d_cw_prev[j][a] / v_cw[j][a];
                eta_cw[j][a] = MIN(e, eta_cw[j][a] * MAX(r, f));
                // LOG(ERROR) << strprintf("eta_cw[%i][%i][%i]=%7.2g", k, j, a, eta_cw[j][a]);
            }
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::InitGradientAverage() {
    LOG(ERROR) << "Initializing average ...";

    int k,l = 0;
    for (ConstTransitionIter ti = crf_.transitions_begin();
         ti != crf_.transitions_end(); ++ti) {
        k = ti->source;
        l = ti->target;
        v_tr_[k][l] = d_tr_[k][l] * d_tr_[k][l];
    }

    for (int k = 0; k < num_states_; ++k)
        for (int a = 0; a < alphabet_size_; ++a)
            v_pc_[k][a] = d_pc_[k][a] * d_pc_[k][a];

    for (int k = 0; k < num_states_; ++k)     {
        Matrix& d_cw = *d_cw_[k];
        Matrix& v_cw = *v_cw_[k];

        for (int j = 0; j < num_cols_; ++j)
            for (int a = 0; a < alphabet_size_; ++a)
                v_cw[j][a] = d_cw[j][a] * d_cw[j][a];
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::AddRegularizersToLikelihood(float w) {
    const float sigma_squared = opts_->sigma * opts_->sigma;
    float x = 0.0f;

    assert(crf_.transitions_logspace());
    for (ConstTransitionIter ti = crf_.transitions_begin();
         ti != crf_.transitions_end(); ++ti) {
        logp_ -= w * ti->weight * ti->weight / (2 * sigma_squared);
    }

    for (int k = 0; k < num_states_; ++k) {
        CrfState<Alphabet>& s_k = crf_[k];

        for (int a = 0; a < alphabet_size_; ++a) {
            const float log_fa = fast_log2(sm_->f(a));

            x = s_k(a) - log_fa;
            logp_ -= w * x * x / (2 * sigma_squared);

            for (int j = 0; j < num_cols_; ++j) {
                x = s_k[j][a] - log_fa;
                logp_ -= w * x * x / (2 * sigma_squared);
            }
        }
    }
}

template< class Alphabet, template<class> class Subject >
void Sgd<Alphabet, Subject>::Init() {
    // Create matrices for context weight derivatives and learning rates
    for (int k = 0; k < num_states_; ++k) {
        MatrixPtr md(new Matrix(num_cols_, alphabet_size_, 0.0));
        d_cw_.push_back(md);
        MatrixPtr md_prev(new Matrix(num_cols_, alphabet_size_, 0.0));
        d_cw_prev_.push_back(md_prev);
        MatrixPtr mp(new Matrix(num_cols_, alphabet_size_, opts_->eta));
        eta_cw_.push_back(mp);
        MatrixPtr mv(new Matrix(num_cols_, alphabet_size_, 0.0));
        v_cw_.push_back(mv);
    }

    // Compute total number of data columns
    num_data_cols_ = 0;
    for (DataIter di = data_.begin(); di != data_.end(); ++di)
        num_data_cols_ += (**di).length();

    // Set the total amount of work to do per scan
    if (progress_table_)
        progress_table_->set_total_work(num_data_cols_ * num_states_);
}

template< class Alphabet, template<class> class Subject >
inline bool Sgd<Alphabet, Subject>::IsDone() const {
    if (scan_ < opts_->min_scans)
        return false;
    else if (scan_ >= opts_->max_scans)
        return true;
    else
        return fabs(rel_diff()) <= opts_->epsilon;
}


template< class Alphabet, template<class> class Subject >
SgdProgressTable<Alphabet, Subject>::SgdProgressTable(
        const Sgd<Alphabet, Subject>* sgd,
        FILE* fout,
        int width)
        : ProgressTable(fout, width),
          sgd_(sgd) {}

template< class Alphabet, template<class> class Subject >
void SgdProgressTable<Alphabet, Subject>::PrintHeader() {
    fprintf(fout_, "%-4s %4s %6s %3s %8s  %-30s  %9s  %8s\n", "Scan", "Itrs",
            "Conn", "B", "Eta", "Gradient calculation", "log(L)", "+/-");
    fputs(std::string(82, '-').c_str(), fout_);
    fputc('\n', fout_);
}

template< class Alphabet, template<class> class Subject >
void SgdProgressTable<Alphabet, Subject>::PrintRowBegin() {
    Reset();
    fprintf(fout_, "%-4i %4i %6.1f %3i %8.2g  ", sgd_->scan_, sgd_->iters_,
            sgd_->crf_.connectivity(), static_cast<int>(sgd_->batches_.size()),
            sgd_->eta_);
    fflush(fout_);
}

template< class Alphabet, template<class> class Subject >
void SgdProgressTable<Alphabet, Subject>::PrintRowEnd() {
    if (sgd_->scan() == 1)
        fprintf(fout_, "  %9.5f\n", sgd_->logp_);
    else
        fprintf(fout_, "  %9.5f  %+8.5f\n", sgd_->logp_, sgd_->rel_diff());
    fflush(fout_);
}

}  // namespace cs

#endif  // CS_SGD_INL_H_
