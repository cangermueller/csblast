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
// Foo debug

#ifndef CS_SGD_H_
#define CS_SGD_H_

#include "crf-inl.h"
#include "func.h"
#include "progress_bar.h"
#include "crf_pseudocounts-inl.h"

namespace cs {

using std::string;

// Parameter wrapper for various SGD parameters
struct SgdParams {
    SgdParams()
            : nblocks(1000),
              eta_mode(ETA_MODE_FUNC),
              eta_init(0.1),
              eta_decay(2.0),
              eta_reinit(0.1),
              eta_reinit_num(0),
              eta_reinit_delta(0.05),
              mu(0.01),
              rho(0.5),
              max_eta(1.0),
              gamma(0.9),
              toll(0.001),
              early_delta(0.01),
              min_epochs(1),
              max_epochs(150),
              min_ll(0.0),
              min_ll_repeats(0),
              sigma_bias(1.0),
              sigma_context(0.6),
              sigma_context_pos_min(0.1),
              sigma_context_pos_max(0.6),
              sigma_decay(0.8),
              sigma_pc_min(0.1),
              sigma_pc_max(1.0),
              sigma_relax_epoch(0),
              sigma_relax_steps(3),
              context_penalty(0.0),
              context_penalty_epoch(0),
              context_penalty_steps(3),
              seed(0) {}

    size_t nblocks;               // number of training blocks
    size_t eta_mode;              // mode for updating the learning rate eta
    double eta_init;              // initial learning rate eta
    double eta_decay;             // decay of the function for updating the learning rate eta
    double eta_reinit;            // learning rate eta for reinitialization
    size_t eta_reinit_num;        // number of reinitializations
    double eta_reinit_delta;      // threshold of LL change for reinitializing eta
    double mu;                    // meta learning rate
    double rho;                   // lower bound for multiplier in learning rate adaption
    double max_eta;               // upper bound for learning rates
    double gamma;                 // parameter governing exponential average of derivatives
    double toll;                  // LL change for convergence
    double early_delta;           // Deviation from the maximal LL on the validation set for early-stopping
    size_t min_epochs;            // minimal number of epochs
    size_t max_epochs;            // maximal number of epochs
    double min_ll;                // minimal LL on training set
    size_t min_ll_repeats;        // maximal number of repetitions of the first epoch
    double sigma_bias;            // sigma in prior of bias weights
    double sigma_context;         // sigma in prior of context weights
    double sigma_context_pos_min; // minimum sigma in unsymmetric prior of positive context weights
    double sigma_context_pos_max; // maximum sigma in unsymmetric prior of positive context weights
    double sigma_decay;           // exponential decay of sigma in prior of context weights
    double sigma_pc_min;          // minimum sigma in prior of pseudocounts weights
    double sigma_pc_max;          // maximum sigma in prior of pseudocounts weights
    size_t sigma_relax_epoch;     // epoch for beginning to relax the prior
    size_t sigma_relax_steps;     // number of steps for relaxing the prior
    double context_penalty;       // penalty that each context weight column refers to a density distribution
    size_t context_penalty_epoch; // epoch for relaxing the context penalty
    size_t context_penalty_steps; // number of steps for relaxing the context penalty
    unsigned int seed;            // seed for rng that shuffles training set after each epoch

    static const size_t ETA_MODE_ALAP = 1;  // Use ALAP3 for updating the learning rate
    static const size_t ETA_MODE_FUNC = 2;  // Use a function for updateing the learning rate
};


template<class Abc>
struct SgdState : public DerivCrfFuncIO<Abc> {
    SgdState(const Crf<Abc>& c)
            : DerivCrfFuncIO<Abc>(c),
              steps(0),
              eta(c.nweights(), 0.0),
              avg(c.nweights(), 0.0),
              grad_prev(c.nweights(), 0.0) {}

    size_t steps;              // number of SGD steps already performed
    Vector<double> eta;        // learning rates eta for each CRF weight
    Vector<double> avg;        // exponential average of learning rates
    Vector<double> grad_prev;  // previous gradient of likelihood and prior combined
};


template<class Abc, class TrainingPair>
struct Sgd {
    Sgd(const DerivCrfFunc<Abc, TrainingPair>& tf,
        const SgdParams& p)
            : func(tf),
              params(p),
              eta_fac(static_cast<double>((p.eta_decay - 1) * tf.trainset.size()) / 
                  (1e6 * p.nblocks)),
              eta_reinit(false),
              ran(p.seed) {}

    // Shuffles training set and then runs one epoche of stochastic gradient descent
    // comprising of 'nblocks' individual gradient descent steps.
    void operator() (SgdState<Abc>& s, ProgressBar* prog_bar = NULL) {
        s.loglike = 0.0;
        s.prior = 0.0;
        // Shuffle training set before each epoch
        random_shuffle(func.shuffle.begin(), func.shuffle.end(), ran);

        for (size_t b = 0; b < params.nblocks; ++b) {
            // Save previous gradient of likelihood and prior as combined vector
            for (size_t i = 0; i < s.grad_prev.size(); ++i)
                s.grad_prev[i] = s.grad_loglike[i] + s.grad_prior[i];
            // Calculate gradient and increment likelihood based on training block 'b'
            func.df(s, b, params.nblocks, prog_bar);
            // Udpate learning rates
            UpdateLearningRates(s);
            // Updates CRF weights based on learning rates and gradient
            UpdateCRF(s);

            //printf("%02zu prior: %.4g\n", s.steps, func.prior(s.crf) / s.crf.nweights());
            s.steps++;
        }
    }

    // Moves CRF weights along the gradient direction with parameter-specific
    // step sizes given in 'eta' vector.
    void UpdateCRF(SgdState<Abc>& s) {
        // Calculate delta vector
        Vector<double> delta(s.eta.size());
        double delta_max = 0.0;
        for (size_t i = 0; i < s.eta.size(); ++i) {
            delta[i] = s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
            delta_max = MAX(delta_max, fabs(delta[i]));
        }
        // Scale delta vector to the maximum parameter change threshold
        if (delta_max > kDeltaMax) {
            for (size_t i = 0; i < delta.size(); ++i)
                delta[i] = kDeltaMax * delta[i] / delta_max;
        }
        // Update CRF parameters
        for (size_t k = 0, i = 0; k < s.crf.size(); ++k) {
            s.crf[k].bias_weight += delta[i]; 
            ++i;
            for (size_t j = 0; j < s.crf.wlen(); ++j)
                for (size_t a = 0; a < Abc::kSize; ++a, ++i)
                    s.crf[k].context_weights[j][a] += delta[i];
            for (size_t a = 0; a < Abc::kSize; ++a, ++i) 
                s.crf[k].pc_weights[a] += delta[i];
            UpdatePseudocounts(s.crf[k]);
        }
    }

    // Updates exponential averages with new gradient
    void UpdateAverages(SgdState<Abc>& s) {
        double tmp = s.steps == 0 ? 1.0 : 1.0 - params.gamma;
        for (size_t i = 0; i < s.avg.size(); ++i)
            s.avg[i] = params.gamma * s.avg[i] + tmp * SQR(s.grad_loglike[i] + s.grad_prior[i]);
    }

    // Updates all learning rates
    void UpdateLearningRates(SgdState<Abc>& s) {
        if (s.steps == 0 || eta_reinit) {
            eta_init = eta_reinit ? params.eta_reinit : params.eta_init;
            eta_init_step = s.steps;
            Assign(s.eta, eta_init);
            eta_reinit = false;
        } else {
            if (params.eta_mode == SgdParams::ETA_MODE_ALAP) {                                  
                UpdateAverages(s);
                double tmp = 0.0;
                for (size_t i = 0; i < s.eta.size(); ++i) {
                    tmp = (s.grad_loglike[i] + s.grad_prior[i]) * s.grad_prev[i] / s.avg[i];
                    s.eta[i] = MIN(params.max_eta, s.eta[i] * MAX(params.rho, 1.0 + params.mu * tmp));
                }
            } else {
                double eta = eta_init / ((s.steps - eta_init_step) * eta_fac + 1);
                Assign(s.eta, eta);
            }
        }
    }

    static const double kDeltaMax = 1000; // Maximum parameter change per SGD iteraton
    DerivCrfFunc<Abc, TrainingPair> func; // training set function
    const SgdParams& params;              // SGD parameter    
    const double eta_fac;                 // Parameter for calculating the decay of the learning rate
    bool eta_reinit;                      // Indicates whether eta is to be reinitialized
    double eta_init;                      // Eta used for (re)initialization
    size_t eta_init_step;                 // Step when eta was (re)initialized
    Ran ran;                              // RNG for shuffling of training set
};


template<class Abc, class TrainingPairT, class TrainingPairV>
struct SgdOptimizer {
    SgdOptimizer(const DerivCrfFunc<Abc, TrainingPairT>& tf,
                 const CrfFunc<Abc, TrainingPairV>& vf,
                 const SgdParams& p,
                 const vector<CountProfile<Abc> >& ns = vector<CountProfile<Abc> >(),
                 const double npc = 0.0,
                 CrfInit<Abc>* ci = NULL) 
            : sgd(tf, p),
              func(vf),
              params(p),
              neff_samples(ns),
              neff_pc(npc),
              crf_init(ci) { } 

    double Optimize(Crf<Abc>& crf, FILE* fout = NULL) { 

        scoped_ptr<ProgressBar> prog_bar;
        SgdState<Abc> s(crf);
        size_t epoch = 1, max_vepoch = 1, nconv = 0, nearly = 0, neta = 0, neta_reinit = 0, nmin_ll = 0;
        double max_tloglike = -DBL_MAX, max_vloglike = -DBL_MAX;
        double val_loglike, old_loglike, init_loglike;
        std::string best_line;

        size_t status_len = 80;
        if (fout) {
            prog_bar.reset(new ProgressBar(fout, 16));
            fprintf(fout, "%-5s %-16s %9s %9s %9s %9s %9s %7s\n", 
                "Epoch", "Gradient descent", "LL-Train", "+/-", "Prior", "LL-Val", "Neff", "Eta");
            fprintf(fout, "%s\n", std::string(status_len, '-').c_str());
        }

        // Compute the initial likelihood
        init_loglike =  sgd.func(s.crf) / sgd.func.trainset.size();
        s.loglike = init_loglike;
        while (((nconv < kMaxConvBumps && nearly < kMaxConvBumps) || epoch <= params.min_epochs) && epoch <= params.max_epochs) {
            // Print first part of table row
            if (fout) {
                fprintf(fout, "%-4zu  ", epoch); fflush(fout);
                prog_bar->Init((sgd.func.trainset.size() + 1) * crf.size());
            }
            
            // Save last likelihood for calculation of delta
            old_loglike = s.loglike;
            // Update sigma in prior for pseudocounts weights
            DerivCrfFuncPrior<Abc>& prior = sgd.func.prior;
            UnsymmetricDerivCrfFuncPrior<Abc>* uprior = dynamic_cast<UnsymmetricDerivCrfFuncPrior<Abc>* >(&sgd.func.prior);
            if (params.sigma_relax_epoch == 0 || epoch >= params.sigma_relax_epoch + params.sigma_relax_steps) {
                prior.sigma_pc = params.sigma_pc_max;
                if (uprior) uprior->sigma_context_pos = params.sigma_context_pos_max;
            } else if (epoch >= params.sigma_relax_epoch) {
                prior.sigma_pc = params.sigma_pc_min + 
                    (params.sigma_pc_max - params.sigma_pc_min) / params.sigma_relax_steps * 
                    (1 + epoch - params.sigma_relax_epoch);
                if (uprior) 
                    uprior->sigma_context_pos = params.sigma_context_pos_min + 
                        (params.sigma_context_pos_max - params.sigma_context_pos_min) / params.sigma_relax_steps * 
                        (1 + epoch - params.sigma_relax_epoch);
            } else {
                prior.sigma_pc = params.sigma_pc_min;
                if (uprior)
                    uprior->sigma_context_pos = params.sigma_context_pos_min;
            }
            // Update context penalty
            if (params.context_penalty_epoch == 0 || epoch < params.context_penalty_epoch) {
              prior.context_penalty = params.context_penalty;
            } else if (epoch >= params.context_penalty_epoch + params.context_penalty_steps - 1) {
              prior.context_penalty = 0.0;
            } else {
              prior.context_penalty = params.context_penalty - (1 + epoch - params.context_penalty_epoch) * 
                                      params.context_penalty / params.context_penalty_steps;
            }

            // Run on epoche of SGD
            sgd(s, prog_bar.get());

            // Normalize likelihood and prior to user friendly scale
            s.loglike /= sgd.func.trainset.size();
            s.prior /= crf.nweights();
            // Calculate likelihood on validation set
            val_loglike = func(s.crf) / func.trainset.size();
            // Calculate delta for convergence
            double delta = s.loglike - old_loglike;
            // Keep track of how many times we were under convergence threshold
            if (delta > params.toll) nconv = 0;
            else ++nconv;
            // Keep track of how many times we were under the maximal likelihood
            if (val_loglike > max_vloglike - params.early_delta) nearly = 0;
            else ++nearly;
            // Keep track of how many times we were under threshold for reinitializing eta
            if (delta > params.eta_reinit_delta) neta = 0;
            else {
                if (++neta >= kMaxConvBumps && neta_reinit < params.eta_reinit_num) {
                    sgd.eta_reinit = true;
                    ++neta_reinit;
                    nconv = 0;
                    nearly = 0;
                    neta = 0;
                }
            }

            // Save CRF with the maximum likelihood on the validation set
            if (val_loglike >= max_vloglike) {
                max_vloglike = val_loglike;
                max_vepoch = epoch;
                crf = s.crf;
                if (!crffile_vset.empty()) {
                    FILE* fout = fopen(crffile_vset.c_str(), "w");
                    if (!fout) throw Exception("Can't write to file '%s'!", crffile_vset.c_str());
                    crf.Write(fout);
                    fclose(fout);
                }
            }

            // Save CRF with the maximum likelihood on the validation set
            if (s.loglike >= max_tloglike) {
                max_tloglike = s.loglike;
                if (!crffile_tset.empty()) {
                    FILE* fout = fopen(crffile_tset.c_str(), "w");
                    if (!fout) throw Exception("Can't write to file '%s'!", crffile_tset.c_str());
                    s.crf.Write(fout);
                    fclose(fout);
                }
            }

            // Compute the Neff for the given samples
            ConstantAdmix admix(neff_pc);
            CrfPseudocounts<Abc> pc(s.crf);
            double neff = 0.0;
            if (neff_samples.size() > 0) {
              int nsamples = static_cast<int>(neff_samples.size());
#pragma omp parallel for schedule(static)
              for (int i = 0; i < nsamples; ++i) {
                double n = Neff(pc.AddTo(neff_samples[i], admix));
#pragma omp atomic
                neff += n;
              }
              neff /= neff_samples.size();
            }
            if (fout) prog_bar->Complete();

            // Print second part of table row
            if (fout) {
                double eta = s.eta[0];
                if (params.eta_mode == SgdParams::ETA_MODE_ALAP) {
                    double sum = 0.0;
                    for (size_t i = 0; i < s.eta.size(); ++i)
                        sum += s.eta[i];
                    eta /= s.eta.size();
                }
                char line[1000];
                sprintf(line, " %9.4f %+9.4f %9.4f %9.4f %9.4f %7.2g",
                        s.loglike, delta, s.prior, val_loglike, neff, eta);
                fprintf(fout, "%s\n", line);
                if (val_loglike == max_vloglike) best_line = line;

            }
	
            // Repeat the first epoch if LL is to low
            if (epoch == 1 && s.loglike < params.min_ll) {
              if (nmin_ll < params.min_ll_repeats) {
                  s = SgdState<Abc>(crf);
                  if (crf_init != NULL) (*crf_init)(s.crf);
                  nconv = 0; 
                  nearly = 0;
                  neta = 0;
                  s.loglike = init_loglike;
                  ++nmin_ll;
              } else {
                break;
              }
            } else {
              epoch++;
            }
        }
        if (fout && !best_line.empty()) {
            fprintf(fout, "%s\n", std::string(status_len, '-').c_str());
            fprintf(fout, "%-4zu %16s %s\n", max_vepoch, "", best_line.c_str());
        }
        return max_vloglike;
    }
	

    static const size_t kMaxConvBumps = 5;

    Sgd<Abc, TrainingPairT> sgd;      // SGD algorithm encapsulation
    CrfFunc<Abc, TrainingPairV> func; // validation set function
    const SgdParams& params;          // SGD parameters
    const vector<CountProfile<Abc> >&
        neff_samples;                 // samples for computing the Neff
    const double neff_pc;             // pseudocounts admix for computing the Neff
    CrfInit<Abc>* crf_init;           // Object for reinitializing the CRF
    string crffile_vset;              // Output file for best CRF on the validation set
    string crffile_tset;              // Output file for last CRF on the training set
};


}  // namespace cs

#endif  // CS_SGD_H_
