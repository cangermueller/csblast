// Copyright 2009, Andreas Biegert

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
              eta_init(0.001),
              eta_mode(ETA_MODE_ALAP),
              eta_decay(2),
              mu(0.01),
              rho(0.5),
              max_eta(1.0),
              gamma(0.9),
              toll(1e-3),
              early_delta(0.01),
              min_epochs(10),
              max_epochs(150),
              sigma_pc_min(0.01),
              sigma_pc_max(1.0),
              sigma_pc_epoch(0),
              sigma_pc_delta(0.002),
              sigma_pc_steps(20),
              seed(0) {}

    size_t nblocks;        // number of training blocks
    double eta_init;       // initial learning rate eta
    size_t eta_mode;       // mode for updating the learning rate eta
    double eta_decay;      // decay of the function for updating the learning rate eta
    double mu;             // meta learning rate
    double rho;            // lower bound for multiplier in learning rate adaption
    double max_eta;        // upper bound for learning rates
    double gamma;          // parameter governing exponential average of derivatives
    double toll;           // LL change for convergence
    double early_delta;    // Deviation from the maximal LL on the validation set for early-stopping
    size_t min_epochs;     // minimal number of epochs
    size_t max_epochs;     // maximal number of epochs
    double sigma_pc_min;   // Minimum sigma in prior of pseudocounts weights
    double sigma_pc_max;   // Maximum sigma in prior of pseudocounts weights
    size_t sigma_pc_epoch; // Epoch for beginning to relax the prior of pseudocounts weights
    double sigma_pc_delta; // Gradient for beginning to relax the prior of pseudocount weights
    size_t sigma_pc_steps; // Number of steps for relaxing the prior of pseudocounts weights
    unsigned int seed;     // seed for rng that shuffles training set after each epoch

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

            s.steps++;
        }
    }

    // Moves CRF weights along the gradient direction with parameter-specific
    // step sizes given in 'eta' vector.
    void UpdateCRF(SgdState<Abc>& s) {
        for (size_t k = 0, i = 0; k < s.crf.size(); ++k) {
            s.crf[k].bias_weight += s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
            ++i;
            for (size_t j = 0; j < s.crf.wlen(); ++j)
                for (size_t a = 0; a < Abc::kSize; ++a, ++i)
                    s.crf[k].context_weights[j][a] += s.eta[i] * (s.grad_loglike[i] +
                                                                  s.grad_prior[i]);
            for (size_t a = 0; a < Abc::kSize; ++a, ++i) 
                s.crf[k].pc_weights[a] += s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
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
        if (s.steps == 0) {
            Assign(s.eta, params.eta_init); // start with constant eta
        } else {
            if (params.eta_mode == SgdParams::ETA_MODE_ALAP) {
                UpdateAverages(s);
                double tmp = 0.0;
                for (size_t i = 0; i < s.eta.size(); ++i) {
                    tmp = (s.grad_loglike[i] + s.grad_prior[i]) * s.grad_prev[i] / s.avg[i];
                    s.eta[i] = MIN(params.max_eta, s.eta[i] * MAX(params.rho, 1.0 + params.mu * tmp));
                }
            } else {
                double eta = params.eta_init / (s.steps * (params.eta_decay - 1) / params.nblocks + 1);
                Assign(s.eta, eta);
            }
        }
    }

    DerivCrfFunc<Abc, TrainingPair> func; // training set function
    const SgdParams& params;              // SGD parameter
    Ran ran;                              // RNG for shuffling of training set
};


template<class Abc, class TrainingPairT, class TrainingPairV>
struct SgdOptimizer {
    SgdOptimizer(const DerivCrfFunc<Abc, TrainingPairT>& tf,
                 const CrfFunc<Abc, TrainingPairV>& vf,
                 const SgdParams& p,
                 const vector<CountProfile<Abc> >& ns = vector<CountProfile<Abc> >(),
                 const double npc = 0.0) 
            : sgd(tf, p),
              func(vf),
              params(p),
              neff_samples(ns),
              neff_pc(npc) { } 

    double Optimize(Crf<Abc>& crf, FILE* fout = NULL) { 

        scoped_ptr<ProgressBar> prog_bar;
        SgdState<Abc> s(crf);
        size_t epoch = 1, best_epoch = 1, nconv = 0, nearly = 0, nsigma_pc = 0;
        double best_loglike = -DBL_MAX, delta = DBL_MAX, old_loglike, val_loglike;
        size_t sigma_pc_epoch = params.sigma_pc_epoch;
        std::string best_line;

        if (fout) {
            prog_bar.reset(new ProgressBar(fout, 16));
            fprintf(fout, "%-5s %-16s %8s %8s %9s %8s %8s\n", 
                "Epoch", "Gradient descent", "LL-Train", "+/-", "Prior", "LL-Val", "Neff");
            fprintf(fout, "%s\n", std::string(68, '-').c_str());
        }

        // Compute the initial likelihood
        s.loglike =  sgd.func(s.crf) / sgd.func.trainset.size();
        while (((nconv < kMaxConvBumps && nearly < kMaxEarlyBumps) || epoch <= params.min_epochs) 
            && epoch <= params.max_epochs) {
            // Print first part of table row
            if (fout) {
                fprintf(fout, "%-4zu  ", epoch); fflush(fout);
                prog_bar->Init((sgd.func.trainset.size() + 1) * crf.size());
            }
            
            // Save last likelihood for calculation of delta
            old_loglike = s.loglike;
            // Update sigma in prior for pseudocounts weights
            if (sigma_pc_epoch == 0 || epoch < sigma_pc_epoch) {
                sgd.func.prior.sigma_pc = params.sigma_pc_min;
            } else if (epoch >= sigma_pc_epoch + params.sigma_pc_steps) {
                sgd.func.prior.sigma_pc = params.sigma_pc_max;
            } else {
                sgd.func.prior.sigma_pc = params.sigma_pc_min + 
                    (params.sigma_pc_max - params.sigma_pc_min) / params.sigma_pc_steps * (epoch - sigma_pc_epoch);
            }
            // sgd.func.prior.sigma_pc = params.sigma_pc_min + (params.sigma_pc_max - params.sigma_pc_min) / 
            //     (1 + exp((static_cast<int>(params.sigma_pc_epoch - epoch)) / params.sigma_pc_delta));
            // Run on epoche of SGD
            sgd(s, prog_bar.get());
            // Normalize likelihood and prior to user friendly scale
            s.loglike /= sgd.func.trainset.size();
            s.prior /= crf.nweights();
            // Calculate likelihood on validation set
            val_loglike = func(s.crf) / func.trainset.size();
            // Calculate delta for convergence
            delta = s.loglike - old_loglike;
            // Keep track of how many times we were under convergence threshold
            if (delta > params.toll) nconv = 0;
            else ++nconv;
            // Keep track of how many times we were under the gradient for relaxing the pseudocount weights
            if (sigma_pc_epoch == 0) {
                if (fabs(delta) > params.sigma_pc_delta) nsigma_pc = 0;
                else if (++nsigma_pc >= kMaxConvBumps) sigma_pc_epoch = epoch;
            }
            // Keep track of how many times we were under the maximal likelihood
            if (val_loglike > best_loglike - params.early_delta) nearly = 0;
            else ++nearly;

            // Save CRF with the maximum likelihood on the validation set
            if (val_loglike >= best_loglike) {
                best_loglike = val_loglike;
                best_epoch = epoch;
                crf = s.crf;
                if (!crffile_vset.empty()) {
                    FILE* fout = fopen(crffile_vset.c_str(), "w");
                    if (!fout) throw Exception("Can't write to file '%s'!", crffile_vset.c_str());
                    crf.Write(fout);
                    fclose(fout);
                }
            }

            // Save CRF of the current epoch
            if (!crffile_tset.empty()) {
                FILE* fout = fopen(crffile_tset.c_str(), "w");
                if (!fout) throw Exception("Can't write to file '%s'!", crffile_tset.c_str());
                s.crf.Write(fout);
                fclose(fout);
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
            prog_bar->Complete();

            // Print second part of table row
            if (fout) {
                char line[100];
                sprintf(line, " %8.4f %+8.4f %9.4f %8.4f %8.4f",
                        s.loglike, delta, s.prior, val_loglike, neff);
                fprintf(fout, "%s\n", line);
                if (epoch == best_epoch) best_line = line;

            }
		
            epoch++;
        }
        if (fout && !best_line.empty()) {
            fprintf(fout, "%s\n", std::string(68, '-').c_str());
            fprintf(fout, "%-4zu %16s %s\n", best_epoch, "", best_line.c_str());
        }
        return best_loglike;
    }
	

    static const size_t kMaxConvBumps = 3;
    static const size_t kMaxEarlyBumps = 5;

    Sgd<Abc, TrainingPairT> sgd;      // SGD algorithm encapsulation
    CrfFunc<Abc, TrainingPairV> func; // validation set function
    const SgdParams& params;          // SGD parameters
    const vector<CountProfile<Abc> >&
        neff_samples;                 // samples for computing the Neff
    const double neff_pc;             // pseudocounts admix for computing the Neff
    string crffile_vset;              // Output file for best CRF on the validation set
    string crffile_tset;              // Output file for last CRF on the training set
};


}  // namespace cs

#endif  // CS_SGD_H_
