// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-30		

#ifndef CS_LR_SGD_INL_H_
#define CS_LR_SGD_INL__H_

#include "lr_func-inl.h"
#include "progress_bar.h"


namespace cs {


// LrSgdParams


// Stores sgd parameters.
struct LrSgdParams {
    LrSgdParams()
            : nblocks(100),
              eta(0.001),
              mu(0.01),
              rho(0.5),
              max_eta(1.0),
              gamma(0.9),
              toll(1e-3),
              min_epochs(10),
              max_epochs(500),
              seed(0) {}


	// Variables
    size_t nblocks;     // number of training blocks
    double eta;         // learning rate eta for context and pseudocount weights
    double mu;          // meta learning rate
    double rho;         // lower bound for multiplier in learning rate adaption
    double max_eta;     // upper bound for learning rates
    double gamma;       // parameter governing exponential average of derivatives
    double toll;        // LL change for convergence
    int min_epochs;     // minimal number of epochs
    int max_epochs;     // maximal number of epochs
    unsigned int seed;  // seed for rng that shuffles training set after each epoch
};


// LrSgdState


// Store a sgd state.
template<class Abc>
struct LrSgdState : public DerivLrFuncIO<Abc> {
    LrSgdState(const LrParams<Abc>& p)
            : DerivLrFuncIO<Abc>(p),
              steps(0),
              eta(p.cfg().nparams, static_cast<float>(0)),
              avg(p.cfg().nparams, static_cast<float>(0)),
              grad_prev(p.cfg().nparams, static_cast<float>(0)) {}


	// Variables
    size_t steps;             // number of SGD steps already performed
    Vector<float> eta;        // learning rates eta for each parameter 
    Vector<float> avg;        // exponential average of learning rates
    Vector<float> grad_prev;  // previous gradient of likelihood and prior combined
};


// LrSgd


// Performs a sgd step.
template<class Abc, class TrainingPair>
struct LrSgd {
    LrSgd(const DerivLrFunc<Abc, TrainingPair>& func_, const LrSgdParams& params)
            : func(func_),
              nblocks(params.nblocks),
              eta(params.eta),
              mu(params.mu),
              rho(params.rho),
              max_eta(params.max_eta),
              gamma(params.gamma),
              ran(params.seed) {}

    // Shuffles training set and then runs one epoche of stochastic gradient descent
    // comprising of 'nblocks' individual gradient descent steps.
    void operator() (LrSgdState<Abc>& s, ProgressBar* prog_bar = NULL) {
        s.loglike = 0.0;
        s.prior = 0.0;
        // Shuffle training set before each epoch
        random_shuffle(func.shuffle.begin(), func.shuffle.end(), ran);
		if (prog_bar) prog_bar->Init(nblocks);

        for (size_t b = 0; b < nblocks; b++) {
            // Save previous gradient of likelihood and prior as combined vector
            for (size_t i = 0; i < s.grad_prev.size(); i++)
                s.grad_prev[i] = s.grad_loglike[i] + s.grad_prior[i];
            // Calculate gradient and increment likelihood based on training block 'b'
            func.df(s, b, nblocks);
            // Udpate averages of partial derivatives
            UpdateAverages(s);
            // Udpate learning rates
            UpdateLearningRates(s);
            // Updates parameters based on learning rates and gradient
            UpdateLrParams(s);

            s.steps++;
			if (prog_bar) prog_bar->Advance();
        }
    }

    // Moves LrParams weights along the gradient direction with parameter-specific
    // step sizes given in 'eta' vector.
    void UpdateLrParams(LrSgdState<Abc>& s) {
		float* p = s.params.params().begin();
		size_t nparams = s.params.cfg().nparams;
		for (size_t i = 0; i < nparams; i++, p++)
			*p += s.eta[i] * (s.grad_loglike[i] + s.grad_prior[i]);
	}

    // Updates all learning rates based on the current gradient, the previous
    // gradient and the exponential average vector.
    void UpdateLearningRates(LrSgdState<Abc>& s) {
        if (s.steps == 0) {
            Assign(s.eta, static_cast<float>(eta)); // start with constant eta
        } else {
            double tmp = 0.0;
            for (size_t i = 0; i < s.eta.size(); ++i) {
                tmp = (s.grad_loglike[i] + s.grad_prior[i]) * s.grad_prev[i] / s.avg[i];
                s.eta[i] = MIN(max_eta, s.eta[i] * MAX(rho, 1.0 + mu * tmp));
            }
        }
    }

    // Updates exponential averages with new gradient
    void UpdateAverages(LrSgdState<Abc>& s) {
        double tmp = s.steps == 0 ? 1.0 : 1.0 - gamma;
        for (size_t i = 0; i < s.avg.size(); ++i)
            s.avg[i] = gamma * s.avg[i] + tmp * SQR(s.grad_loglike[i] + s.grad_prior[i]);
    }


	// Variables
    DerivLrFunc<Abc, TrainingPair> func; 	// training set function
    size_t nblocks;     // number of training blocks
    double eta;         // learning rate eta for context and pseudocount weights
    double mu;          // meta learning rate
    double rho;         // lower bound for multiplier in learning rate adaption
    double max_eta;     // upper bound for learning rates
    double gamma;       // parameter governing exponential average of derivatives
    Ran ran;            // RNG for shuffling of training set
};	// LrSgd


// LrSgdOptimizer


// Optimizes the likelihood the validation set by sgd.
template<class Abc, class TrainingPair>
struct LrSgdOptimizer {
    LrSgdOptimizer(const DerivLrFunc<Abc, TrainingPair>& dfunc,
                 const LrFunc<Abc, TrainingPair>& func_,
                 const LrSgdParams& params)
            : sgd(dfunc, params),
              func(func_),
              toll(params.toll),
              min_epochs(params.min_epochs),
              max_epochs(params.max_epochs) {}

	// Optimizes the likelihood by altering params.
	// Returns the best likelihood value achieved on the validation set.
    double Optimize(LrParams<Abc>& params, FILE* fout = NULL) {
        scoped_ptr<ProgressBar> prog_bar;
		LrParams<Abc> p(params);
        LrSgdState<Abc> s(p);
        int epoch = 0, nconv = 0;
        double best_loglike = -DBL_MAX;
		double delta = DBL_MAX;
		double old_loglike;
		double val_loglike;


        if (fout) {
            prog_bar.reset(new ProgressBar(fout, 16));
            fprintf(fout, "%-5s %-16s  %10s %10s %10s %10s\n", "Epoch",
                    "Gradient descent", "LL-Train", "+/-", "Prior", "LL-Val");
            fprintf(fout, "%s\n", std::string(67, '-').c_str());
        }

        while ((nconv < kMaxConvBumps || epoch < min_epochs) && epoch < max_epochs) {
            // Print first part of table row
            if (fout) {
				fprintf(fout, "%-4d  ", epoch + 1); 
				fflush(fout);
			}

            // Save last likelihood for calculation of delta
            old_loglike = s.loglike;
            // Run on epoche of SGD
            sgd(s, prog_bar.get());
            // Normalize likelihood and prior to user friendly scale
            s.loglike /= sgd.func.tset.size();
            s.prior /= params.cfg().nparams;
            delta = s.loglike - old_loglike;
            // Keep track of how many times we were under convergence threshold
            if (fabs(delta) > toll) nconv = 0;
            else ++nconv;
            // Calculate likelihood on validation set
            val_loglike = func(s.params) / func.tset.size();
            if (val_loglike > best_loglike) {
                best_loglike = val_loglike;
                params = s.params;
            }
            // Print second part of table row
            if (fout) {
                fprintf(fout, "  %10.4f %10.4f %10.4f %10.4f\n",
                        s.loglike, epoch == 0 ? 0.0 : delta, s.prior, val_loglike);
            }
            epoch++;
        }

        return best_loglike;
    }


	// Variables
    static const int kMaxConvBumps = 3;

    LrSgd<Abc, TrainingPair> sgd;       // SGD algorithm encapsulation
    LrFunc<Abc, TrainingPair> func;		// validation set function
    double toll;                      	// LL change for convergence
    int min_epochs;                   	// minimal number of epochs
    int max_epochs;                   	// maximal number of epochs
};	// LrSgdOptimizer


}  // namespace cs

#endif  // CS_LR_SGD_INL_H_
