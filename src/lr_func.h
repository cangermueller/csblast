// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-24		

#ifndef CS_LR_FUNC_H_
#define CS_LR_FUNC_H_

#include "lr_params-inl.h"
#include "substitution_matrix-inl.h"
#include "training_sequence.h"
	

namespace cs {

using std::vector;


// LrKmerBuilder


// Calculates unique kmer identifiers.
template<class Abc>
struct LrKmerBuilder {
	LrKmerBuilder(const LrCfg<Abc>& cfg_) : cfg(cfg_), pow(new long[cfg.klen + 1]) {
		pow[0] = 1;
		for (size_t i = 1; i <= cfg.klen; ++i) pow[i] = pow[i - 1] * Abc::kSize;
	}

	~LrKmerBuilder() { delete[] pow; }

	// Calculates kmer identifier for kmer at position kpos in window at position wpos.
	long at(const Sequence<Abc>& s, size_t kpos, int wpos = 0) const;

	// Calculates all kmer identifiers for window at position wpos.
	vector<long> operator() (const Sequence<Abc>& s, int wpos = 0) const;


	// Variables
	const LrCfg<Abc> cfg;	// window config
	long* const pow;		// pow[i] = Abc::kSize^i, i = [0, klen] 
};


// TrainingBlock


// Stores information of a training block.
struct TrainingBlock {
    TrainingBlock() : beg(0), end(0), size(0), frac(0) {}

    TrainingBlock(size_t b, size_t e, size_t s, float f)
            : beg(b), end(e), size(s), frac(f) {}


	// Variables	
    size_t beg;		// index of the first sequence
    size_t end;		// index of the last sequence
    size_t size;	// number of sequences
    float frac;		// size / total number of sequences
};


// LrFunc


/*
Likelihood function
L = sum{n = 1 to N} ( sum{a = A} tn(a)*wa*phi(xn) ) -
    sum{n = 1 to N} ( tn*log( sum{a = 1 to A} exp(wa*phi(xn)) ) ) - 
    N*(sum{a = 1 to A} log(p(a)))	// Normalization term

	n		: index of a training sequence
	N		: total number of training sequences
	xn		: training sequence n
	phi(xn)	: returns kmer vector of the given sequence
	    	: phi(xn)[k] = 1 <=> kmer k occurs in sequence xn
    tn(a)	: pseudo count of training sequence n
	tn		: tn = sum {a = 1 to A} tn(a)
	a		: index of a letter 
	A		: total number of letters of the alphabet
	wa		: kmer weight vector 
	p(a)	: equilibrium frequency of a letter
*/
template<class Abc, class TrainingPair>
struct LrFunc {
    typedef vector<TrainingPair> TrainingSet;

    LrFunc(
			const LrCfg<Abc>& cfg_,
			const TrainingSet& tset_, 
			const SubstitutionMatrix<Abc>& sm_);

	LrFunc(const LrFunc<Abc, TrainingPair>& func);

    virtual ~LrFunc() { delete tset_kmers; }

	// Calculates likelihood for the given parameters.
    double operator() (const LrParams<Abc>& params) const;


	// Variables
	const LrCfg<Abc> cfg;				// window configs
    const TrainingSet& tset;			// training set
	Matrix<long>* const tset_kmers;		// kmers of each training sequence
	const SubstitutionMatrix<Abc>& sm;	// substitution matrix
	double ll_norm;						// likelihood normalization term
};


// DerivLrFuncIO


// Input output parameters for DerivLrFuncIO
template<class Abc>
struct DerivLrFuncIO {
    DerivLrFuncIO(const LrParams<Abc>& params_) : 
		params(params_), 
		grad_loglike(params_.cfg().nparams, 0.0), 
		grad_prior(params_.cfg().nparams, 0.0), 
		loglike(-DBL_MAX), 
		prior(-DBL_MAX) {}

    virtual ~DerivLrFuncIO() {}


	// variables
    const LrParams<Abc>& params;	// current parameters
    Vector<float> grad_loglike;  	// partial derivatives of LL
    Vector<float> grad_prior;    	// partial derivatives of prior
    double loglike;               	// log of conditional likelihood
    double prior;                 	// log of prior probability
};


// DerivLrFunc


// Derivation of the likelihood function
template<class Abc, class TrainingPair>
struct DerivLrFunc : public LrFunc<Abc, TrainingPair> {
    typedef vector<TrainingPair> TrainingSet;

    DerivLrFunc(
			const LrCfg<Abc>& cfg_,
			const TrainingSet& tset_, 
			const SubstitutionMatrix<Abc>& sm_, 
			float sigma_ = 1.0) :  
		LrFunc<Abc, TrainingPair>::LrFunc(cfg_, tset_, sm_),
		shuffle(tset_.size()), 
		sigma(sigma_) { 
		for (size_t i = 0; i < tset.size(); i++) shuffle[i] = i;
    }

	// Returns a training block.
    TrainingBlock GetBlock(size_t b, size_t nblocks) const;

	// Calculates partial derivation on the given training block.
    void df(DerivLrFuncIO<Abc>& s, size_t b = 0, size_t nblocks = 1) const;

	/*
	Calculates parital derivation of the likelihood function.
	dL/dwaf = sum { 1 to N} (tn(a)-tn*p(a|xn,wa))*phi(xn)[f]

		n		: index of a training sequence
		N		: total number of training sequences
		xn		: training sequence n
		f		: index of a feature (kmer)
		F		: total number of features (kmers)
		phi(xn)	: returns kmer vector of the given sequence
				: phi(xn)[f] = 1 <=> kmer f occurs in sequence xn
		tn(a)	: pseudo count of training sequence n
		tn		: tn = sum {a = 1 to A} tn(a)
		a		: index of a letter 
		A		: total number of letters of the alphabet
		wa		: kmer weight vector 
	*/
    void CalculateLikelihoodGradient(
			const TrainingBlock& block, DerivLrFuncIO<Abc>& s) const;

	// Calculates parital derivation of the prior.
	// dP/dwaf = -(1/sigma^2)*waf*block.frac
    void CalculatePriorGradient(
			const TrainingBlock& block, DerivLrFuncIO<Abc>& s) const;

	// Calculates the prior.
	// P = -1/(2*sigma^2)*sum {a = 1 to A} ( sum {f = 1 to F} wa[f]^2 )
    double CalculatePrior(const LrParams<Abc>& params) const;


	// Variables 
    using LrFunc<Abc, TrainingPair>::cfg;			// window config
    using LrFunc<Abc, TrainingPair>::tset;			// training set
    using LrFunc<Abc, TrainingPair>::tset_kmers;	// kmers of each training sequence
    using LrFunc<Abc, TrainingPair>::sm;			// substitution matrix
    using LrFunc<Abc, TrainingPair>::ll_norm;		// likelihood function normalization term
    vector<int> shuffle;							// indices of the shuffled training set
	const float sigma;								// sigma for calculating the prior
};	// DerivLrFunc


}	// namespace cs

#endif	// CS_LR_FUNC_H_
