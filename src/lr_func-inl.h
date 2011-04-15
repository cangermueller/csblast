// Copyright 2010, Andreas Biegert
// Author	: Angermueller Christof 
//			: angermueller@lmb.uni-muenchen.de
// Date		: 2010-08-24		

#ifndef CS_LR_FUNC_INL_H_
#define CS_LR_FUNC_INL_H_

#include "lr_func.h"


namespace cs {


// LrKmerBuilder


template<class Abc>
long LrKmerBuilder<Abc>::at(const Sequence<Abc>& s, size_t kpos, int wpos) const {
	long kmer = kpos * pow[cfg.klen];
	for (size_t i = 0; i < cfg.klen; i++)
		kmer += s[wpos + kpos + i] * pow[cfg.klen - i - 1];
	return kmer;
}

template<class Abc>
vector<long> LrKmerBuilder<Abc>::operator() (const Sequence<Abc>& s, int wpos) const {
	size_t kfirst = wpos < 0 ? fabs(wpos) : 0;
	int nkmers = MIN(wpos + cfg.wlen, s.length()) - MAX(0, wpos) - cfg.klen + 1;
	if (nkmers < 0) nkmers = 0;
	vector<long> kmers(nkmers);
	for (int i = 0; i < nkmers; i++)
		kmers[i] = at(s, static_cast<int>(kfirst + i), wpos);
	return kmers;
}


// LrFunc


template<class Abc, class TrainingPair>
LrFunc<Abc, TrainingPair>::LrFunc(
		const LrCfg<Abc>& cfg_, 
		const TrainingSet& tset_, 
		const SubstitutionMatrix<Abc>& sm_) : 
	cfg(cfg_), tset(tset_), tset_kmers(new Matrix<long>(tset.size(), cfg.nkmers)), sm(sm_) {

	// Calculate kmers for each training sequence.
	LrKmerBuilder<Abc> builder(cfg);
	long* ts_kmers = tset_kmers->begin();
	typename vector<TrainingPair>::const_iterator t = tset.begin();
	int nkmers_size = cfg.nkmers * sizeof(long);
	while (t != tset.end()) {
		vector<long> t_kmers = builder(t->x);
		memcpy(ts_kmers, &t_kmers[0], nkmers_size);
		t++;
		ts_kmers += cfg.nkmers;
	}
	// Calculate the likelihood normalization term.
	ll_norm = 0.0;
	for (size_t a = 0; a < Abc::kSize; a++) {
		double pc_sum = 0.0;
		typename vector<TrainingPair>::const_iterator t = tset.begin();
		while (t != tset.end()) {
			pc_sum += t->y[a];
			t++;
		}
		ll_norm += log(sm.p(a)) * pc_sum;
	}	
	ll_norm = -ll_norm;
}

template<class Abc, class TrainingPair>
LrFunc<Abc, TrainingPair>::LrFunc(const LrFunc<Abc, TrainingPair>& func) :
	cfg(func.cfg), tset(func.tset), tset_kmers(new Matrix<long>(tset.size(), cfg.nkmers)),
	sm(func.sm), ll_norm(func.ll_norm) { *tset_kmers = *func.tset_kmers; }


template<class Abc, class TrainingPair>
double LrFunc<Abc, TrainingPair>::operator() (const LrParams<Abc>& params) const {
	double ll = 0;
	size_t ntset = tset.size();
#pragma omp parallel for schedule(static)
	for (size_t n = 0; n < ntset; n++) {
		double ll_n = 0;
		double scalar_exp_sum = 0;
		double neff = 0; 
		long* kmers = (*tset_kmers)[n];
		const ProfileColumn<Abc>& pc = tset[n].y;
		for (size_t a = 0; a < Abc::kSize; a++) {
			double scalar = 0;
			for (size_t k = 0; k < cfg.nkmers; k++)
				scalar += params[a][kmers[k]];
			ll_n += pc[a] * scalar;
			assert(fabs(scalar) < 500);
			scalar_exp_sum += exp(scalar);
			neff += pc[a];
		}
		ll_n -= neff * log(scalar_exp_sum);
#pragma omp atomic
		ll += ll_n;
	}
	return ll + ll_norm;
}


// DerivLrFunc


template<class Abc, class TrainingPair>
TrainingBlock DerivLrFunc<Abc, TrainingPair>::GetBlock(size_t b, size_t nblocks) const {
	assert(b < nblocks);
	size_t block_size = iround(static_cast<float>(tset.size()) / nblocks);
	size_t beg = b * block_size;
	size_t end = (b == nblocks - 1) ? tset.size() : (b+1) * block_size;
	size_t size = end - beg;
	float frac = static_cast<double>(size) / tset.size();
	return TrainingBlock(beg, end, size, frac);
}    

template<class Abc, class TrainingPair>
void DerivLrFunc<Abc, TrainingPair>::df(DerivLrFuncIO<Abc>& s, size_t b, size_t nblocks) const {

	assert(b < nblocks);
	const TrainingBlock block(GetBlock(b, nblocks));
	const size_t n_beg = block.beg;
	const size_t n_end = block.end;

	double ll = 0;
#pragma omp parallel for schedule(static)
	for (size_t n = n_beg; n < n_end; n++) {
		long tset_i = shuffle[n];
		double ll_n = 0;
		double scalar_exp_sum = 0;
		double neff = 0; 
		long* kmers = (*tset_kmers)[tset_i];
		const ProfileColumn<Abc>& pc = tset[tset_i].y;
		for (size_t a = 0; a < Abc::kSize; a++) {
			double scalar = 0;
			for (size_t k = 0; k < cfg.nkmers; k++)
				scalar += s.params[a][kmers[k]];
			ll_n += pc[a] * scalar;
			assert(fabs(scalar) < 500);	// overflow if scalar > 720
			scalar_exp_sum += exp(scalar);
			neff += pc[a];
		}
		ll_n -= neff * log(scalar_exp_sum);
#pragma omp atomic
		ll += ll_n;
	}
	s.loglike += ll + ll_norm * block.frac;	// normalization depends on the block fraction
	s.prior   += block.frac * CalculatePrior(s.params);

	CalculateLikelihoodGradient(block, s);
	CalculatePriorGradient(block, s);
}

template<class Abc, class TrainingPair>
void DerivLrFunc<Abc, TrainingPair>::CalculateLikelihoodGradient(
		const TrainingBlock& block, DerivLrFuncIO<Abc>& s) const {

	memset(&s.grad_loglike[0], 0, cfg.nparams * sizeof(float));
#pragma omp parallel for schedule(static)
	for (size_t n = block.beg; n < block.end; n++) {
		long tset_i = shuffle[n];
		double scalar_exp[Abc::kSize];
		double scalar_exp_sum = 0;
		double neff = 0;
		long* kmers = (*tset_kmers)[tset_i];
		const ProfileColumn<Abc>& pc = tset[tset_i].y;
		for (size_t a = 0; a < Abc::kSize; a++) {
			scalar_exp[a] = 0;
			for (size_t k = 0; k < cfg.nkmers; k++) 
				scalar_exp[a] += s.params[a][kmers[k]];
			assert(fabs(scalar_exp[a]) < 500);
			scalar_exp[a] = exp(scalar_exp[a]);
			scalar_exp_sum += scalar_exp[a];
			neff += pc[a];
		}
		double factor = neff / scalar_exp_sum;
		for (size_t a = 0; a < Abc::kSize; a++) {
			long offset = a * cfg.nparams_let;
			float grad = static_cast<float>(pc[a] - scalar_exp[a] * factor);
			for (size_t k = 0; k < cfg.nkmers; k++) 
#pragma omp atomic
				s.grad_loglike[offset + kmers[k]] += grad;
		}
	}
}

template<class Abc, class TrainingPair>
void DerivLrFunc<Abc, TrainingPair>::CalculatePriorGradient(
		const TrainingBlock& block, DerivLrFuncIO<Abc>& s) const {

	float factor = -block.frac / SQR(sigma);
	float* p = s.params.params().begin();
	for (size_t i = 0; i < cfg.nparams; i++, p++)
		s.grad_prior[i] = *p * factor;
}


template<class Abc, class TrainingPair>
double DerivLrFunc<Abc, TrainingPair>::CalculatePrior(const LrParams<Abc>& params) const {

	double prior = 0;
	float* p = params.params().begin();
	for (size_t i = 0; i < cfg.nparams; i++, p++)
		prior += SQR(*p);
	return prior * -0.5 / SQR(sigma);
}



}	// namespace cs

#endif	// CS_LR_FUNC_INL_H_
