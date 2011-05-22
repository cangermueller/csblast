#include <gtest/gtest.h>
#include "cs.h"
#include "lr_func-inl.h"
#include "sequence-inl.h"
#include "blosum_matrix.h"


namespace cs {

using std::cout;
using std::vector;

float err = 0.0001;

TEST(LrFuncTest, LrKmerBuilderTest) {
	LrCfg<AA> cfg(5, 3);
	LrKmerBuilder<AA> b(cfg);
	Sequence<AA> s("ahers");
	vector<long> kmers = b(s);
	ASSERT_EQ(kmers.size(), 3);
	ASSERT_EQ(kmers[0], 166);
	ASSERT_EQ(kmers[0], b.at(s, 0)); 
	ASSERT_EQ(kmers[1], 11321); 
	ASSERT_EQ(kmers[1], b.at(s, 1)); 
	ASSERT_EQ(kmers[2], 18435); 
	ASSERT_EQ(kmers[2], b.at(s, 2)); 

	size_t center = cfg.wlen / 2;
	for (size_t i = 0; i < cfg.wlen; i++) {
		vector<long> kmers_win = b(s, i - center);
		if (i == 0) {
			ASSERT_EQ(kmers_win.size(), 1);
			ASSERT_EQ(kmers_win[0], 16166);
		} else if (i == 1) {
			ASSERT_EQ(kmers_win.size(), 2);
			ASSERT_EQ(kmers_win[0], 8166);
			ASSERT_EQ(kmers_win[1], 19321);
		} else if (i == 2) {
			ASSERT_EQ(kmers_win.size(), kmers.size());
			for (size_t j = 0; j < kmers.size(); j++)
				ASSERT_EQ(kmers_win[j], kmers[j]);
		} else if (i == 3) {
			ASSERT_EQ(kmers_win.size(), 2);
			ASSERT_EQ(kmers_win[0], 3321);
			ASSERT_EQ(kmers_win[1], 10435);
		} else if (i == 4) {
			ASSERT_EQ(kmers_win.size(), 1);
			ASSERT_EQ(kmers_win[0], 2435);
		}
	}
	vector<long> kmers_win;
	kmers_win = b(s, -3);
	ASSERT_EQ(kmers_win.size(), 0);
	kmers_win = b(s, 3);
	ASSERT_EQ(kmers_win.size(), 0);
	kmers_win = b(s, 10);
	ASSERT_EQ(kmers_win.size(), 0);
	
}

TEST(LrFuncTest, DISABLED_LrFuncInitTest) {
	vector< TrainingSequence<AA> > tset;
	ReadAll(fopen("../data/test/lr_func_init.tsq", "r"), tset);
	LrCfg<AA> cfg(13, 4);
	BlosumMatrix sm(BLOSUM62);
	LrKmerBuilder<AA> builder(cfg);

	LrFunc<AA, TrainingSequence<AA> > func(cfg, tset, sm);
	for (size_t t = 0; t < tset.size(); t++) {
		vector<long> kmers = builder(tset[t].x);
		for (size_t k = 0; k < cfg.nkmers; k++) {
			ASSERT_EQ(kmers[k], (*func.tset_kmers)[t][k]);
		}
	}
}

TEST(LrFuncTest, DISABLED_LrFuncTest) {
	vector< TrainingSequence<AA> > tset;
	ReadAll(fopen("../data/test/lr_func.tsq", "r"), tset);
	LrCfg<AA> cfg(13, 4);
	BlosumMatrix sm(BLOSUM62);
	LrKmerBuilder<AA> builder(cfg);
	LrParams<AA> params(cfg);
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) 
			params[a][p] = static_cast<long>(((a + 1) * p) % 5) * (p % 2 ? 1 : -1);
	}
	LrFunc<AA, TrainingSequence<AA> > func(cfg, tset, sm);

	double ll = 0;
	for (size_t n = 0; n < tset.size(); n++) {
		vector<long> kmers = builder(tset[n].x);
		double scalar_exp_sum = 0;
		double neff = 0;
		for (size_t a = 0; a < AA::kSize; a++) {
			double scalar = 0;
			for (size_t i = 0; i < kmers.size(); i++)
				scalar += static_cast<long>(((a + 1) * kmers[i]) % 5) * 
					(kmers[i] % 2 ? 1 : -1);
			assert(scalar < 50);
			ll += tset[n].y[a] * scalar;
			scalar_exp_sum += exp(scalar);
			neff += tset[n].y[a];
		}
		ll -= neff * log(scalar_exp_sum);
	}
	ll += func.ll_norm;
	ASSERT_NEAR(ll, func(params), err);

	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) 
			params[a][p] = a + 1;
	}
	ll = 0;
	for (size_t n = 0; n < tset.size(); n++) {
		double scalar_exp_sum = 0;
		double neff = 0;
		for (size_t a = 0; a < 20; a++) {
			ll += tset[n].y[a] * 10 * (a + 1);
			scalar_exp_sum += exp(10 * (a + 1));
			neff += tset[n].y[a];
		}
		ll -= neff * log(scalar_exp_sum);
	}
	ll += func.ll_norm;
	ASSERT_NEAR(ll, func(params), err);
}

TEST(LrFuncTest, DISABLED_LrFuncDerivTest) {
	vector< TrainingSequence<AA> > tset;
	ReadAll(fopen("../data/test/lr_func_deriv.tsq", "r"), tset);
	LrCfg<AA> cfg(13, 4);
	BlosumMatrix sm(BLOSUM62);
	LrKmerBuilder<AA> builder(cfg);
	LrParams<AA> params(cfg);
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) 
			params[a][p] = a + 1;
	}
	float sigma = 10;
	DerivLrFunc<AA, TrainingSequence<AA> > func(cfg, tset, sm, sigma);
	DerivLrFuncIO<AA> s(params);
	s.loglike = 0;
	s.prior = 0;
	func.df(s, 0, 1);	

	double ll = 0;
	for (size_t n = 0; n < tset.size(); n++) {
		double scalar_exp_sum = 0;
		double neff = 0;
		for (size_t a = 0; a < 20; a++) {
			ll += tset[n].y[a] * 10 * (a + 1);
			scalar_exp_sum += exp(10 * (a + 1));
			neff += tset[n].y[a];
		}
		ll -= neff * log(scalar_exp_sum);
	}
	ll += func.ll_norm;
	ASSERT_NEAR(ll, s.loglike, err);

	double scalar_exp[AA::kSize];
	double scalar_exp_sum = 0;
	double neff = 0;
	for (size_t a = 0; a < AA::kSize; a++) {
		scalar_exp[a] = exp((a + 1) * 10);
		scalar_exp_sum += scalar_exp[a];
		neff += tset[0].y[a];
	}
	size_t i = 0;
	for (size_t a = 0; a < AA::kSize; a++) {
		float grad = tset[0].y[a] - neff * scalar_exp[a] / scalar_exp_sum;
		for (size_t p = 0; p < cfg.nparams_let; p++) {
			if (s.grad_loglike[i] != 0)
				ASSERT_NEAR(grad, s.grad_loglike[i], err);
			i++;
		}
	}

	for (size_t a = 0; a < AA::kSize; a++) {
		size_t offset = a * cfg.nparams_let;
		float ref = s.grad_prior[offset];
		for (size_t p = 1; p < cfg.nparams_let; p++)
			ASSERT_FLOAT_EQ(ref, s.grad_prior[offset + p]);
	}

	double prior = 0;
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) prior += SQR(params[a][p]);
	}
	prior /= -(2 * SQR(sigma));
	ASSERT_NEAR(prior, s.prior, err);

}		


}  // namespace cs

