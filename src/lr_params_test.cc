#include <gtest/gtest.h>
#include "cs.h"
#include "lr_params-inl.h"

namespace cs {

const char* const TMPFILE = "../data/test/lr_params.lr";

TEST(LrParamsTest, AccTest) {
	LrCfg<AA> cfg(3, 2);

	LrCfg<AA> cfg2(3, 3);
	LrCfg<AA> cfg3(cfg);
	ASSERT_TRUE(cfg != cfg2);
	ASSERT_TRUE(cfg == cfg3);
	
	LrParams<AA> params(cfg);
	float f[AA::kSize][cfg.nparams_let];
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) {
			f[a][p] = rand();
			params[a][p] = f[a][p];
		}
	}
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) {
			ASSERT_EQ(params[a][p], f[a][p]);
		}
	}
}


TEST(LrParamsTest, IOTest) {
	LrCfg<AA> cfg(5, 3);
	LrParams<AA> params(cfg, GaussianLrParamsInit<AA>(1));
	FILE* fout = fopen(const_cast<const char*>(TMPFILE), "w");
	params.Write(fout);
	fclose(fout);
	LrParams<AA> params_r(fopen(TMPFILE, "r"));
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) 
			ASSERT_NEAR(params[a][p], params_r[a][p], 0.000001);
	}
}

TEST(LrParamsTest, CopyTest) {
	LrCfg<AA> cfg(13, 4);
	GaussianLrParamsInit<AA> init(10);
	LrParams<AA> params(cfg, init);
	LrParams<AA> params_c = params;
	LrParams<AA> params_cc(cfg);
	params_cc = params;

	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++)
			params_c[a][p]++;
	}
	for (size_t a = 0; a < AA::kSize; a++) {
		for (size_t p = 0; p < cfg.nparams_let; p++) {
			ASSERT_EQ(params[a][p] + 1, params_c[a][p]);
			ASSERT_EQ(params[a][p], params_cc[a][p]);
		}
	}

}

}  // namespace cs

int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


