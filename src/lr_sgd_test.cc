#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "lr_sgd-inl.h"
#include "training_sequence.h"
#include "lr_pseudocounts-inl.h"

namespace cs {

class LrSgdTest : public testing::Test {
 protected:
  typedef std::vector<TrainingSequence<AA> > TrainingSet;

  static const size_t kWlen = 13;
  static const size_t kKlen = 4;
  static const size_t kTsetSize = 100;

  LrSgdTest() : cfg_(kWlen, kKlen) ,file_(new char[100]) {}
  ~LrSgdTest() { delete[] file_; }

  virtual void SetUp() {
	sprintf(file_, "../data/test/lr_sgd/lr_sgd_%zu_%zu_%zu.lr", 
			kWlen, kKlen, kTsetSize);
    srand(0);
    const size_t c = (kWlen - 1) / 2;
    for (size_t n = 0; n < kTsetSize; ++n) {
      // Generate input sequence by sampling from Blosum background freqs
      Sequence<AA> seq(kWlen);
      seq.set_header(strprintf("seq%u", n + 1));
      for (size_t i = 0; i < kWlen; ++i) {
        double r = rand() / double(RAND_MAX);
        double sum = 0.0f;
        for (size_t a = 0; a < AA::kSize; ++a) {
          sum += sm_.p(a);
          if (r <= sum) { seq[i] = a; break; }
        }
      }
      // Set pseudocounts to conditional probabilities P(a|b) given aa 'b' in
      // central position of input sequence window.
      ProfileColumn<AA> col;
      for (size_t a = 0; a < AA::kSize; ++a)
        col[a] = sm_.r(a, seq[c]);

      tset_.push_back(TrainingSequence<AA>(seq, col));
    }
	sparams_.max_epochs = 150;
  }

  LrCfg<AA> cfg_;
  TrainingSet tset_;
  BlosumMatrix sm_;
  LrSgdParams sparams_;
  char* const file_;
};

TEST_F(LrSgdTest, Learning) {
  GaussianLrParamsInit<AA> init(0.1);
  LrParams<AA> params(cfg_, init);
  DerivLrFunc<AA, TrainingSequence<AA> > func(cfg_, tset_, sm_);
  LrSgdOptimizer<AA, TrainingSequence<AA> > sgd_opt(func, func, sparams_);
  double loglike = sgd_opt.Optimize(params, stdout);
  printf("loglike: %1.4f\n", loglike);
  FILE* fout = fopen(file_, "w");
  params.Write(fout);
  fclose(fout);
}

TEST_F(LrSgdTest, DISABLED_Prediction) {
	//FILE* fin = fopen(file_, "r");
	FILE* fin = fopen("../data/lr/nr30_neff2.5_hhblast_1round_neff6_min2_W13_N5000000_s1.lr", "r");
	LrParams<AA> params(fin);
	fclose(fin);
	ConstantAdmix adm(1.0);
	LrPseudocounts<AA> pc(params);
	Sequence<AA> seq("ARNDCQEGHILKMFPSTWYV");
	Profile<AA> profile = pc.AddTo(seq, adm);

	printf("Profile for sequence '%s':\n", seq.ToString().c_str());
	printf("%7sexpected actual  diff\n", "");  
	for (size_t i = 0; i < seq.length(); i++) {
		for (size_t a = 0; a < AA::kSize; a++)
			printf("%c -> %c: %5.3f   %5.3f   %5.3f\n",
					seq.chr(i), AA::kIntToChar[a],
					sm_.r(a, seq[i]), profile[i][a], 
					fabs(sm_.r(a, seq[i]) - profile[i][a]));
	}
	size_t center = kWlen / 2;
	for (size_t i = center; i < seq.length() - center; i++) {
		for (size_t a = 0; a < AA::kSize; a++)
			EXPECT_NEAR(sm_.r(a, seq[i]), profile[i][a], 0.02);
	}
}



}  // namespace cs
