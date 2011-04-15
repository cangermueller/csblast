#include <gtest/gtest.h>
#include "cs.h"
#include "lr_pseudocounts-inl.h"


namespace cs {


TEST(LrPseudocountsTest, LrPseudocountsTest) {
	FILE* fin = fopen("../data/test/lr_pseudocounts.lr", "r");
	LrParams<AA> params(fin);
	LrPseudocounts<AA> pc(params);
	ConstantAdmix adm(0);
	Sequence<AA> seq("ahahlaksasdnceinc");
	Profile<AA> p(seq.length());
	pc.AddToSequence(seq, adm, p);
	for (size_t i = 0; i < p.length(); i++) {
		double sum = 0.0;
		for (size_t a = 0; a < AA::kSize; a++) 
			sum += p[i][a];
		ASSERT_NEAR(sum, 1.0, 0.00001);
	}
	std::cout << p << "\n";
}		







}  // namespace cs

int main(int argc, char** argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}


