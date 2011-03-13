#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "context_library-inl.h"
#include "count_profile-inl.h"
#include "em_clustering.h"
#include "matrix_pseudocounts-inl.h"

namespace cs {

const double kDelta = 0.001;

class EMClusteringTestGCN4 : public testing::Test {
 protected:
  typedef std::vector<CountProfile<AA> > TrainingSet;

  virtual void SetUp() {
    // Read full length training data profiles
    FILE* fp = fopen("../data/gcn4_W13.prf", "r");
    ReadAll(fp, trainset_);
    fclose(fp);
  }

  TrainingSet trainset_;
};

TEST_F(EMClusteringTestGCN4, InitBySampling) {
  BlosumMatrix m;
  MatrixPseudocounts<AA> pc(m);
  ConstantAdmix admix(1.0);
  SamplingLibraryInit<AA> init(trainset_, pc, admix);
  ContextLibrary<AA> lib(7, 13, init);

  EMClustering<AA> em(trainset_, lib, m, 1.6, 0.85, 1e-5);
  for (int i = 0; i < 2000; ++i) {
    em.EStep();
    em.MStep();
  }

  EXPECT_NEAR(0.154798, em.loglike, kDelta);
}

TEST_F(EMClusteringTestGCN4, InitByGaussian) {
  BlosumMatrix m;
  MatrixPseudocounts<AA> pc(m);
  ConstantAdmix admix(1.0);
  GaussianLibraryInit<AA> init(0.3, m, 123);
  ContextLibrary<AA> lib(7, 13, init);

  EMClustering<AA> em(trainset_, lib, m, 1.6, 0.85, 1e-5);
  int iters = RunClustering(em, 1e-4, 1000, stdout);
  EXPECT_NEAR(0.155049, em.loglike, kDelta);
  EXPECT_EQ(22, iters);
}

}  // namespace cs
