#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "context_library-inl.h"
#include "func.h"
#include "matrix_pseudocounts-inl.h"
#include "pdf_writer-inl.h"
#include "tamura_nei_matrix.h"
#include "training_sequence.h"

namespace cs {

const double kDelta = 0.01;

class ContextLibraryTest : public testing::Test {
  protected:

    virtual void SetUp() {
        FILE* fin = fopen("../data/context_profile.prf", "r");
        p1_.Read(fin);
        rewind(fin);
        p2_.Read(fin);
        rewind(fin);
        p3_.Read(fin);
        fclose(fin);
    }

    ContextProfile<AA> p1_;
    ContextProfile<AA> p2_;
    ContextProfile<AA> p3_;
};

TEST_F(ContextLibraryTest, SimpleConstruction) {
    ContextLibrary<AA> lib(3, 13);

    EXPECT_EQ((size_t)3, lib.size());
    EXPECT_EQ((size_t)13, lib.wlen());
    lib.SetProfile(0, p1_);
    lib.SetProfile(1, p2_);
    lib.SetProfile(2, p3_);

    EXPECT_NEAR(-6.20, lib[0].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[0].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[0].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_NEAR(-6.20, lib[1].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[1].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[1].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_NEAR(-6.20, lib[2].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[2].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[2].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_TRUE(lib[0].is_log);
    EXPECT_TRUE(lib[1].is_log);
    EXPECT_TRUE(lib[2].is_log);

    // FILE* fp = fopen("../data/context_library.lib", "w");
    // lib.Write(fp);
    // fclose(fp);
}

TEST_F(ContextLibraryTest, SortByEntropy) {
    FILE* fin = fopen("../data/K62.lib", "r");
    ContextLibrary<AA> lib(fin);
    fclose(fin);

    EXPECT_EQ((size_t)62, lib.size());
    EXPECT_EQ((size_t)1, lib.wlen());
    EXPECT_NEAR(0.0056, lib[0].probs[0][AA::kCharToInt['C']], kDelta);
    lib.SortByEntropy();
    EXPECT_NEAR(0.8196, lib[0].probs[0][AA::kCharToInt['C']], kDelta);
}

TEST_F(ContextLibraryTest, ConstructionFromSerializedLibrary) {
    FILE* fin = fopen("../data/context_library.lib", "r");
    ContextLibrary<AA> lib(fin);
    fclose(fin);

    EXPECT_EQ((size_t)3, lib.size());
    EXPECT_EQ((size_t)13, lib.wlen());
    EXPECT_NEAR(-6.20, lib[0].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[0].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[0].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_NEAR(-6.20, lib[1].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[1].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[1].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_NEAR(-6.20, lib[2].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(-4.37, lib[2].probs[2][AA::kCharToInt['N']], kDelta);
    EXPECT_NEAR( 0.00, lib[2].probs[0][AA::kCharToInt['X']], kDelta);
    EXPECT_TRUE(lib[0].is_log);
    EXPECT_TRUE(lib[1].is_log);
    EXPECT_TRUE(lib[2].is_log);
}

TEST(ContextLibraryTestInitialization, RandomSampleInitializer) {
    FILE* fin = fopen("../data/1Q7L.fas", "r");
    Alignment<AA> ali(fin, FASTA_ALIGNMENT);
    fclose(fin);
    ali.AssignMatchColumnsByGapRule();

    CountProfile<AA> p_full(ali, true);
    CountProfile<AA> p_win(p_full, 0, 13);
    std::vector< CountProfile<AA> > profiles(100, p_win);

    BlosumMatrix m;
    MatrixPseudocounts<AA> pc(m);
    ConstantAdmix admix(1.0);
    SamplingLibraryInit<AA> init(profiles, pc, admix);
    ContextLibrary<AA> lib(10, 13, init);

    EXPECT_EQ((size_t)10, lib.size());
    EXPECT_EQ((size_t)13, lib.wlen());
    EXPECT_EQ((size_t)13, lib[0].probs.length());
    EXPECT_NEAR(0.0680, lib[0].probs[0][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(0.0621, lib[0].probs[1][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(0.0976, lib[0].probs[2][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(0.1090, lib[0].probs[3][AA::kCharToInt['A']], kDelta);
    EXPECT_NEAR(0.1309, lib[0].probs[4][AA::kCharToInt['A']], kDelta);
}

TEST(ContextLibraryTestColorMapper, LearningDrosophilaPromotorContexts) {
    const size_t nshades = 5;  // shades per color channel

    FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
    ContextLibrary<Dna> lib(fin);
    fclose(fin);
    TamuraNeiMatrix sm;
    GaussianLibraryInit<Dna> init(0.1, sm, 1238087456);
    ContextLibrary<Dna> som(nshades * nshades * nshades, lib.wlen(), init);
    CoEmission<Dna> co_emission(lib.wlen(), sm);
    TransformToLin(lib);
    TransformToLin(som);

    LearnContextMap(lib, som, co_emission, 20000, 2.0, 0.1);

    AssignContextColors(lib, som, co_emission);
    AssignContextNames(lib, som, co_emission);

    ContextLibraryPdfWriter<Dna> lib_writer(lib);
    lib_writer.WriteToFile("/tmp/drosophila_promotors_K200.pdf");

    // fp = fopen("/tmp/som_learned.lib", "w");
    // som.Write(fp);
    // fclose(fp);
}

TEST(ContextLibraryTestColorMapper, LearningAbstractStates) {
    const size_t nshades = 4;  // shades per color channel
    const size_t ncolors = nshades * nshades * nshades;  // total number of colors
    BlosumMatrix sm;

    FILE* fin = fopen("../data/K62.lib", "r");
    ContextLibrary<AA> lib(fin);
    fclose(fin);
    TransformToLin(lib);
    ContextLibrary<AA> lib_pc(lib);

    GaussianLibraryInit<AA> init(0.1, sm, 3286671908);
    ContextLibrary<AA> som(ncolors, lib.wlen(), init);
    TransformToLin(som);

    CoEmission<AA> co_emission(lib.wlen(), sm, 1.0);
    LearnContextMap(lib, som, co_emission, 20000, 1.0, 0.1, 20000, 20000, 123);

    FILE* fp = fopen("../data/K62_som.lib", "w");
    som.Write(fp);
    fclose(fp);

    for (size_t k = 0; k < lib.size(); ++k)
        lib[k].name = strprintf("%c", AS62::kIntToChar[k]);

    AssignContextColors(lib, som, co_emission);
    LOG(ERROR) << lib;

    ContextLibraryPdfWriter<AA> lib_writer(lib);
    lib_writer.ncols = AS62::kSize;
    lib_writer.margin = 0.5;
    lib_writer.WriteToFile("/tmp/K62.pdf");
}

TEST(ContextLibraryTestPdfWriter, K4000) {
    FILE* fin = fopen("../data/K4000.lib", "r");
    ContextLibrary<AA> contexts(fin);
    fclose(fin);
    TransformToLin(contexts);

    ContextLibraryPdfWriter<AA> writer(contexts);
    writer.WriteToFile("/home/andreas/tmp/K4000.pdf");
}

TEST(ContextLibraryTestPdfWriter, FindContextsOfDrosophilaMotifs) {
    std::vector<Sequence<Dna> > seqs;
    FILE* fin = fopen("../data/drosophila_motifs.fas", "r");
    ReadAll(fin, seqs);
    fclose(fin);

    fin = fopen("../data/drosophila_promotors_K1000.lib", "r");
    ContextLibrary<Dna> contexts(fin);
    fclose(fin);
    TransformToLog(contexts);

    Emission<Dna> emission(contexts.wlen());
    Vector<double> pp(contexts.size());
    Vector<int> kmax(seqs.size(), 0);
    double pmax = 0.0;

    for (size_t n = 0; n < seqs.size(); ++n) {
        const size_t idx = (seqs[n].length() - 1) / 2;
        CalculatePosteriorProbs(contexts, emission, seqs[n], idx, &pp[0]);
        kmax[n] = 0;
        pmax = pp[0];
        for (size_t k = 1; k < contexts.size(); ++k) {
            if (pp[k] > pmax) {
                kmax[n] = k;
                pmax = pp[k];
            }
        }
        LOG(ERROR) << strprintf("pp[%s] = %.4f", seqs[n].header().c_str(), pmax);
        LOG(ERROR) << strprintf("\\definecolor{col%s}{rgb}{%.2f,%.2f,%.2f}",
                                seqs[n].header().c_str(),
                                contexts[kmax[n]].color.red,
                                contexts[kmax[n]].color.green,
                                contexts[kmax[n]].color.blue);
    }
    LOG(ERROR) << strprintf("\\definecolor{col%s}{rgb}{%.2f,%.2f,%.2f}",
                            "background",
                            contexts[982].color.red,
                            contexts[982].color.green,
                            contexts[982].color.blue);
    LOG(ERROR) << strprintf("\\definecolor{col%s}{rgb}{%.2f,%.2f,%.2f}",
                            "A-rich",
                            contexts[955].color.red,
                            contexts[955].color.green,
                            contexts[955].color.blue);

    TransformToLin(contexts);
    for (size_t n = 0; n < seqs.size(); ++n) {
        std::string file = strprintf("/tmp/motifs/%s.pdf", seqs[n].header().c_str());
        LOG(ERROR) << file;
        ProfilePdfWriter<Dna> writer(contexts[kmax[n]].probs);
        writer.WriteToFile(file);
    }

    ProfilePdfWriter<Dna> writer_982(contexts[982].probs);
    writer_982.WriteToFile("/tmp/motifs/background.pdf");

    ProfilePdfWriter<Dna> writer_955(contexts[955].probs);
    writer_955.WriteToFile("/tmp/motifs/A-rich.pdf");
}

}  // namespace cs
