#include <gtest/gtest.h>

#include "cs.h"
#include "abstract_state_matrix-inl.h"
#include "blosum_matrix.h"
#include "context_library-inl.h"
#include "pdf_writer-inl.h"

namespace cs {

TEST(AbstractStatesTest, TranslateIntoAbstractStateProfile) {
    AbstractStateMatrix<AS62> matrix("../data/nr20f_151208_neff2.5_K62_r0.mat");

    FILE* fin = fopen("../data/K4000.lib", "r");
    ContextLibrary<AA> contexts(fin);
    TransformToLog(contexts);
    fclose(fin);

    Emission<AA> emission(contexts.wlen());

    fin = fopen("../data/zinc_finger.fas", "r");
    Alignment<AA> alignment(fin, FASTA_ALIGNMENT);
    fclose(fin);
    CountProfile<AA> counts(alignment);

    ProfilePdfWriter<AA> profile_writer(counts.counts);
    profile_writer.WriteToFile("/home/andreas/tmp/zinc_finger_aa.pdf");

    Profile<AS62> profile(counts.length());
    profile = TranslateIntoStateProfile<AS62>(counts, contexts, emission, matrix);

    fin = fopen("../data/K62.lib", "r");
    ContextLibrary<AA> alphabet(fin);
    fclose(fin);
    StateProfilePdfWriter<AS62, AA> as_profile_writer(profile, alphabet);
    as_profile_writer.WriteToFile("/home/andreas/tmp/zinc_finger_as.pdf");
}

TEST(AbstractStatesTest, AbstractStateMatrix) {
    AbstractStateMatrix<AS62> matrix("../data/nr20f_151208_neff2.5_K62_r2.mat");

    FILE* fin = fopen("../data/K4000.lib", "r");
    ContextLibrary<AA> contexts(fin);
    fclose(fin);

    fin = fopen("../data/K62.lib", "r");
    ContextLibrary<AA> alphabet(fin);
    fclose(fin);

    // Determine score range
    double score_max = -1000.0;
    double score_min =  1000.0;
    for (size_t k = 0; k < contexts.size(); ++k) {
        for (size_t a = 0; a < alphabet.size(); ++a) {
            score_max = MAX(score_max, matrix.s(k,a));
            score_min = MIN(score_min, matrix.s(k,a));
        }
    }
    printf("score_max=%6.2f\n", score_max);
    printf("score_min=%6.2f\n", score_min);

    // Generate PDF
    AbstractStateMatrixPdfWriter<AS62, AA> matrix_writer(matrix, contexts, alphabet);
    matrix_writer.WriteToFile("/home/andreas/tmp/nr20f_151208_neff2.5_K62_r2.pdf");
}

TEST(BlosumMatrixAbstractStatesTest, Blosum45Blosum62Blosum80) {
    const size_t nshades = 4;  // shades per color channel
    const size_t ncolors = nshades * nshades * nshades;  // total number of colors
    const double kPriorScaleFac = 0.322580645;
    BlosumMatrix sm;

    FILE* fin = fopen("../data/K62.lib", "r");
    ContextLibrary<AA> alphabet(fin);
    fclose(fin);
    TransformToLin(alphabet);

    double sum = 0.0;
    BlosumMatrix blosum45(BLOSUM45);
    for (size_t a = 0; a < AA::kSize; ++a) {
        alphabet[a].prior = kPriorScaleFac * blosum45.p(a);
        sum += alphabet[a].prior;
        for (size_t b = 0; b < AA::kSize; ++b) {
            alphabet[a].probs[0][b] = blosum45.r(b,a);
            alphabet[a].pc[b] = blosum45.r(b,a);
        }
    }

    BlosumMatrix blosum62(BLOSUM62);
    for (size_t a = 0; a < AA::kSize; ++a) {
        alphabet[AA::kSize + a].prior = kPriorScaleFac * blosum62.p(a);
        sum += alphabet[AA::kSize + a].prior;
        for (size_t b = 0; b < AA::kSize; ++b) {
            alphabet[AA::kSize + a].probs[0][b] = blosum62.r(b,a);
            alphabet[AA::kSize + a].pc[b] = blosum62.r(b,a);
        }
    }

    BlosumMatrix blosum80(BLOSUM80);
    for (size_t a = 0; a < AA::kSize; ++a) {
        alphabet[2 * AA::kSize + a].prior = kPriorScaleFac * blosum80.p(a);
        sum += alphabet[2 * AA::kSize + a].prior;
        for (size_t b = 0; b < AA::kSize; ++b) {
            alphabet[2 * AA::kSize + a].probs[0][b] = blosum80.r(b,a);
            alphabet[2 * AA::kSize + a].pc[b] = blosum80.r(b,a);
        }
    }

    alphabet[60].prior = kPriorScaleFac * 0.05;
    sum += alphabet[60].prior;
    for (size_t b = 0; b < AA::kSize; ++b) {
        alphabet[60].probs[0][b] = blosum45.p(b);
        alphabet[60].pc[b] = blosum45.p(b);
    }

    alphabet[61].prior = kPriorScaleFac * 0.05;
    sum += alphabet[61].prior;
    for (size_t b = 0; b < AA::kSize; ++b) {
        alphabet[61].probs[0][b] = blosum80.p(b);
        alphabet[61].pc[b] = blosum80.p(b);
    }

    EXPECT_NEAR(1.000, sum, 0.0001);
    alphabet.SortByEntropy();

    for (size_t k = 0; k < AS62::kSize; ++k)
      alphabet[k].name = strprintf("%c", AS62::kIntToChar[k]);

    GaussianLibraryInit<AA> init(0.1, sm, 3286671908);
    ContextLibrary<AA> som(ncolors, alphabet.wlen(), init);
    TransformToLin(som);

    CoEmission<AA> co_emission(alphabet.wlen(), sm, 1.0);
    LearnContextMap(alphabet, som, co_emission, 20000, 1.0, 0.1, 20000, 20000, 123);
    AssignContextColors(alphabet, som, co_emission);

    FILE* fout = fopen("../data/K62_blosum45-62-80.lib", "w");
    alphabet.Write(fout);
    fclose(fout);

    ContextLibraryPdfWriter<AA> alphabet_writer(alphabet);
    alphabet_writer.ncols = 20;
    alphabet_writer.WriteToFile("/home/andreas/tmp/K62_blosum45-62-80.pdf");
}


}  // namespace cs
