#include <gtest/gtest.h>

#include "cs.h"
#include "blosum_matrix.h"
#include "pairwise_aligner-inl.h"
#include "alignment-inl.h"
#include "context_library-inl.h"
#include "emission.h"
#include "pdf_writer-inl.h"
#include "po_hmm-inl.h"
#include "po_hmm_merger-inl.h"
#include "sequence-inl.h"
#include "sparse_profile.h"
#include "tamura_nei_matrix.h"

namespace cs {

const double kDelta = 0.01;

// TEST(POHmmTest, ConstructionFromSequence) {
//   Sequence<AA> seq("ARNDCQEGHILKMFPSTWYV\n", "header");
//   POHmm<AA> hmm(seq);

//   EXPECT_EQ(AA::kSize, hmm.size());
// }

TEST(POHmmTest, ConstructionFromSmallAlignment) {
    FILE* fin = fopen("../data/truealn2_part.fas", "r");
    Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
    fclose(fin);
    ASSERT_EQ(0, static_cast<int>(ali.ninsert()));

    LOG(ERROR) << ali;

    POHmm<Dna> hmm(ali);

    LOG(ERROR) << hmm.neff;

    EXPECT_EQ(ali.ncols(), hmm.size());
    EXPECT_EQ(8, static_cast<int>(hmm.seqs.size()));
}

// TEST(POHmmTest, ConstructionFromAlignmentWithLongInsert) {
//   FILE* fin = fopen("../data/ali_with_long_insert.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   ASSERT_EQ(0, static_cast<int>(ali.ninsert()));

//   LOG(ERROR) << ali;

//   Matrix<double> wi;
//   Vector<double> wg;
//   Vector<double> neffs;
//   neffs = PositionSpecificWeightsAndDiversity(ali, wi);
//   double neff_global = GlobalWeightsAndDiversity(ali, wg);
//   POHmm<Dna> hmm(trans, ali, wi, wg, neffs);

//   LOG(ERROR) << neff_global;
//   LOG(ERROR) << StringifyRange(&wg[0], &wg[0] + wg.size());
//   LOG(ERROR) << hmm;

//   EXPECT_EQ(ali.ncols(), hmm.size());
//   EXPECT_EQ(2, static_cast<int>(hmm.seqs.size()));
// }

// TEST(POHmmTest, EmissionsBasedOnSequence) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);

//   Sequence<Dna> seq("TTCAGTTCGGTTT\n", "header");
//   POHmm<Dna> hmm(trans, seq);
//   Emission<Dna> emission(13);

//   // Find central vertex in POG
//   POHmm<Dna>::VertexIndex index;
//   POHmm<Dna>::VertexIter vi, vertex_end;
//   tie(vi, vertex_end) = vertices(hmm.g);
//   for (; vi != vertex_end; ++vi) if (index[*vi] == 7) break;

//   double seq_score = emission(lib[0].probs, seq, 6);
//   double hmm_score = emission(lib[0].probs, hmm.g, *vi);
//   EXPECT_NEAR( -9.715, seq_score, kDelta);
//   EXPECT_NEAR( -9.715, hmm_score, kDelta);
// }

// TEST(POHmmTest, EmissionsBasedOnAlignment) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);

//   fin = fopen("../data/context_ali.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   Matrix<double> wi;
//   Vector<double> wg;
//   Vector<double> neffs;
//   neffs = PositionSpecificWeightsAndDiversity(ali, wi);
//   GlobalWeightsAndDiversity(ali, wg);

//   POHmm<Dna> hmm(trans, ali, wi, wg, neffs);
//   Emission<Dna> emission(13);

//   ali.AssignMatchColumnsBySequence();
//   CountProfile<Dna> cp(ali);

//   // Find central vertex in POG
//   POHmm<Dna>::VertexIndex index;
//   POHmm<Dna>::VertexIter vi, vertex_end;
//   tie(vi, vertex_end) = vertices(hmm.g);
//   for (; vi != vertex_end; ++vi) if (index[*vi] == 7) break;

//   double ali_score = emission(lib[0].probs, cp, 6);
//   double hmm_score = emission(lib[0].probs, hmm.g, *vi);
//   EXPECT_NEAR( -9.715, ali_score, kDelta);
//   EXPECT_NEAR(-13.527, hmm_score, kDelta);
// }

// TEST(POHmmTest, ConstructionFromMalT) {
//   FILE* fin = fopen("../data/MalT_diverse.fas", "r");
//   Alignment<AA> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   ASSERT_EQ(0, static_cast<int>(ali.ninsert()));

//   Matrix<double> wi;
//   Vector<double> wg;
//   Vector<double> neffs;
//   neffs = PositionSpecificWeightsAndDiversity(ali, wi);
//   GlobalWeightsAndDiversity(ali, wg);
//   POHmm<AA> hmm(trans, ali, wi, wg, neffs);

//   LOG(ERROR) << hmm;

//   EXPECT_EQ(ali.ncols(), hmm.size());
//   EXPECT_EQ(21, static_cast<int>(hmm.seqs.size()));
// }

// TEST(POHmmTest, PosteriorProbs) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);
//   Emission<Dna> emission(13);

//   fin = fopen("../data/ali1000.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   ASSERT_EQ(0, static_cast<int>(ali.ninsert()));

//   Matrix<double> wi;
//   Vector<double> wg;
//   Vector<double> neffs;
//   neffs = PositionSpecificWeightsAndDiversity(ali, wi);
//   GlobalWeightsAndDiversity(ali, wg);
//   POHmm<Dna> hmm(trans, ali, wi, wg, neffs);

//   ASSERT_EQ(ali.ncols(), hmm.size());
//   ASSERT_EQ(9, static_cast<int>(hmm.seqs.size()));

//   POHmm<Dna>::VertexIndex index;
//   POHmm<Dna>::VertexIter vi, vertex_end;
//   Matrix<double> pp(hmm.size(), lib.size());
//   for (tie(vi, vertex_end) = vertices(hmm.g); vi != vertex_end; ++vi) {
//     if (index[*vi] != kStartEndVertex) {
//       CalculatePosteriorProbs(lib, emission, hmm.g, *vi, &pp[index[*vi] - 1][0]);
//     }
//   }

//   EXPECT_NEAR(0.8771, pp[155][39], kDelta);
//   EXPECT_NEAR(0.5667, pp[52][115], kDelta);
// }

// TEST(POHmmTest, AbstractStateProbs) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);
//   Emission<Dna> emission(13);

//   fin = fopen("../data/ali1000.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   ASSERT_EQ(0, static_cast<int>(ali.ninsert()));

//   Matrix<double> wi;
//   Vector<double> wg;
//   Vector<double> neffs;
//   neffs = PositionSpecificWeightsAndDiversity(ali, wi);
//   GlobalWeightsAndDiversity(ali, wg);
//   POHmm<Dna> hmm(trans, ali, wi, wg, neffs);

//   ASSERT_EQ(ali.ncols(), hmm.size());
//   ASSERT_EQ(9, static_cast<int>(hmm.seqs.size()));

//   hmm.AssignContextStateProbs(lib, emission);

//   EXPECT_NEAR(0.8771, hmm.g[156].contexts[ 39], kDelta);
//   EXPECT_NEAR(0.5667, hmm.g[ 53].contexts[115], kDelta);
// }

// TEST(POHmmTest, ForwardBackwardSeqVersusSeq) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);

//   Sequence<Dna> seq_x("CGCGCGCGCGAAAGTCCGGACTCAACAAC\n", "seq1");
//   Sequence<Dna> seq_y("CGAAAGTCCGGACTCAACAACAAC\n", "seq2");
//   POHmm<Dna> hmm_x(trans, seq_x);
//   POHmm<Dna> hmm_y(trans, seq_y);

//   LOG(ERROR) << hmm_x;
//   LOG(ERROR) << hmm_y;

//   Emission<Dna> emission(13);
//   hmm_x.AssignContextStateProbs(lib, emission, 0.01);
//   hmm_y.AssignContextStateProbs(lib, emission, 0.01);

//   AlignmentMatrices<Dna> matrices(hmm_x, hmm_y);
//   PairwiseAligner<Dna> aligner(lib);
//   PairAlignment pairali1 = aligner.Align(matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali2 = aligner.Realign(pairali1, matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali3 = aligner.Realign(pairali2, matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali4 = aligner.Realign(pairali3, matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali5 = aligner.Realign(pairali4, matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali6 = aligner.Realign(pairali5, matrices);
//   LOG(ERROR) << matrices;
//   PairAlignment pairali7 = aligner.Realign(pairali6, matrices);
//   LOG(ERROR) << matrices;

//   LOG(ERROR) << pairali1;
//   LOG(ERROR) << pairali2;
//   LOG(ERROR) << pairali3;
//   LOG(ERROR) << pairali4;
//   LOG(ERROR) << pairali5;
//   LOG(ERROR) << pairali6;
//   LOG(ERROR) << pairali7;
// }

// TEST(POHmmTest, ForwardBackwardAliVersusSeq) {
//   FILE* fin = fopen("../data/drosophila_promotors_K1000.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);

//   fin = fopen("../data/bugfix_alignment.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);
//   ASSERT_EQ(0, static_cast<int>(ali.ninsert()));
//   POHmm<Dna> hmm_x(ali);

//   fin = fopen("../data/bugfix_seq.seq", "r");
//   Sequence<Dna> seq(fin);
//   fclose(fin);
//   POHmm<Dna> hmm_y(seq);

//   LOG(ERROR) << hmm_x;
//   LOG(ERROR) << hmm_y;

//   Emission<Dna> emission(13);
//   hmm_x.AssignContextStateProbs(lib, emission, 0.01);
//   hmm_y.AssignContextStateProbs(lib, emission, 0.01);
//   TransformToLin(lib);

//   AlignmentMatrices<Dna> matrices(hmm_x, hmm_y);
//   PairwiseAligner<Dna> aligner(lib);

//   PairAlignment pairali1 = aligner.Align(matrices);

//   PosteriorMatrixPdfWriter<Dna> matrix_writer(matrices, NULL, &pairali1);
//   matrix_writer.WriteToFile("/tmp/posterior_matrix.pdf");
//   IndexPath path = GetPathInX(pairali1);
//   POHmmPdfWriter<Dna> hmm_writer(hmm_x, NULL, &path);
//   hmm_writer.WriteToFile("/tmp/po_hmm_x.pdf");

//   LOG(ERROR) << matrices;
//   // PairAlignment pairali2 = aligner.Realign(pairali1, matrices);
//   // LOG(ERROR) << matrices;
//   // PairAlignment pairali3 = aligner.Realign(pairali2, matrices);
//   // LOG(ERROR) << matrices;
//   // PairAlignment pairali4 = aligner.Realign(pairali3, matrices);
//   // LOG(ERROR) << matrices;

//   LOG(ERROR) << pairali1;
//   // LOG(ERROR) << pairali2;
//   // LOG(ERROR) << pairali3;
//   // LOG(ERROR) << pairali4;

//   POHmmMerger<Dna> merger(hmm_x, hmm_y);
//   merger.AddAlignment(pairali1);
//   // merger.AddAlignment(pairali2);
//   // merger.AddAlignment(pairali3);
//   // merger.AddAlignment(pairali4);
//   POHmm<Dna> hmm_z = merger.Finalize();
//   LOG(ERROR) << hmm_z;

//   LOG(ERROR) << seq;
//   LOG(ERROR) << ali;

//   LOG(ERROR) << hmm_z.GetAlignment(0);
//   // LOG(ERROR) << hmm_z.GetAlignment(1);
//   // LOG(ERROR) << hmm_z.GetAlignment(2);
//   // LOG(ERROR) << hmm_z.GetAlignment(3);

//   // PairAlignment pairali3 = pairali1;
//   // pairali3.path.clear();
//   // pairali3.path.push_back(Step( 1,  1, MM, 1.0));
//   // pairali3.path.push_back(Step( 2,  2, MM, 1.0));
//   // pairali3.path.push_back(Step( 3,  3, MM, 1.0));
//   // pairali3.path.push_back(Step( 4,  4, MM, 1.0));
//   // pairali3.path.push_back(Step( 8,  5, MM, 1.0));
//   // pairali3.path.push_back(Step( 9,  6, MM, 1.0));
//   // pairali3.path.push_back(Step(10,  7, MM, 1.0));
//   // pairali3.path.push_back(Step(11,  8, MM, 1.0));
//   // pairali3.path.push_back(Step(12,  9, MM, 1.0));
//   // pairali3.path.push_back(Step(13, 10, MM, 1.0));
//   // pairali3.path.push_back(Step(14, 11, MM, 1.0));

//   // POHmmMerger<Dna> merger2(hmm_z, hmm_y);
//   // merger2.AddAlignment(pairali3);
//   // POHmm<Dna> hmm_q = merger2.Finalize();
//   // LOG(ERROR) << hmm_q;
//   // LOG(ERROR) << hmm_q.GetAlignment(0);
// }

// TEST(POHmmTest, MergeSeqHmmWithSeqHmm) {
//   FILE* fin = fopen("../data/drosophila_promotors_K200.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   TransformToLog(lib);
//   fclose(fin);

//   Sequence<Dna> seq_x("TCCGAAAGT  CCGGAC\n", "seq1");
//   Sequence<Dna> seq_y("TC GAAAGTGTCCGGAC\n", "seq2");
//   POHmm<Dna> hmm_x(seq_x);
//   POHmm<Dna> hmm_y(seq_y);

//   Emission<Dna> emission(13);
//   hmm_x.AssignContextStateProbs(lib, emission, 0.01);
//   hmm_y.AssignContextStateProbs(lib, emission, 0.01);

//   AlignmentMatrices<Dna> matrices(hmm_x, hmm_y);
//   PairwiseAligner<Dna> aligner(lib);

//   PairAlignment pairali1 = aligner.Align(matrices);
//   LOG(INFO) << matrices;
//   PairAlignment pairali2 = aligner.Realign(pairali1, matrices);
//   LOG(INFO) << matrices;
//   PairAlignment pairali3 = aligner.Realign(pairali2, matrices);
//   LOG(INFO) << matrices;
//   PairAlignment pairali4 = aligner.Realign(pairali3, matrices);
//   LOG(INFO) << matrices;

//   LOG(INFO) << pairali1;
//   LOG(INFO) << pairali2;
//   LOG(INFO) << pairali3;
//   LOG(INFO) << pairali4;

//   POHmmMerger<Dna> merger(hmm_x, hmm_y);

//   merger.AddAlignment(pairali1);
//   merger.AddAlignment(pairali2);
//   merger.AddAlignment(pairali3);
//   merger.AddAlignment(pairali4);

//   POHmm<Dna> hmm_z = merger.Finalize();
//   LOG(INFO) << hmm_z;

//   LOG(INFO) << hmm_z.GetAlignment(0);
//   LOG(INFO) << hmm_z.GetAlignment(1);
//   LOG(INFO) << hmm_z.GetAlignment(2);
//   LOG(INFO) << hmm_z.GetAlignment(3);
// }

// TEST(POHmmTest, SimpleProgressiveAlignment) {
//   static size_t kNumAlis = 3;
//   FILE* fin = fopen("../data/drosophila_promotors_K1000.lib", "r");
//   ContextLibrary<Dna> lib(fin);
//   fclose(fin);

//   fin = fopen("../data/FBgn0029167_mult.fas", "r");
//   //fin = fopen("../data/ali2566_part.fas", "r");
//   Alignment<Dna> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);

//   LOG(ERROR) << ali;

//   Emission<Dna> emission(lib.wlen());
//   POHmm<Dna> hmm_x(ali.GetSequence(0));
//   PairwiseAligner<Dna> aligner(lib, 0.1);
//   std::vector<PairAlignment> alis;

//   LOG(ERROR) << hmm_x;

//   // hmm_x.AssignContextStateProbs(lib, emission, 0.0001);
//   // TransformToLin(lib);
//   // POHmmPdfWriter<Dna> hmm_x_writer(hmm_x, &lib);
//   // hmm_x_writer.WriteToFile("/tmp/hmm_x.pdf");
//   TransformToLog(lib);

//   for (size_t n = 1; n < ali.nseqs(); ++n) {
//     POHmm<Dna> hmm_y(ali.GetSequence(n));
//     LOG(ERROR) << std::string(100, '-');
//     LOG(ERROR) << "Aligning with " << hmm_y.seqs.front().header();

//     hmm_x.AssignContextStateProbs(lib, emission, 0.0001);
//     hmm_y.AssignContextStateProbs(lib, emission, 0.0001);

//     // POHmm<Dna> hmm_x_bak(hmm_x);
//     // if (n == 2) {
//     //   hmm_x = hmm_y;
//     //   hmm_y = hmm_x_bak;
//     // }

//     POHmmMerger<Dna> merger(hmm_x, hmm_y);
//     AlignmentMatrices<Dna> matrices(hmm_x, hmm_y);
//     PairAlignment pairali = aligner.Align(matrices);

//     //std::string basename = strprintf("/tmp/ali2566_part_%zu", n);
//     std::string basename = strprintf("/tmp/eve_stripe_%zu", n);
//     // if (false && n >= 0) {
//     //   IndexPath path_x(GetPathInX(pairali));
//     //   POHmmPdfWriter<Dna> hmm_writer(hmm_x, &lib, &path_x);
//     //   hmm_writer.WriteToFile(basename + "_0_hmm.pdf");

//     //   PosteriorMatrixPdfWriter<Dna> matrix_writer(matrices, &lib, &pairali);
//     //   matrix_writer.WriteToFile(basename + "_0.pdf");
//     // }

//     merger.AddAlignment(pairali);
//     alis.clear();
//     alis.push_back(pairali);

//     for (size_t a = 0; a < kNumAlis; ++a) {
//       PairAlignment subpairali = aligner.Realign(alis.back(), matrices);
//       if (subpairali.pmin < 0.05) break;

//       if (false && n >= 0) {
//         PosteriorMatrixPdfWriter<Dna> matrix_writer(matrices, &lib, &subpairali);
//         matrix_writer.WriteToFile(basename + strprintf("_%zu.pdf", a + 1));
//       }

//       merger.AddAlignment(subpairali);
//       alis.push_back(subpairali);
//     }

//     hmm_x = merger.Finalize();
//     hmm_x.AssignContextStateProbs(lib, emission, 0.0001);
//     LOG(ERROR) << hmm_x;

//     if (false && n >= 1) {
//       POHmmPdfWriter<Dna> hmm_writer(hmm_x, &lib);
//       hmm_writer.WriteToFile(basename + "_final_hmm.pdf");
//     }

//     hmm_x.SanityCheck();

//     for (size_t a = 0; a < hmm_x.num_alis; ++a) {
//       Alignment<Dna> poa(hmm_x.GetAlignment(a));
//       LOG(DEBUG) << "aliw=" << hmm_x.aliw[a];
//       LOG(DEBUG) << poa;
//       for (size_t n = 0; n < hmm_x.seqs.size(); ++n) {
//         Sequence<Dna> seq_ali(ali.GetSequence(n));
//         Sequence<Dna> seq_poa(poa.GetSequence(n));
//         LOG(DEBUG) << seq_ali;
//         LOG(DEBUG) << seq_poa;
//         assert_eq(seq_ali.length(), seq_poa.length());
//       }
//     }
//   }
//   Alignment<Dna> final_ali(hmm_x.GetAlignment(0));
//   AlignmentPdfWriter<Dna> final_ali_writer(final_ali);
//   final_ali_writer.WriteToFile("/tmp/FBgn0029167_final_alignment_gapf0.1.pdf");

//   AlignmentPdfWriter<Dna> ref_ali_writer(ali);
//   ref_ali_writer.WriteToFile("/tmp/FBgn0029167_ref_alignment.pdf");

//   LOG(ERROR) << ali;
//   LOG(ERROR) << hmm_x;
// }

// TEST(DummyTest, GarfielSeq1dAlignmentPdfWriter) {
//   FILE* fin = fopen("../data/garfield_seq1.seq", "r");
//   Alignment<AA> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);

//   AlignmentPdfWriter<AA> ali_writer(ali);
//   ali_writer.WriteToFile("/tmp/garfield_seq1.pdf");
// }

// TEST(DummyTest, GarfielSeq2dAlignmentPdfWriter) {
//   FILE* fin = fopen("../data/garfield_seq2.seq", "r");
//   Alignment<AA> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);

//   AlignmentPdfWriter<AA> ali_writer(ali);
//   ali_writer.WriteToFile("/tmp/garfield_seq2.pdf");
// }

// TEST(DummyTest, GarfielSeq3dAlignmentPdfWriter) {
//   FILE* fin = fopen("../data/garfield_seq3.seq", "r");
//   Alignment<AA> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);

//   AlignmentPdfWriter<AA> ali_writer(ali);
//   ali_writer.WriteToFile("/tmp/garfield_seq3.pdf");
// }

// TEST(DummyTest, GarfielSeq4dAlignmentPdfWriter) {
//   FILE* fin = fopen("../data/garfield_seq4.seq", "r");
//   Alignment<AA> ali(fin, FASTA_ALIGNMENT);
//   fclose(fin);

//   AlignmentPdfWriter<AA> ali_writer(ali);
//   ali_writer.WriteToFile("/tmp/garfield_seq4.pdf");
// }

TEST(DummyTest, DummyTestPdfWriter) {
    FILE* fin = fopen("../data/zinc_finger_hhblits.fas", "r");
    Alignment<AA> ali(fin, FASTA_ALIGNMENT);
    fclose(fin);

    AlignmentPdfWriter<AA> ali_writer(ali);
    ali_writer.WriteToFile("/tmp/zinc_finger.pdf");

    CountProfile<AA> profile(ali, true);
    Normalize(profile.counts, 1.0);

    BlosumMatrix sm;
    Sequence<AA> cons(ConsensusSequence(profile, sm));
    Alignment<AA> cons_ali(cons);

    AlignmentPdfWriter<AA> cons_writer(cons_ali);
    cons_writer.WriteToFile("/tmp/zinc_finger_cons.pdf");

    ProfilePdfWriter<AA> profile_writer(profile.counts);
    profile_writer.WriteToFile("/tmp/zinc_finger_profile.pdf");
}

TEST(DummyTestEve, DummyTestEvePdfWriter) {
    FILE* fin = fopen("/home/andreas/tmp/eve0.fas", "r");
    Alignment<Dna> eve0(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<Dna> eve0_writer(eve0);
    eve0_writer.WriteToFile("/home/andreas/tmp/eve0.pdf");

    fin = fopen("/home/andreas/tmp/eve1.fas", "r");
    Alignment<Dna> eve1(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<Dna> eve1_writer(eve1);
    eve1_writer.WriteToFile("/home/andreas/tmp/eve1.pdf");

    fin = fopen("/home/andreas/tmp/eve2.fas", "r");
    Alignment<Dna> eve2(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<Dna> eve2_writer(eve2);
    eve2_writer.WriteToFile("/home/andreas/tmp/eve2.pdf");

    fin = fopen("/home/andreas/tmp/eve3.fas", "r");
    Alignment<Dna> eve3(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<Dna> eve3_writer(eve3);
    eve3_writer.WriteToFile("/home/andreas/tmp/eve3.pdf");

    fin = fopen("/home/andreas/tmp/eve4.fas", "r");
    Alignment<Dna> eve4(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<Dna> eve4_writer(eve4);
    eve4_writer.WriteToFile("/home/andreas/tmp/eve4.pdf");
}


TEST(DummyTestGarfield, DummyTestGarfieldPdfWriter) {
    FILE* fin = fopen("../data/garfield_seq1_seq2_corr.fas", "r");
    Alignment<AA> ali0(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali0_writer(ali0);
    ali0_writer.WriteToFile("/home/andreas/tmp/garfield_seq1_seq2_corr.pdf");

    fin = fopen("../data/garfield_seq1_seq2.fas", "r");
    Alignment<AA> ali1(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali1_writer(ali1);
    ali1_writer.WriteToFile("/home/andreas/tmp/garfield_seq1_seq2.pdf");

    fin = fopen("../data/garfield_seq1_seq3.fas", "r");
    Alignment<AA> ali2(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali2_writer(ali2);
    ali2_writer.WriteToFile("/home/andreas/tmp/garfield_seq1_seq3.pdf");

    fin = fopen("../data/garfield_seq1_seq4.fas", "r");
    Alignment<AA> ali3(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali3_writer(ali3);
    ali3_writer.WriteToFile("/home/andreas/tmp/garfield_seq1_seq4.pdf");

    fin = fopen("../data/garfield_seq3_seq2.fas", "r");
    Alignment<AA> ali4(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali4_writer(ali4);
    ali4_writer.WriteToFile("/home/andreas/tmp/garfield_seq3_seq2.pdf");

    fin = fopen("../data/garfield_seq4_seq2.fas", "r");
    Alignment<AA> ali5(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali5_writer(ali5);
    ali5_writer.WriteToFile("/home/andreas/tmp/garfield_seq4_seq2.pdf");
}

TEST(DummyTestGarfield, DummyZincWriter) {
    FILE* fin = fopen("../data/zinc_finger_cropped.fas", "r");
    Alignment<AA> ali0(fin, FASTA_ALIGNMENT);
    fclose(fin);
    AlignmentPdfWriter<AA> ali0_writer(ali0);
    ali0_writer.WriteToFile("/home/andreas/Desktop/zinc_finger.pdf");
}

}  // namespace cs
