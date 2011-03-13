#include <gtest/gtest.h>

#include "cs.h"
#include "alignment-inl.h"
#include "pdf_writer-inl.h"
#include "matrix_pseudocounts-inl.h"
#include "progressive_aligner.h"
#include "sequence-inl.h"
#include "tamura_nei_matrix.h"

using std::vector;

namespace cs {

TEST(GuidedProgressiveAlignerTest, PrankAlignment) {
  TamuraNeiMatrix sm;
  MatrixPseudocounts<Dna> pc(sm);

  vector<Sequence<Dna> > seqs;
  FILE* fin = fopen("../data/iYDL154W.seq", "r");
  // FILE* fin = fopen("../data/TEL3R.seq", "r");
  ReadAll(fin, seqs);
  fclose(fin);

  // fin = fopen("../data/drosophila_promotors_K1000.lib", "r");
  // ContextLibrary<Dna> lib(fin);
  // fclose(fin);
  // TransformToLog(lib);

  ProgressiveAlignerParams params;  
  ProgressiveAligner<Dna> prog_aligner(&params, &seqs, &sm, &pc, NULL, true);

  Alignment<Dna> ali = prog_aligner.Align();
  LOG(ERROR) << ali;

  // FILE* fout = fopen("/home/andreas/tmp/csalign_test/iYDL154W.fas", "w");
  // FILE* fout = fopen("/home/andreas/tmp/csalign_test/TEL3R.fas", "w");
  // ali.Write(fout, FASTA_ALIGNMENT);
  // fclose(fout);

  // AlignmentPdfWriter<Dna> ali_writer(ali);
  // ali_writer.WriteToFile("/tmp/TEL3R.pdf");
}

}  // namespace cs
