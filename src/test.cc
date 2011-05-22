#include "cs.h"
#include "blosum_matrix.h"

namespace cs {


void test(int argc, char** argv) {
  /*
  if (argc != 4) {
    puts("test CRF CRF-OUT WEIGHT");
    exit(1);
  }
  FILE* fin = fopen(argv[1], "r");
  if (!fin) {
    puts("Error reading CRF!");
    exit(1);
  }
  Crf<AA> crf(fin);
  fclose(fin);
  double weight = atof(argv[3]);
  for (size_t k = 0; k < crf.size(); ++k) {
    for (size_t i = 0; i < crf.wlen(); ++i) {
      for (size_t a = 0; a < AA::kSize; ++a) {
        crf[k].context_weights[i][a] *= weight;
      }
    }
  }
  FILE* fout = fopen(argv[2], "w");
  crf.Write(fout);
  fclose(fout);
  puts("Done!");
  */
  BlosumMatrix sm(BLOSUM62);
  for (size_t a = 0; a < AA::kSize; ++a) {
    printf("%f\n", sm.r(a, 0));
  }

}


} // namespace cs


int main(int argc, char** argv) {
  cs::test(argc, argv);
}

