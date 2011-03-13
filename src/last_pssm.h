// Copyright 2009, Andreas Biegert

#ifndef CS_LAST_PSSM_H_
#define CS_LAST_PSSM_H_

#include "profile-inl.h"
#include "sequence-inl.h"
#include "substitution_matrix-inl.h"

namespace cs {

// Container class for a LAST position specific scoring matrix (PSSM) and its
// associated query sequence.
class LastPssm {
  public:
    // Constructor to create a PSSM from query string and a sequence profile.
    LastPssm(const Sequence<Dna>& seq, const Profile<Dna>& prof, const SubstitutionMatrix<Dna>& mat)
            : query_(seq), pssm_(prof.length()) {

        for (size_t i = 0; i < pssm_.length(); ++i) {
            for (size_t a = 0; a < Dna::kSize; ++a) {
                pssm_[i][a] = log(prof[i][a] / mat.p(a));
                // std::cerr << prof[i][a] << "\t" << mat.p(a) << "\t" << pssm_[i][a] << std::endl;
            }
        }
    }

    // Writes PSSM in LAST format
    void Write(FILE* fout) const {
        // Print header section
        fputs("Last\n", fout);

        // Print alphabet description line
        for (size_t a = 0; a < Dna::kSize; ++a)
            fprintf(fout, "\t%c", Dna::kIntToChar[a]);
        fputs("\n", fout);

        // Print PSSM in Last's human-readable PSSM format
        for (size_t i = 0; i < pssm_.length(); ++i) {
            fprintf(fout, "%zu %c", i+1, query_.chr(i));
            for (size_t a = 0; a < Dna::kSize; ++a) {
                fprintf(fout, "\t%i", iround(pssm_[i][a]));
            }
            fputs("\n", fout);
        }
    }

  private:
    // Query sequence with which search was started
    Sequence<Dna> query_;
    // The PSSM including as log[p(a) / f(a)]
    Profile<Dna> pssm_;
};  // class LastPssm

}  // namespace cs

#endif  // CS_LAST_PSSM_H_
