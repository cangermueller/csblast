// Copyright 2009, Andreas Biegert

#include "library_pseudocounts.h"

#include <cassert>

#include "amino_acid.h"
#include "emitter-inl.h"
#include "count_profile-inl.h"
#include "matrix.h"
#include "profile-inl.h"
#include "sequence-inl.h"
#include "utils-inl.h"

using std::valarray;

namespace cs {

template<>
void LibraryPseudocounts<AminoAcid>::add_to_sequence(
    const Sequence<AminoAcid>& seq,
    const Admixture& pca,
    Profile<AminoAcid>* profile) const {
  LOG(DEBUG1) << "Adding context-specific library pseudocounts to sequence ...";
  LOG(DEBUG1) << seq;

  if (seq.length() != profile->num_cols())
    throw Exception("Cannot add context-specific pseudocounts: "
                    "sequence and profile have different length!");

  const int length        = seq.length();
  const int num_profiles  = lib_.num_profiles();
  const int alphabet_size = AminoAcid::instance().size();
  const int center        = lib_.center();
  const float tau         = pca(1.0f);  // number of effective seqs is one

  // Profile probabilities P(p_k|X_i) at position i
  valarray<double> prob(0.0, num_profiles);
  // Pseudocount vector P(a|X_i) at position i
  valarray<double> pc(0.0, alphabet_size);
  // Output profile with pseudocounts
  Profile<AminoAcid>& p = *profile;

  for (int i = 0; i < length; ++i) {
    // Calculate profile probabilities P(p_k|X_i)
    for (int k = 0; k < num_profiles; ++k) {
      prob[k] = fast_pow2(emitter_(lib_[k], seq, i)) * lib_[k].prior();
    }
    prob /= prob.sum();  // normalization

    // Calculate pseudocount vector P(a|X_i)
    pc = 0.0;
    for (int k = 0; k < num_profiles; ++k) {
      pc[0]  += prob[k] * fast_pow2(lib_[k][center][0]);
      pc[1]  += prob[k] * fast_pow2(lib_[k][center][1]);
      pc[2]  += prob[k] * fast_pow2(lib_[k][center][2]);
      pc[3]  += prob[k] * fast_pow2(lib_[k][center][3]);
      pc[4]  += prob[k] * fast_pow2(lib_[k][center][4]);
      pc[5]  += prob[k] * fast_pow2(lib_[k][center][5]);
      pc[6]  += prob[k] * fast_pow2(lib_[k][center][6]);
      pc[7]  += prob[k] * fast_pow2(lib_[k][center][7]);
      pc[8]  += prob[k] * fast_pow2(lib_[k][center][8]);
      pc[9]  += prob[k] * fast_pow2(lib_[k][center][9]);
      pc[10] += prob[k] * fast_pow2(lib_[k][center][10]);
      pc[11] += prob[k] * fast_pow2(lib_[k][center][11]);
      pc[12] += prob[k] * fast_pow2(lib_[k][center][12]);
      pc[13] += prob[k] * fast_pow2(lib_[k][center][13]);
      pc[14] += prob[k] * fast_pow2(lib_[k][center][14]);
      pc[15] += prob[k] * fast_pow2(lib_[k][center][15]);
      pc[16] += prob[k] * fast_pow2(lib_[k][center][16]);
      pc[17] += prob[k] * fast_pow2(lib_[k][center][17]);
      pc[18] += prob[k] * fast_pow2(lib_[k][center][18]);
      pc[19] += prob[k] * fast_pow2(lib_[k][center][19]);
    }
    pc /= pc.sum();  // normalization

    // Add pseudocounts to sequence
    for(int a = 0; a < alphabet_size; ++a) {
      assert(pc[a] > 0.0);  // should be always true because of log-scaling
      const float pa = (1.0f - tau) *
        (static_cast<int>(seq[i]) == a ? 1.0f : 0.0f) + tau * pc[a];
      p[i][a] = p.logspace() ? fast_log2(pa) : pa;
    }
  }
  normalize(profile);

  LOG(DEBUG1) << *profile;
}

}  // namespace cs
