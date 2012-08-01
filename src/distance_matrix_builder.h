/*
  Copyright 2009 Andreas Biegert

  This file is part of the CS-BLAST package.

  The CS-BLAST package is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  The CS-BLAST package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CS_DISTANCE_MATRIX_BUILDER_H_
#define CS_DISTANCE_MATRIX_BUILDER_H_

#include "sequence-inl.h"
#include "progress_bar.h"

namespace cs {

// Abstract builder for pairwise distance matrices.
class DistanceMatrixBuilder {
  public:
    // Derived classes need to
    DistanceMatrixBuilder(size_t dim, bool verbose) : dim_(dim), verbose_(verbose) {}
    virtual ~DistanceMatrixBuilder() {}

    // Constructs the distance matrix and optionally track the progress
    Matrix<float> Build() {
        if (verbose_) puts("Computing distance matrix ...");
        Matrix<float> dist(dim_, dim_, 0.0);

        // Compute all pairwise distances
        for (size_t i = 0; i < dim_ - 1; ++i) {
            for (size_t j = i + 1; j < dim_; ++j) {
                dist[i][j] = dist[j][i] = ComputeDistance(i,j);
                if (verbose_) printf("  distance betwenn sequence %zu and %zu: %.4f\n", i, j, dist[i][j]);
            }
        }
        return dist;
    }

  protected:
    // Virtual function for calculating pairwise distance between entity i and j.
    // To be implemented by derived classes.
    virtual float ComputeDistance(size_t i, size_t j) = 0;

    size_t dim_;    // dimension of matrix
    bool verbose_;  // print detailed distance info
};


template<class Abc>
class KmerDistanceMatrixBuilder : public DistanceMatrixBuilder {
  public:
    typedef typename std::vector<Sequence<Abc> > SeqVec;
    typedef std::map<size_t, size_t> KmerTable;

    // Constructs a distance functor over given sequences. If k-mer length 0 is
    // provided we choose it as k = 1 + log_|A|(L)  as in FSA.
    KmerDistanceMatrixBuilder(const SeqVec* seqs, size_t k, bool verbose = false)
            : DistanceMatrixBuilder(seqs->size(), verbose),
              seqs_(seqs),
              kmer_length_(k),
              num_kmers_(0),
              counts_(seqs->size()) {
        const size_t nseqs = seqs_->size();

        // Choose best k-mer size based on median seq length and alphabet size?
        if (kmer_length_ == 0) {
            // Sort lengths of sequences
            std::vector<size_t> lengths;
            for (typename SeqVec::const_iterator it = seqs_->begin(); it != seqs_->end(); ++it)
                lengths.push_back(it->length());
            std::sort(lengths.begin(), lengths.end());

            // Calculate median sequence length
            size_t median = 0;
            if (nseqs % 2 == 1) median = lengths[(nseqs - 1) / 2];
            else median = 0.5 * (lengths[nseqs / 2] + lengths[nseqs / 2 - 1]);

            // This is our best k-mer length
            kmer_length_ = 1 + iround(log(median) / log(Abc::kSize));
            if (verbose_) printf("Set optimal k-mer length for distance calculation to k=%zu\n", kmer_length_);
        }
        num_kmers_ = static_cast<size_t>(pow(Abc::kSizeAny, kmer_length_));

        // Precompute powers of alphabet size
        powers_.Resize(kmer_length_);
        for (size_t k = 0; k < kmer_length_; ++k)
            powers_[k] = static_cast<size_t>(pow(Abc::kSizeAny, k));
    }

    size_t kmer_length() const { return kmer_length_; }

    virtual ~KmerDistanceMatrixBuilder() {}

  private:
    // Calculates distanace from fractional identiy between seq i and j based on
    // k-mer counts.
    virtual float ComputeDistance(size_t i, size_t j) {
        // Compute k-mer counts if not already cached
      if (counts_[i].empty()) CountKmers(seqs_->at(i), &counts_[i]);
      if (counts_[j].empty()) CountKmers(seqs_->at(j), &counts_[j]);

        // Count common k-mers
        size_t count = 0;
        if (counts_[i].size() < counts_[j].size()) {
            for (KmerTable::iterator it = counts_[i].begin(); it != counts_[i].end(); ++it) {
                KmerTable::iterator it2 = counts_[j].find(it->first);
                if (it2 != counts_[j].end()) count += MIN(it->second, it2->second);
            }
        } else {
            for (KmerTable::iterator it = counts_[j].begin(); it != counts_[j].end(); ++it) {
                KmerTable::iterator it2 = counts_[i].find(it->first);
                if (it2 != counts_[i].end()) count += MIN(it->second, it2->second);
            }
        }

        // Normalize by length of shorter sequence
        return 1 - static_cast<float>(count) / (MIN(seqs_->at(i).length(), seqs_->at(j).length()) - kmer_length_ + 1);
    }

    // Counts k-mers in sequence i and stores them in counts_ matrix.
    void CountKmers(const Sequence<Abc>& seq, KmerTable* counts) {
        for (size_t i = 0; i < seq.length() - kmer_length_; ++i) {
            // TODO: we can do this more efficiently by computing k-mer as
            // function of previous k-mer
            size_t kmer = 0;
            for (size_t k = 0; k < kmer_length_; ++k)
                kmer += seq[i + k] * powers_[kmer_length_ - k - 1];
            assert(kmer < num_kmers_);
            (*counts)[kmer]++;
        }
    }

    const SeqVec* seqs_;         // unaligned sequences to be compared
    size_t kmer_length_;         // length of k-mers
    size_t num_kmers_;           // total number of k-mers
    Vector<KmerTable> counts_;   // k-mer counts of each sequence
    Vector<size_t> powers_;      // precomputed powers of alphabet size
};



template<class Abc>
class KimuraDistanceMatrixBuilder : public DistanceMatrixBuilder {
  public:
    // Constructs a distance functor over given alignment.
    KimuraDistanceMatrixBuilder(const Alignment<Abc>* ali, bool verbose = false)
            : DistanceMatrixBuilder(ali->nseqs(), verbose),
              ali_(ali)
    { }

    virtual ~KimuraDistanceMatrixBuilder() {}

  private:
    // Calculates distanace from fractional identiy in alignment of seq i and j
    virtual float ComputeDistance(size_t i, size_t j) {
        const Alignment<Abc>& ali = *ali_;
        size_t idents = 0, matches = 0;

        for (size_t k = 0; k < ali.nmatch(); ++k) {
            if (ali[k][i] < Abc::kGap && ali[k][j] < Abc::kGap) {
                matches++;
                if (ali[k][i] == ali[k][j]) idents++;
            }
        }
        // Assume the worst if no residues were aligned
        if (matches == 0) return 10.0;
        double p = 1.0 - static_cast<float>(idents) / static_cast<float>(matches);

        // Use Kimura's formula
        if (p < 0.75)
            return -log(1 - p - (p*p) / 5);

        // Per ClustalW, return 10.0 for anything over 93% (taken from MUSCLE)
        if (p > 0.93)
            return 10.0;

        assert(p <= 1 && p >= 0.75);
        int table_index = static_cast<int>((p - 0.75)*1000 + 0.5);
        if (table_index < 0 || table_index >= kNumEntries) {
            std::cout << std::endl << *ali_ << std::endl;
            throw Exception("Error in Kimura distance calculation of sequences i=%zu and j=%zu (index=%d, p=%.2f)!", i, j, table_index, p);
        }

        return kDayhoffPams[table_index] / 100.0;
    }

    // Dayhoff Pams table
    static const int kDayhoffPams[];
    // Number of Dayhoff Entries
    static const int kNumEntries = 181;

    const Alignment<Abc>* ali_;  // aligned sequences to be compared
};


// ClustalW uses a table lookup for values > 0.75.
// The following table was copied from the ClustalW file dayhoff.h.
template<class Abc>
const int KimuraDistanceMatrixBuilder<Abc>::kDayhoffPams[] = {
    195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
    196,   /* 75.1% observed d; 196 PAMs estimated */
    197,    198,    199,    200,    200,    201,    202,  203,
    204,    205,    206,    207,    208,    209,    209,    210,    211,  212,
    213,    214,    215,    216,    217,    218,    219,    220,    221,  222,
    223,    224,    226,    227,    228,    229,    230,    231,    232,  233,
    234,    236,    237,    238,    239,    240,    241,    243,    244,  245,
    246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */
    252,    253,    254,    255,    257,  258,
    260,    261,    262,    264,    265,    267,    268,    270,    271,  273,
    274,    276,    277,    279,    281,    282,    284,    285,    287,  289,
    291,    292,    294,    296,    298,    299,    301,    303,    305,  307,
    309,    311,    313,    315,    317,    319,    321,    323,    325,  328,
    330,    332,    335,    337,    339,    342,    344,    347,    349,  352,
    354,    357,    360,    362,    365,    368,    371,    374,    377,  380,
    383,    386,    389,    393,    396,    399,    403,    407,    410,  414,
    418,    422,    426,    430,    434,    438,    442,    447,    451,  456,
    461,    466,    471,    476,    482,    487,    493,    498,    504,  511,
    517,    524,    531,    538,    545,    553,    560,    569,    577,  586,
    595,    605,    615,    626,    637,    649,    661,    675,    688,  703,
    719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
           /* 92.9% observed; 945 PAMs */
    988    /* 93.0% observed; 988 PAMs */
};

}  // namespace cs

#endif  // CS_DISTANCE_MATRIX_H_
