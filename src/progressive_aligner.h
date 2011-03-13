// Copyright 2009, Andreas Biegert

#ifndef CS_PROGRESSIVE_ALIGNER_H_
#define CS_PROGRESSIVE_ALIGNER_H_

#include "guide_tree.h"
#include "pairwise_aligner-inl.h"
#include "po_hmm-inl.h"
#include "po_hmm_merger-inl.h"
#include "distance_matrix_builder.h"

namespace cs {

struct ProgressiveAlignerParams {
    ProgressiveAlignerParams(float gf,
                             float cs,
                             float mcp,
                             float min_prob,
                             size_t max_alis,
                             size_t min_alis,
                             size_t k,
                             float pca,
                             float pcb,
                             float wcenter,
                             float b,
                             float gapo,
                             float gape,
                             size_t s)
            : gap_factor(gf),
              context_score(cs),
              min_context_prob(mcp),
              min_ali_prob(min_prob),
              max_sub_alis(max_alis),
              min_sub_alis(min_alis),
              kmer_length(k),
              pc_admix(pca),
              pc_redn(pcb),
              weight_center(wcenter),
              beta(b),
              gap_open(gapo),
              gap_extend(gape),
              max_skip_count(s)
    { }

    ProgressiveAlignerParams() { Reset(); }

    // Resets parameters to default values
    void Reset() {
        gap_factor       = 0.5;
        context_score    = 0.1;
        min_context_prob = 0.01;
        min_ali_prob     = 0.01;
        max_sub_alis     = 60;
        min_sub_alis     = 5;
        kmer_length      = 0;
        pc_admix         = 0.6;
        pc_redn          = 5.0;
        weight_center    = 1.6;
        beta             = 0.80;
        gap_open         = kGapOpen;
        gap_extend       = kGapExt;
        max_skip_count   = 0;
    }

    float gap_factor;        // gap factor in PO-HMM pairwise alignment
    float context_score;     // weight of context score in pairwise alignment
    float min_context_prob;  // minimal posterior prob for storing context in sparse profile columns
    float min_ali_prob;      // minimum pmin prob for suboptimal alignments
    size_t max_sub_alis;     // maximal number of suboptimal alignments per pairwise alignment
    size_t min_sub_alis;     // minimal number of suboptimal alignments per pairwise alignment
    size_t kmer_length;      // k-mer length in distance matrix calculation
    float pc_admix;          // maximal pseudocount admix
    float pc_redn;           // giverns pseudocount admix reduction
    float weight_center;     // central window weight
    float beta;              // context weight reduction
    float gap_open;          // gap open probability
    float gap_extend;        // gap extension probability
    size_t max_skip_count;   // maximal skip-count before vertex gets disabled
};


// A progressive aligner that first constructs a guide tree by UPGMA. This guide
// tree determines the alignment order during the progressive alignment step.
template<class Abc>
class ProgressiveAligner {
  public:
    // Typedefs
    typedef typename std::vector<Sequence<Abc> > SeqVec;
    typedef typename std::vector<Alignment<Abc>* > AliVec;
    typedef typename std::vector<POHmm<Abc>* > HmmVec;
    typedef typename std::vector<std::vector<Sequence<Abc> > > SeqSets;

    // Constructs a progressive aligner for given sequences.
    ProgressiveAligner(const ProgressiveAlignerParams* params,
                       const SubstitutionMatrix<Abc>* mat = NULL,
                       const Pseudocounts<Abc>* pc = NULL,
                       ContextLibrary<Abc>* lib = NULL,
                       bool verbose = false,
                       bool sanity = false)
            : params_(params),
              mat_(mat),
              pc_(pc),
              lib_(lib),
              verbose_(verbose),
              sanity_(sanity) {
        // Setup aligner and emission functor
        if (lib_ != NULL)
            emission_.reset(new Emission<Abc>(lib_->wlen(), params_->weight_center, params_->beta));
    }

    ~ProgressiveAligner() { }

    // Constructs a guide tree and then progressively aligns all sequences
    // starting at the leaf nodes.
    Alignment<Abc>* AlignSeqs(const SeqVec& seqs, const GuideTree& tree, std::string pdf_basename = "", bool verbose = false) {
        TreeNodeIndex index;                            // guide tree nodes
        HmmVec hmms;                                    // pointers to PO-HMMs as they are built
        scoped_ptr<TamuraNeiMatrix> matx, maty;         // substitution matrices
        scoped_ptr<PairwiseAligner<Abc> > aligner;      // PO-HMM aligner who does all the work

        // Convert initial sequences into PO-HMMs
        if (verbose) puts("Converting input sequences to PO-HMMs ...");
        for (size_t k = 0; k < seqs.size(); ++k) {
            hmms.push_back(new POHmm<Abc>(seqs[k], params_->gap_open, params_->gap_extend));

            if (mat_ != NULL)
                pc_->AddTo(CSBlastAdmix(params_->pc_admix, params_->pc_redn), hmms.back());
            if (lib_ != NULL)
                hmms[k]->AssignContextStateProbs(*lib_, *emission_, params_->min_context_prob);

            if (false && !pdf_basename.empty()) {
                Alignment<Abc> ali(seqs[k]);
                std::string filename = strprintf("%s-seq-%zu.pdf", pdf_basename.c_str(), k);
                AlignmentPdfWriter<Abc> ali_writer(ali);
                ali_writer.WriteToFile(filename);
            }
        }

        // Progressive alignment loop
        if (verbose) puts("Progressive alignment ...");
        for (size_t i = seqs.size(); i < num_vertices(tree); ++i) {
            // Determine PO-HMMs k and l to be aligned for forming PO-HMM i
            TreeAdjacencyIter ai = adjacent_vertices(i, tree).first;
            size_t k = index[*ai++];
            size_t l = index[*ai++];

            double bk = tree[*in_edges(k, tree).first].length;
            double bl = tree[*in_edges(l, tree).first].length;

            if (mat_ == NULL)  {  // phylogeny aware alignment scheme
                matx.reset(new TamuraNeiMatrix(bk));
                maty.reset(new TamuraNeiMatrix(bl));

                if (lib_ != NULL)
                    aligner.reset(new PairwiseAligner<Abc>(matx.get(), maty.get(), *lib_, params_->gap_factor, params_->context_score));
                else
                    aligner.reset(new PairwiseAligner<Abc>(matx.get(), maty.get(), params_->gap_factor));

            } else {  // co-emission based alignment scheme
                if (lib_ != NULL)
                    aligner.reset(new PairwiseAligner<Abc>(mat_, NULL, *lib_, params_->gap_factor, params_->context_score));
                else
                    aligner.reset(new PairwiseAligner<Abc>(mat_, NULL, params_->gap_factor));
            }

            POHmmMerger<Abc> merger(*hmms[k], *hmms[l], params_->gap_open, params_->gap_extend);
            AlignmentMatrices<Abc> matrices(*hmms[k], *hmms[l]);
            std::vector<PairAlignment> pairalis;

            // Calculate optimal alignment
            if (verbose) printf("  aligning nodes %zu and %zu to construct node %zu (%.4f, %.4f)\n", k, l, i, bk, bl);
            PairAlignment pairali = aligner->Align(matrices);

            merger.AddAlignment(pairali);
            pairalis.push_back(pairali);
            if (verbose) printf("    best alignment:  sum-probs = %.4f  probs-min = %.4f\n", pairali.sum_probs, pairali.pmin);

            if (false && !pdf_basename.empty()) {
                std::string filename = strprintf("%s-hmm-%zu.pdf", pdf_basename.c_str(), k);
                IndexPath path(GetPathInX(pairali));
                POHmmPdfWriter<Abc> hmm_writer(*hmms[k], lib_, &path);
                hmm_writer.WriteToFile(filename);
            }

            if (false && !pdf_basename.empty()) {
                std::string filename = strprintf("%s-alignment-%zu-0.pdf", pdf_basename.c_str(), i);
                PosteriorMatrixPdfWriter<Abc> matrix_writer(matrices, lib_, &pairali);
                matrix_writer.WriteToFile(filename);
            }

            // Calculate suboptimal alignments (except for last node)
            for (size_t n = 0; n < params_->max_sub_alis && i != num_vertices(tree) - 1; ++n) {
                pairali = aligner->Realign(pairalis.back(), matrices);
                if ((n + 1 > params_->min_sub_alis && pairali.pmin < params_->min_ali_prob) || (true && pairali.path == pairalis.back().path)) break;
                merger.AddAlignment(pairali);
                pairalis.push_back(pairali);
                if (verbose) printf("    subalignment %zu:  sum-probs = %.4f  probs-min = %.4f\n", n + 1, pairali.sum_probs, pairali.pmin);

                if (false && !pdf_basename.empty()) {
                    std::string filename = strprintf("%s-alignment-%zu-%zu.pdf", pdf_basename.c_str(), i, n + 1);
                    PosteriorMatrixPdfWriter<Abc> matrix_writer(matrices, lib_, &pairali);
                    matrix_writer.WriteToFile(filename);
                }
            }

            // Construct new PO-HMM by merging k and l based on pairwise alignments
            hmms.push_back(new POHmm<Abc>(merger.Finalize()));
            if (verbose) printf("    merged PO-HMM contains %zu sequences with Neff = %.2f\n", hmms.back()->seqs.size(), hmms.back()->neff);
            if (sanity_) hmms[i]->SanityCheck();

            if (mat_ != NULL)
                pc_->AddTo(CSBlastAdmix(params_->pc_admix, params_->pc_redn), hmms.back());
            else
                hmms.back()->AssignAncestralProbs(*hmms[k], *hmms[l], *matx, *maty);

            hmms.back()->DisableSkipped(params_->max_skip_count);

            if (lib_ != NULL)
                hmms[i]->AssignContextStateProbs(*lib_, *emission_, params_->min_context_prob);

            if (!pdf_basename.empty()) {
                std::string filename = strprintf("%s-hmm-%zu.pdf", pdf_basename.c_str(), i);
                POHmmPdfWriter<Abc> hmm_writer(*hmms[i], lib_);
                hmm_writer.hide_profiles = true;
                hmm_writer.WriteToFile(filename);
            }
        }

        // Fetch final alignment
        Alignment<Abc>* ali = new Alignment<Abc>(hmms.back()->GetAlignment());
        ali->Rearrange(seqs);

        // Free memory
        for (size_t i = 0; i < hmms.size(); ++i) delete hmms[i];

        // Print alignment to PDF if needed
        if (!pdf_basename.empty()) {
            std::string filename = strprintf("%s-alignment-final.pdf", pdf_basename.c_str());
            AlignmentPdfWriter<Abc> ali_writer(*ali);
            ali_writer.WriteToFile(filename);
        }

        return ali;
    }

    // Align the whole set of sequences
    void Align(const SeqSets& seq_sets, AliVec& io_alis, const std::vector<std::string>* pdf_basenames = NULL) {
        const int num_sets = seq_sets.size();
        Matrix<float> dist_matrix = BuildDistanceMatrix(seq_sets, io_alis);
        GuideTree tree = BuildGuideTree(dist_matrix, verbose_);
        scoped_ptr<ProgressBar> prog_bar;
        if (io_alis.empty()) io_alis.assign(num_sets, NULL);

        if (verbose_) prog_bar.reset(new ProgressBar(stdout, 70, num_sets));

        if (verbose_) printf("Progressive alignment on %d sequence sets ...\n", num_sets);
        #pragma omp parallel for schedule(static)
        for (int n = 0; n < num_sets; ++n) {
            try {
                // std::cout << pdf_basenames->at(n) << std::endl;
                if (io_alis[n] != NULL) { // free previous alignment
                    delete io_alis[n];
                    io_alis[n] = NULL;
                }
                io_alis[n] = AlignSeqs(seq_sets[n], tree, pdf_basenames != NULL ? pdf_basenames->at(n) : "", seq_sets.size() == 1);
            } catch (std::exception& e) {
                // Swallow exceptions
                io_alis[n] = NULL;
                if (num_sets == 1) throw e;
            }

            if (verbose_) {
#pragma omp critical (advance_progress)
                prog_bar->Advance(1);
            }
        }
        if (verbose_) prog_bar->Complete();
        if (verbose_) puts("");
    }

    // Build distance matrix
    Matrix<float> BuildDistanceMatrix(const SeqSets& seq_sets, const AliVec& io_alis) {
        // Calculate weight of each set
        Vector<float> weights(seq_sets.size());
        float sum_weights = 0.0;

        for (size_t n = 0; n < seq_sets.size(); ++n) {
            if (!io_alis.empty() && io_alis[n] == NULL) continue;

            const size_t nseqs = seq_sets[n].size();
            std::vector<size_t> lengths;
            for (typename SeqVec::const_iterator it = seq_sets[n].begin(); it != seq_sets[n].end(); ++it)
                lengths.push_back(it->length());
            std::sort(lengths.begin(), lengths.end());

            size_t median = 0;
            if (nseqs % 2 == 1) median = lengths[(nseqs - 1) / 2];
            else median = 0.5 * (lengths[nseqs / 2] + lengths[nseqs / 2 - 1]);

            weights[n] = median;
            sum_weights += median;
        }
        for (size_t n = 0; n < seq_sets.size(); ++n)
            weights[n] /= sum_weights;

        // Calculate distance matrix as weighted sum of all individual distances
        const size_t num_seqs = seq_sets.front().size();
        Matrix<float> dist(num_seqs, num_seqs, 0.0);

        for (size_t n = 0; n < seq_sets.size(); ++n) {
            if (!io_alis.empty() && io_alis[n] == NULL) continue;

            scoped_ptr<DistanceMatrixBuilder> builder;
            if (io_alis.empty())
                builder.reset(new KmerDistanceMatrixBuilder<Abc>(&seq_sets[n], params_->kmer_length, false));
            else
                builder.reset(new KimuraDistanceMatrixBuilder<Abc>(io_alis[n], false));

            Matrix<float> dist_n = builder->Build();
            for (size_t i = 0; i < num_seqs; ++i)
                for (size_t j = 0; j < num_seqs; ++j)
                    dist[i][j] += dist_n[i][j] * weights[n];
        }

        return dist;
    }

  private:
    // Parameter wrapper
    const ProgressiveAlignerParams* params_;
    // Substitution matrix for PCs and scoring
    const SubstitutionMatrix<Abc>* mat_;
    // Pseudocount factory for profile probs in PO-HMMs
    const Pseudocounts<Abc>* pc_;
    // Context library for context scoring
    ContextLibrary<Abc>* lib_;
    // Emission functor needed for context assignments
    scoped_ptr<Emission<Abc> > emission_;
    // Be talkative
    bool verbose_;
    // Perform expensive sanity check for each PO-HMM
    bool sanity_;
};

}  // namespace cs

#endif  // CS_PROGRESSIVE_ALIGNER_H_
