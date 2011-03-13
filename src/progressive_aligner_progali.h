// Copyright 2009, Andreas Biegert

#ifndef CS_PROGRESSIVE_ALIGNER_H_
#define CS_PROGRESSIVE_ALIGNER_H_

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>

#include "pairwise_aligner-inl.h"
#include "po_hmm-inl.h"
#include "po_hmm_merger-inl.h"
#include "distance_matrix_builder.h"

using boost::adjacency_list;
using boost::vecS;
using boost::bidirectionalS;
using boost::graph_traits;
using boost::property_map;
using boost::vertex_index_t;
using boost::vertex_index;

namespace cs {

struct ProgressiveAlignerParams {
    ProgressiveAlignerParams(float gf,
                             float cs,
                             float mcp,
                             float map,
                             size_t msa,
                             size_t k,
                             float pca,
                             float pcb,
                             float wcenter,
                             float b,
                             float gapo,
                             float gape)
            : gap_factor(gf),
              context_score(cs),
              min_context_prob(mcp),
              min_ali_prob(map),
              max_sub_alis(msa),
              kmer_length(k),
              pc_admix(pca),
              pc_redn(pcb),
              weight_center(wcenter),
              beta(b),
              gap_open(gapo),
              gap_extend(gape)
    { }

    ProgressiveAlignerParams() { Reset(); }

    // Resets parameters to default values
    void Reset() {
        gap_factor       = 0.5;
        context_score    = 0.1;
        min_context_prob = 0.01;
        min_ali_prob     = 0.001;
        max_sub_alis     = 60;
        kmer_length      = 0;
        pc_admix         = 0.5;
        pc_redn          = 8.0;
        weight_center    = 5.0;
        beta             = 0.85;
        gap_open         = kGapOpen;
        gap_extend       = kGapExt;
    }

    float gap_factor;        // gap factor in PO-HMM pairwise alignment
    float context_score;     // weight of context score in pairwise alignment
    float min_context_prob;  // minimal posterior prob for storing context in sparse profile columns
    float min_ali_prob;      // minimum pmin prob for suboptimal alignments
    size_t max_sub_alis;     // maximal number of suboptimal alignments per pairwise alignment
    size_t kmer_length;      // k-mer length in distance matrix calculation
    float pc_admix;          // maximal pseudocount admix
    float pc_redn;           // giverns pseudocount admix reduction
    float weight_center;     // central window weight
    float beta;              // context weight reduction
    float gap_open;          // gap open probability
    float gap_extend;        // gap extension probability
};


// POD class representing a tree node.
struct TreeNodeProperties {
    // Constructs an empty node
    TreeNodeProperties(size_t d = 0, size_t n = 1, double t = 0.0, double h = 0.0)
            : depth(d), size(n), dist(t), height(h) {}

    size_t depth;   // depth of leafs is 0
    size_t size;    // number of seqs in this cluster
    double dist;    // distance time of child nodes
    double height;  // decimal height of this node in tree
};

// POD class representing a tree edge.
struct TreeEdgeProperties {
    // Constructs an edge with given length
    TreeEdgeProperties(double l = 0.0) : length(l) {}

    double length;  // edge length
};


// A progressive aligner that first constructs a guide tree by UPGMA. This guide
// tree determines the alignment order during the progressive alignment step.
// TODO: add some sort of parameter wrapper for progressive aligner
template<class Abc>
class ProgressiveAligner {
  public:
    // Typedefs
    typedef adjacency_list<vecS, vecS, bidirectionalS, TreeNodeProperties, TreeEdgeProperties> Tree;
    typedef graph_traits<Tree> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor Node;
    typedef typename property_map<Tree, vertex_index_t>::type NodeIndex;
    typedef typename GraphTraits::adjacency_iterator AdjacencyIter;
    typedef typename std::vector<Sequence<Abc> > SeqVec;
    typedef typename std::vector<POHmm<Abc>* > HmmVec;

    // Constructs a progressive aligner for given sequences. The context library
    // is optional. If it is provided, the aligner will use context scoring instead
    // of standard substitution matrix based alignment.
    ProgressiveAligner(const ProgressiveAlignerParams* params,
                       const SeqVec* seqs,
                       const SubstitutionMatrix<Abc>* mat,
                       const Pseudocounts<Abc>* pc,
                       ContextLibrary<Abc>* lib,
                       bool verbose = false,
                       bool sanity = false,
                       std::string pdf_basename = "")
            : params_(params),
              seqs_(seqs),
              mat_(mat),
              pc_(pc),
              lib_(lib),
              verbose_(verbose),
              sanity_(sanity),
              pdf_basename_(pdf_basename) {
        // Setup aligner and emission
        if (lib_ != NULL) emission_.reset(new Emission<Abc>(lib_->wlen(), params_->weight_center, params_->beta));
    }

    virtual ~ProgressiveAligner() { Reset(); }

    // Constructs a guide tree and then progressively aligns all sequences
    // starting at the leaef nodes.
    Alignment<Abc> Align() {
        Reset();
        Tree tree = BuildGuideTree();
        NodeIndex index;

        if (verbose_) puts("Progressive alignment ...");
        // Convert initial sequences into PO-HMMs
        for (size_t k = 0; k < seqs_->size(); ++k) {
            hmms_.push_back(new POHmm<Abc>((*seqs_)[k], params_->gap_open, params_->gap_extend));
            pc_->AddTo(CSBlastAdmix(params_->pc_admix, params_->pc_redn), hmms_.back());
            // if (lib_ != NULL) hmms_[k]->AssignContextStateProbs(*lib_, *emission_, params_->min_context_prob);

            if (!pdf_basename_.empty()) {
                Alignment<Abc> ali((*seqs_)[k]);
                std::string filename = strprintf("%s-seq-%zu.pdf", pdf_basename_.c_str(), k);
                AlignmentPdfWriter<Abc> ali_writer(ali);
                ali_writer.WriteToFile(filename);
            }
        }

        for (size_t i = seqs_->size(); i < num_vertices(tree); ++i) {

            // Determine PO-HMMs k and l to be aligned for forming PO-HMM i
            AdjacencyIter ai = adjacent_vertices(i, tree).first;
            size_t k = index[*ai++];
            size_t l = index[*ai++];

            size_t tmp = k;
            k = l;
            l = tmp;

            double bk = tree[*in_edges(k, tree).first].length;
            double bl = tree[*in_edges(l, tree).first].length;

            // TamuraNeiMatrix matx(bk);
            // TamuraNeiMatrix maty(bl);
            // if (lib_ != NULL)
            //    aligner_.reset(new PairwiseAligner<Abc>(&matx, &maty, *lib_, params_->gap_factor, params_->context_score));
            // else
            //    aligner_.reset(new PairwiseAligner<Abc>(&matx, &maty, params_->gap_factor));

            if (lib_ != NULL)
                aligner_.reset(new PairwiseAligner<Abc>(mat_, mat_, *lib_, params_->gap_factor, params_->context_score));
            else
                aligner_.reset(new PairwiseAligner<Abc>(mat_, mat_, (i != 12 || false) ? params_->gap_factor : 0.5));

            POHmmMerger<Abc> merger(*hmms_[k], *hmms_[l], params_->gap_open, params_->gap_extend);
            AlignmentMatrices<Abc> matrices(*hmms_[k], *hmms_[l]);
            std::vector<PairAlignment> pairalis;

            // Calculate optimal alignment
            if (verbose_) printf("  aligning nodes %zu and %zu to construct node %zu (%.4f, %.4f)\n", k, l, i, bk, bl);
            PairAlignment pairali = aligner_->Align(matrices);

            merger.AddAlignment(pairali);
            pairalis.push_back(pairali);
            if (verbose_) printf("    best alignment:  probs-average = %.4f  probs-min = %.4f\n", pairali.avg_probs, pairali.pmin);

            // if (!pdf_basename_.empty() && i == 14 && false) {
            //     std::string filename = strprintf("%s-hmm-%zu.pdf", pdf_basename_.c_str(), k);
            //     IndexPath path(GetPathInX(pairali));
            //     POHmmPdfWriter<Abc> hmm_writer(*hmms_[k], lib_, &path);
            //     hmm_writer.WriteToFile(filename);
            // }

            if (!pdf_basename_.empty()) {
                std::string filename = strprintf("%s-alignment-%zu-0.pdf", pdf_basename_.c_str(), i);
                PosteriorMatrixPdfWriter<Abc> matrix_writer(matrices, lib_, &pairali);
                matrix_writer.WriteToFile(filename);
            }

            // Calculate suboptimal alignments except for last node
            for (size_t n = 0; n < params_->max_sub_alis && i != num_vertices(tree) - 1 && i < 5; ++n) {
                pairali = aligner_->Realign(pairalis.back(), matrices);
                if (pairali.pmin < params_->min_ali_prob) break;
                merger.AddAlignment(pairali);
                pairalis.push_back(pairali);
                if (verbose_) printf("    subalignment %zu:  probs-average = %.4f  probs-min = %.4f\n", n + 1, pairali.avg_probs, pairali.pmin);

                if (!pdf_basename_.empty()) {
                    std::string filename = strprintf("%s-alignment-%zu-%zu.pdf", pdf_basename_.c_str(), i, n + 1);
                    PosteriorMatrixPdfWriter<Abc> matrix_writer(matrices, lib_, &pairali);
                    matrix_writer.WriteToFile(filename);
                }
            }

            // Construct new PO-HMM by merging k and l based on pairwise alignments
            hmms_.push_back(new POHmm<Abc>(merger.Finalize()));
            LOG(ERROR) << *hmms_.back();
            if (verbose_) printf("    merged PO-HMM contains %zu sequences with Neff = %.2f\n", hmms_.back()->seqs.size(), hmms_.back()->neff);
            if (sanity_) hmms_[i]->SanityCheck();

            pc_->AddTo(CSBlastAdmix(params_->pc_admix, params_->pc_redn), hmms_.back());
            // hmms_.back()->AssignAncestralProbs(*hmms_[k], *hmms_[l], *mat_, *mat_);
            // hmms_.back()->AssignAncestralProbs(*hmms_[k], *hmms_[l], matx, maty);
            if (lib_ != NULL) hmms_[i]->AssignContextStateProbs(*lib_, *emission_, params_->min_context_prob);


            if (!pdf_basename_.empty()) {
                std::string filename = strprintf("%s-hmm-%zu.pdf", pdf_basename_.c_str(), i);
                POHmmPdfWriter<Abc> hmm_writer(*hmms_[i], lib_);
                hmm_writer.hide_profiles = true;
                hmm_writer.WriteToFile(filename);
            }
        }

        Alignment<Abc> ali = hmms_.back()->GetAlignment();
        ali.Rearrange(*seqs_);

        if (!pdf_basename_.empty()) {
            std::string filename = strprintf("%s-alignment-final.pdf", pdf_basename_.c_str());
            AlignmentPdfWriter<Abc> ali_writer(ali);
            ali_writer.WriteToFile(filename);
        }

        return ali;
    }

    // Constructs a guide tree by UPGMA using pairwise MAC scores as distances.
    // Tree BuildGuideTree() {
    //     KmerDistanceMatrixBuilder<Abc> dist_builder(seqs_, params_->kmer_length, verbose_);
    //     Matrix<float> dist = dist_builder.Build();
    //     Tree tree;
    //     Vector<Node> nodes(seqs_->size());  // nodes[i] is node of i-th cluster
    //     Bitset disabled(seqs_->size());     // indicates disabled indices

    //     if (verbose_) puts("Building UPGMA guide tree ...");

    //     // Add leaf nodes to guide tree
    //     for (size_t k = 0; k < seqs_->size(); ++k)
    //         nodes[k] = add_vertex(TreeNodeProperties(), tree);

    //     // Using dist[k][l] find the closest 2 clusters, amalgamate, and update
    //     // dist[i][j]'s. Repeatedly do this, until there is only 1 cluster left.
    //     size_t num_clusters = seqs_->size();
    //     while (num_clusters > 1) {
    //         // Find clusters k and l with maximal distance
    //         size_t k = 0, l = 0;
    //         float dmin = 1.0;
    //         for (size_t i = 0; i < seqs_->size() - 1; ++i) {
    //             for (size_t j = i + 1; j < seqs_->size(); ++j) {
    //                 if (!disabled.test(i) && !disabled.test(j) && dist[i][j] < dmin) {
    //                     dmin = dist[i][j];
    //                     k = i;
    //                     l = j;
    //                 }
    //             }
    //         }

    //         // Create new node p with child nodes c[k] and c[l]
    //         size_t depth = 1 + MAX(tree[nodes[k]].depth, tree[nodes[l]].depth);
    //         double height = dmin / 2.0;
    //         size_t sk = tree[nodes[k]].size;
    //         size_t sl = tree[nodes[l]].size;
    //         double bk = height - tree[nodes[k]].height;
    //         double bl = height - tree[nodes[l]].height;
    //         Node p = add_vertex(TreeNodeProperties(depth, sk + sl, dmin, height), tree);
    //         add_edge(p, nodes[k], TreeEdgeProperties(bk), tree);
    //         add_edge(p, nodes[l], TreeEdgeProperties(bl), tree);
    //         nodes[k] = p;
    //         disabled.set(l);
    //         num_clusters--;

    //         // Update distance matrix
    //         float size = sk + sl;
    //         for (size_t j = 0; j < seqs_->size(); ++j) {
    //             if (j != l)
    //                 dist[k][j] = dist[j][k] = (sk * dist[k][j] + sl * dist[l][j]) / size;
    //         }
    //         for (size_t i = 0; i < seqs_->size(); ++i)
    //             dist[i][l] = dist[l][i] = -1;

    //         LOG(ERROR) << strprintf("merging nodes %zu and %zu into node %zu (distance = %.4f  height = %.4f  length[%zu] = %.4f  length[%zu] = %.4f)",
    //                                 k, l, p, dmin, height, k, bk, l, bl);
    //         if (verbose_) printf("  merging nodes %zu and %zu (distance = %.4f)\n", k, l, dmin);
    //     }

    //     return tree;
    // }

    Tree BuildGuideTree() {
        KmerDistanceMatrixBuilder<Abc> dist_builder(seqs_, params_->kmer_length, verbose_);
        Matrix<float> dist = dist_builder.Build();
        Tree tree;
        Vector<Node> nodes(seqs_->size());  // nodes[i] is node of i-th cluster
        Bitset disabled(seqs_->size());     // indicates disabled indices

        if (verbose_) puts("Building UPGMA guide tree ...");

        // Add leaf nodes to guide tree
        for (size_t k = 0; k < seqs_->size(); ++k)
            nodes[k] = add_vertex(TreeNodeProperties(), tree);

        double branch_len = 0.0857;
        Node p, l = 0;
        for (size_t r = 1; r < seqs_->size(); ++r) {
            p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
            add_edge(p, l, TreeEdgeProperties(branch_len), tree);
            add_edge(p, r, TreeEdgeProperties(branch_len * r), tree);
            l = p;
        }

        // // Add node 8
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 0, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 1, TreeEdgeProperties(branch_len), tree);

        // // Add node 9
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 2, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 8, TreeEdgeProperties(branch_len), tree);

        // // Add node 10
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 4, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 5, TreeEdgeProperties(branch_len), tree);

        // // Add node 11
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 6, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 7, TreeEdgeProperties(branch_len), tree);

        // // Add node 12
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 8, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 9, TreeEdgeProperties(branch_len), tree);

        // // Add node 13
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 10, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 11, TreeEdgeProperties(branch_len), tree);

        // // Add node 14
        // p = add_vertex(TreeNodeProperties(0, 0, 0, 0), tree);
        // add_edge(p, 12, TreeEdgeProperties(branch_len), tree);
        // add_edge(p, 13, TreeEdgeProperties(branch_len), tree);

        return tree;
    }

    // Resets state such that we can start a new progressive alignment
    void Reset() {
        // Free any previously constructed PO-HMMs
        for (typename HmmVec::iterator it = hmms_.begin(); it != hmms_.end(); ++it)
            delete *it;
        hmms_.clear();
    }

  private:
    // Parameter wrapper
    const ProgressiveAlignerParams* params_;
    // Sequences to be aligned
    const SeqVec* seqs_;
    // Substitution matrix for PCs and scoring
    const SubstitutionMatrix<Abc>* mat_;
    // Pseudocount factory for profile probs in PO-HMMs
    const Pseudocounts<Abc>* pc_;
    // Context library for context scoring
    ContextLibrary<Abc>* lib_;
    // Emission functor needed for context assignments
    scoped_ptr<Emission<Abc> > emission_;
    // PO-HMM aligner that does all the work
    scoped_ptr<PairwiseAligner<Abc> > aligner_;
    // Pointers to PO-HMMs as they are built
    HmmVec hmms_;
    // Be talkative
    bool verbose_;
    // Perform expensive sanity check for each PO-HMM
    bool sanity_;
    // dir/basename of PDFs to be generated for each internal node
    std::string pdf_basename_;
};

}  // namespace cs

#endif  // CS_PROGRESSIVE_ALIGNER_H_
