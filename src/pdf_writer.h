// Copyright 2010, Andreas Biegert

#ifndef CS_PDF_WRITER_H_
#define CS_PDF_WRITER_H_

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "abstract_state_matrix-inl.h"
#include "pairwise_aligner-inl.h"
#include "po_hmm-inl.h"
#include "context_library-inl.h"
#include "crf-inl.h"

namespace cs {

// Abstract base class for various PDF-writer subclasses using Latex and TikZ
class PdfWriter {
  public:
    PdfWriter() {}
    virtual ~PdfWriter() {}

    // Writes PDF to given outfile
    void WriteToFile(std::string outfile, bool keep = false);

    // TexShade residue width in cm
    static const float kResWidth = 0.21212;

    // Writes header section of latex file
    void WriteHeader(FILE* fp) const;

    // Writes footer section of latex file
    void WriteFooter(FILE* fp) const;

    // Writes latex document body. This method should be implemented by subclasses.
    virtual void WriteBody(FILE* fp) = 0;

    // Create figures in external directory
    std::string external_dir;
};

// PDF-Writer subclass for drawing a colored multiple sequence alignment.
template<class Abc>
class AlignmentPdfWriter : public PdfWriter {
  public:
    AlignmentPdfWriter(const Alignment<Abc>& ali);
    virtual ~AlignmentPdfWriter() { remove(tempfile_.c_str()); }

  private:
    // Writes color-coded alignment using the Texshade package.
    virtual void WriteBody(FILE* fp);

    const Alignment<Abc>& ali_;
    std::string tempfile_;
};


// PDF-Writer subclass for visualizing a DNA or AA profile with profile logos.
template<class Abc>
class ProfilePdfWriter : public PdfWriter {
  public:
    ProfilePdfWriter(const Profile<Abc>& profile);
    virtual ~ProfilePdfWriter() {}

    double col_height;
    bool hide_numbering;
    bool show_chars;
    double pmin;

  private:
    typedef boost::tuple<int, double, uint8_t> GroupProbResTuple;

    // Writes TikZ picture with profile logos
    virtual void WriteBody(FILE* fp);

    // Y-Scale for TikZ
    static const float kScaleY = 0.1;

    Profile<Abc> profile_;
};


// PDF-Writer subclass for visualizing an abstract state profile with profile logos.
template<class AS, class Abc>
class StateProfilePdfWriter : public PdfWriter {
  public:
    StateProfilePdfWriter(const Profile<AS>& profile,
			  const ContextLibrary<Abc>& lib);
    virtual ~StateProfilePdfWriter() {}

    bool hide_numbering;
    double pmin;

  private:
    typedef boost::tuple<Color, double, uint8_t> ColorProbStateTuple;

    // Writes TikZ picture with profile logos
    virtual void WriteBody(FILE* fp);

    // Defines colors needed for profile logos.
    void DefineColors(FILE* fp) const;

    // Y-Scale for TikZ
    static const float kScaleY = 0.1;

    Profile<AS> profile_;
    const ContextLibrary<Abc>& lib_;
};


// PDF-Writer subclass for visualizing a context library with profile logos.
template<class Abc>
class ContextLibraryPdfWriter : public PdfWriter {
  public:
    ContextLibraryPdfWriter(const ContextLibrary<Abc>& lib)
	    : pmin(0.001), margin(1.0), ncols(10), lib_(lib) {}
    virtual ~ContextLibraryPdfWriter() {}

    double pmin;
    double margin;
    size_t ncols;

  private:
    typedef boost::tuple<int, double, uint8_t> GroupProbResTuple;

    // Writes TikZ picture with profile logos
    virtual void WriteBody(FILE* fp);

    // Defines colors needed for profile logos.
    void DefineColors(FILE* fp) const;

    // Y-Scale for TikZ
    static const float kScaleY = 0.1;

    const ContextLibrary<Abc>& lib_;
};


// Helper struct for arranging PO-HMM vertices as a linear chain.
template<class Abc>
struct POHmmGrid {
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::VertexVec VertexVec;

    POHmmGrid(const POHmm<Abc>& hmm) : coords(hmm.size() + 2, 0) {
	// Determine blocks of vertices with same ali-bitset patterns
	const Graph& g = hmm.g;
	VertexVec block(1, kStartEndVertex);  // first block consists of start vertex
	for (size_t i = 1; i < hmm.vertices.size(); ++i) {
	    if (i == 1 || g[hmm.vertices[i]].alis != g[hmm.vertices[i-1]].alis) {
		blocks.push_back(block);
		block.clear();
	    }
	    block.push_back(hmm.vertices[i]);
	}
	blocks.push_back(block);
	// Add last block consisting of virtual end vertex
	blocks.push_back(VertexVec(1, hmm.vertices.size()));

	// Set vertex coordinates
	for (size_t b = 0, i = 0; b < blocks.size(); ++b)
	    for (size_t j = 0; j < blocks[b].size(); ++j, ++i)
		coords[blocks[b][j]] = b + i;
    }

    std::vector<VertexVec> blocks;  // vertex blocks
    std::vector<int> coords;        // grid coordinate for each vertex
};


// PDF-Writer subclass for drawing a complete PO-HMM with vertices, edges, and
// alternative alignments.
template<class Abc>
class POHmmPdfWriter : public PdfWriter {
  public:
    POHmmPdfWriter(const POHmm<Abc>& hmm,
		   const ContextLibrary<Abc>* lib = NULL,
		   const IndexPath* path = NULL);
    virtual ~POHmmPdfWriter();

    bool hide_numbering;
    bool hide_profiles;
    bool hide_alignment;
    bool hide_paths;
    bool reverse;
    double pmin;

  private:
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::Vertex Vertex;
    typedef typename POHmm<Abc>::VertexVec VertexVec;
    typedef typename POHmm<Abc>::ConstVertexVecIter ConstVertexVecIter;
    typedef typename POHmm<Abc>::EdgeIter EdgeIter;
    typedef SparseProfileCol::ElementIter ContextIter;
    typedef std::pair<Color, double> ColorProbPair;
    typedef boost::tuple<int, double, uint8_t> GroupProbResTuple;

    // Writes TikZ picture with PO-HMM vertices and edges.
    virtual void WriteBody(FILE* fp);

    // Generates alignment PDFs for each vertex block
    void GenerateAlignmentBlocks() const;

    // Defines colors needed for contexts state profiles.
    void DefineContextColors(FILE* fp) const;

    // Returns alignment filename for block 'b'
    std::string alignment_filename(size_t b) const {
	return strprintf("%s_%zu.pdf", tempfile_.c_str(), b);
    }

    // Colors for alignment-path boxes
    static const char* kAliColors[];
    // Y-Scale for TikZ
    static const float kScaleY = 0.07029;

    const POHmm<Abc>& hmm_;
    const ContextLibrary<Abc>* lib_;
    const IndexPath* path_;
    POHmmGrid<Abc> grid_;
    std::string tempfile_;
};


// PDF-Writer subclass for drawing posterior probability matrix, optionally with
// highlighted alignment path.
template<class Abc>
class PosteriorMatrixPdfWriter : public PdfWriter {
  public:
    PosteriorMatrixPdfWriter(const AlignmentMatrices<Abc>& data,
			     const ContextLibrary<Abc>* lib = NULL,
			     const PairAlignment* ali = NULL);
    virtual ~PosteriorMatrixPdfWriter();

  private:
    typedef typename POHmm<Abc>::Graph Graph;
    typedef typename POHmm<Abc>::Vertex Vertex;
    typedef typename POHmm<Abc>::VertexVec VertexVec;
    typedef typename POHmm<Abc>::ConstVertexVecIter ConstVertexVecIter;

    // Writes posterior matrix using the Texshade package.
    virtual void WriteBody(FILE* fp);

    const AlignmentMatrices<Abc>& data_;
    const ContextLibrary<Abc>* lib_;
    const PairAlignment* ali_;
    POHmmGrid<Abc> grid_x_;
    POHmmGrid<Abc> grid_y_;
    std::string tempfile_;
};


// PDF-Writer subclass for drawing an abstract state substitution matrix.
template<class AS, class Abc>
class AbstractStateMatrixPdfWriter : public PdfWriter {
  public:
    AbstractStateMatrixPdfWriter(const AbstractStateMatrix<AS>& matrix,
				 const ContextLibrary<Abc>& contexts,
				 const ContextLibrary<Abc>& alphabet);
    virtual ~AbstractStateMatrixPdfWriter();

    double score_min;
    double score_max;

  private:
    // Writes substitution scores of matrix using the Texshade package.
    virtual void WriteBody(FILE* fp);

    // Returns RGB color for substitution score of context k and abstract state a
    Color GetColor(size_t k, size_t a) const;

    // Defines colors needed for context name badges.
    void DefineContextColors(FILE* fp) const;

    const AbstractStateMatrix<AS>& matrix_;
    const ContextLibrary<Abc>& contexts_;
    const ContextLibrary<Abc>& alphabet_;
    std::string tempfile_;
};

inline std::string PairStateToColor(uint8_t xx) {
    std::string rv;
    switch (xx) {
	case MM: rv = "red"; break;
	case MI: rv = "MyGreen!90!white"; break;
	case IM: rv = "MyBlue!90!white"; break;
	default: throw Exception("Unknown pair state '%d'!", xx);
    }
    return rv;
}

// PDF-Writer subclass for drawing an CRF states.
template<class Abc>
class CrfStatePdfWriter : public PdfWriter {
  public:
    CrfStatePdfWriter() : prob_min(0.005) {}

    const CrfState<Abc>* crf_state;
    double prob_min;

    virtual size_t Width() const = 0;
    virtual size_t Height() const = 0;

  protected:
    typedef boost::tuple<int, double, uint8_t> GroupProbResTuple;
  
    std::vector<size_t> SortRes(const double* col) const;
};

template<class Abc>
class ProbCrfStatePdfWriter : public CrfStatePdfWriter<Abc> {
  public:
    ProbCrfStatePdfWriter() 
      : col_height(10.0),
        show_weights(true),
        weights_height(3.0),
        weight_center(0.0),
        weight_decay(0.85),
        show_name(true) {}

    virtual ~ProbCrfStatePdfWriter() {};

    virtual void WriteBody(FILE* fp);
    virtual size_t Width() const;
    virtual size_t Height() const;

    double col_height;
    double show_weights;
    double weights_height;
    double weight_center;
    double weight_decay;
    bool show_name;

  private:
    void DrawColumn(FILE* fp, double* probs, double x, double y) const;
};

template<class Abc>
class WeightCrfStatePdfWriter : public CrfStatePdfWriter<Abc> {
  public:
    WeightCrfStatePdfWriter() 
      : col_height(10.0),
        show_name(true) {}

    virtual ~WeightCrfStatePdfWriter() {};

    virtual void WriteBody(FILE* fp);
    virtual size_t Width() const;
    virtual size_t Height() const;

    double col_height;
    bool show_name;

  private:
    void DrawColumn(FILE* fp, double* probs, double x, double y, bool down) const;
};

template<class Abc>
class CrfPdfWriter : public PdfWriter {
  public:
    CrfPdfWriter(
        const std::vector<CrfState<Abc> >& crf_states, 
        CrfStatePdfWriter<Abc>& crf_state_writer)
      : states_per_row(10),
        crf_states_(crf_states),
        crf_state_writer_(crf_state_writer) {}

    virtual ~CrfPdfWriter() {};

    void WriteBody(FILE* fp);

    size_t states_per_row;
  
  private:
    const std::vector<CrfState<Abc> >& crf_states_;
    CrfStatePdfWriter<Abc>& crf_state_writer_;
};

}  // namespace cs

#endif  // CS_PDF_WRITER_H_
