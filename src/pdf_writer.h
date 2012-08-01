/*
  Copyright 2010-2012 Andreas Biegert, Christof Angermueller

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

#ifndef CS_PDF_WRITER_H_
#define CS_PDF_WRITER_H_

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "abstract_state_matrix-inl.h"
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
