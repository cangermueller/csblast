#ifndef CS_SEQUENCE_ALIGNMENT_H
#define CS_SEQUENCE_ALIGNMENT_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for an alignment of sequences over an alphabet.

#include <cctype>
#include <vector>

#include "smart_ptr.h"
#include "column_major_matrix.h"
#include "sequence.h"
#include "my_exception.h"

namespace cs
{

class SequenceAlignment : private ColumnMajorMatrix<char>
{
  public:
    friend std::istream& operator>> (std::istream& in, SequenceAlignment& alignment);

    SequenceAlignment(int nseqs,
                      int ncols,
                      const SequenceAlphabet* alphabet);
    SequenceAlignment(std::istream& in,
                      const SequenceAlphabet* alphabet);
    virtual ~SequenceAlignment();

    // Returns the integer representation of the character in column j of sequence i (starting at index 0)
    using ColumnMajorMatrix<char>::operator();
    // Returns the character in column j of sequence i (starting at index 0)
    char chr(int i, int j) const;
    // Returns true if the character at position (i,j) is a gap.
    bool gap(int i, int j) const;
    // Returns the number of sequences in the alignment.
    int nseqs() const;
    // Returns the number of alignment columns.
    int ncols() const;
    // Returns the header of sequence i as mutable reference.
    std::string& header(int i);
    // Returns the header of sequence i as const reference.
    const std::string& header(int i) const;
    // Returns the underlying sequence alphabet.
    const SequenceAlphabet& alphabet() const;

  private:
    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void init(std::istream& in);
    // Returns integer representation of gap character.
    int gaptoi() const;

    static const char kGap = '-';
    std::vector< std::string > headers_;
    const SequenceAlphabet* alphabet_;
};//SequenceAlignment



// Initializes a sequence alignment object from FASTA formatted alignment in input stream.
std::istream& operator>> (std::istream& in, SequenceAlignment& alignment);

// Prints the alignment in multi FASTA format to output stream.
std::ostream& operator<< (std::ostream& out, const SequenceAlignment& alignment);

// Calculates column specific sequence weights from subalignments within the global alignment.
ColumnMajorMatrix<float> column_specific_sequence_weights(const SequenceAlignment& alignment);

// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
std::vector<float> global_sequence_weights(const SequenceAlignment& alignment);



inline char SequenceAlignment::chr(int i, int j) const
{ return gap(i,j) ? kGap : alphabet_->itoc((*this)(i,j)); }

inline bool SequenceAlignment::gap(int i, int j) const
{ return (*this)(i,j) == gaptoi(); }

inline int SequenceAlignment::nseqs() const
{ return ColumnMajorMatrix<char>::nrows(); }

inline int SequenceAlignment::ncols() const
{ return ColumnMajorMatrix<char>::ncols(); }

inline std::string& SequenceAlignment::header(int i)
{ return headers_[i]; }

inline const std::string& SequenceAlignment::header(int i) const
{ return headers_[i]; }

inline const SequenceAlphabet& SequenceAlignment::alphabet() const
{ return *alphabet_; }

inline int SequenceAlignment::gaptoi() const
{ return alphabet_->size(); }

}//cs

#endif
