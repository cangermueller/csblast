/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alignment.h"

namespace cs
{

SequenceAlignment::SequenceAlignment(int nseqs,
                                 int ncols,
                                 const SequenceAlphabet* alphabet)
        : ColumnMajorMatrix<char>(nseqs, ncols),
          headers_(nseqs, ""),
          alphabet_(alphabet)
{}

SequenceAlignment::SequenceAlignment(std::istream& in,
                                     const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ in >> *this; }

SequenceAlignment::~SequenceAlignment()
{}

void SequenceAlignment::init(std::istream& in)
{
    std::vector< SmartPtr<Sequence> > sequences(cs::Sequence::read(in, alphabet_));

    if (sequences.empty()) throw MyException("Unable to initialize alignment: no aligned sequences found!");
    const int seqs = sequences.size();
    const int cols = (*sequences[0]).length();
    for (int i = 1; i < seqs; ++i) {
        const int k = (*sequences[i]).length();
        if (k != cols) throw MyException("Bad alignment format: sequence %i has length %i but should have length %i!", i, k, cols);
    }

    headers_.resize(seqs);
    resize(seqs, cols);
    for (int i = 0; i < seqs; ++i) {
        headers_[i] = (*sequences[i]).header();
        for (int j = 0; j < cols; ++j)
            (*this)(i,j) = (*sequences[i])(j);
    }
}

std::istream& operator>> (std::istream& in, SequenceAlignment& alignment)
{
    alignment.init(in);
    return in;
}

std::ostream& operator<< (std::ostream& out, const SequenceAlignment& alignment)
{
    const int kLineLength = 80;
    const int nseqs = alignment.nseqs();
    const int ncols = alignment.ncols();

    for (int i = 0; i < nseqs; ++i) {
        out << '>' << alignment.header(i) << std::endl;
        for (int j = 0; j < ncols; ++j) {
            out << alignment.chr(i,j);
            if ((j+1) % kLineLength == 0) out << std::endl;
        }
    }
    return out;
}

// Calculates column specific sequence weights from subalignments within the global alignment.
ColumnMajorMatrix<float> column_specific_sequence_weights(const SequenceAlignment& alignment)
{
    //TODO
}

// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
std::vector<float> global_sequence_weights(const SequenceAlignment& alignment)
{

    //TODO
}

}//cs


