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
    const int kBufferSize = 1048576; //1MB
    std::string buffer;
    buffer.reserve(kBufferSize);
    std::vector< std::string >().swap(headers_);
    std::vector< std::vector<char> > sequences;

    while (in.good()) {
        //read header
        if (getline(in, buffer)) {
            if (buffer.empty() ||  buffer[0] != '>')
                throw MyException("Bad format: first line of aligned FASTA sequence does not start with '>' character!");
            headers_.push_back(std::string(buffer.begin()+1, buffer.end()));
        } else {
            throw MyException("Failed to read from FASTA formatted input stream!");
        }
        //read sequence
        std::vector<char> sequence;
        while (in.peek() != '>' && getline(in, buffer))
            sequence.insert(sequence.end(), buffer.begin(), buffer.end());
        //strip whitespace and newlines from sequence.
        sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace),
                       sequence.end());
        sequences.push_back(sequence);
    }

    if (sequences.empty()) throw MyException("Unable to initialize alignment: no aligned sequences found!");
    if (headers_.size() != sequences.size()) throw MyException("Unequal number of headers and sequences!");
    const int seqs = sequences.size();
    const int cols = sequences[0].size();
    for (int i = 1; i < seqs; ++i)
        if (static_cast<int>(sequences[i].size()) != cols)
            throw MyException("Bad alignment format: sequence %i has length %i but should have length %i!", i, sequences[i].size(), cols);

    //validate characters and convert to integer representation
    resize(seqs, cols);
    for (int i = 0; i < seqs; ++i)
        for (int j = 0; j < cols; ++j) {
            const char c = sequences[i][j];
            if (c == kGap)
                (*this)(i,j) = gaptoi();
            else if (alphabet_->valid(c))
                (*this)(i,j) = alphabet_->ctoi(c);
            else
                throw MyException("Invalid character %c at position %i of sequence '%s'", c, j, headers_[i].c_str());
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


