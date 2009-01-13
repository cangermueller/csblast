/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence_alignment.h"

#include <cctype>
#include <cstdio>

#include <iostream>
#include <vector>

#include "matrix.h"
#include "my_exception.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "smart_ptr.h"
#include "util.h"

namespace
{

const bool kDebug = false;

} // namespace

namespace cs
{

SequenceAlignment::SequenceAlignment(std::istream& in, const SequenceAlphabet* alphabet)
        : nseqs_(0),
          ncols_(0),
          alphabet_(alphabet)
{ in >> *this; }

void SequenceAlignment::init(std::istream& in)
{
    const int kBufferSize = 1048576; //1MB
    std::string buffer;
    buffer.reserve(kBufferSize);
    std::vector< std::string >().swap(headers_);
    std::vector< std::vector<char> > sequences;
    std::vector<char> sequence;

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
        while (in.peek() != '>' && getline(in, buffer))
            sequence.insert(sequence.end(), buffer.begin(), buffer.end());
        //strip whitespace and newlines from sequence.
        sequence.erase(remove_if(sequence.begin(), sequence.end(), isspace),
                       sequence.end());
        sequences.push_back(sequence);
        std::vector<char>().swap(sequence); //clear for next use
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
                (*this)(i,j) = gap();
            else if (alphabet_->valid(c))
                (*this)(i,j) = alphabet_->ctoi(c);
            else
                throw MyException("Invalid character %c at position %i of sequence '%s'", c, j, headers_[i].c_str());
        }

    set_endgaps();  // Replace gap with endgap for all gaps at either end of a sequence
}

void SequenceAlignment::set_endgaps()
{
    for (int i = 0; i < nseqs_; ++i) {
        for (int j = 0; j < ncols_ && gap(i,j); ++j)   (*this)(i,j) = endgap();
        for (int j = ncols_-1; j >=0 && gap(i,j); --j) (*this)(i,j) = endgap();
    }
}

void SequenceAlignment::resize(int nseqs, int ncols)
{
    if (nseqs == 0 || ncols == 0)
        throw MyException("Bad dimensions for alignment resizing: nseqs=%i ncols=%i", nseqs, ncols);
    nseqs_ = nseqs;
    ncols_ = ncols;
    sequences_.resize(nseqs * ncols);
    headers_.resize(nseqs);
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

std::vector<float> global_weights(const SequenceAlignment& alignment)
{
    const int nseqs = alignment.nseqs();
    const int ncols = alignment.ncols();
    const int nalph = alignment.alphabet()->size()-1;  // alphabet size without ANY character
    const int any   = alignment.alphabet()->any();

    std::vector<float> weights(nseqs, 0.0f);  // weights of sequences
    std::vector<int> adiff(ncols, 0);         // number of different alphabet letters in each column
    std::vector<int> residues(nseqs, 0);      // number of residues in sequence i (excl. ANY)
    Matrix<int> counts(ncols, nalph, 0);      // counts of alphabet letters in each column (excl. ANY)

    // Count number of residues in each column
    for (int i = 0; i < nseqs; ++i) {
        for (int j = 0; j < ncols; ++j) {
            if (alignment(i,j) < any) {
                ++counts(j, alignment(i,j));
                ++residues[i];
            }
        }
    }

    // Count number of different residues in each column
    for (int j = 0; j < ncols; ++j)
        for (int a = 0; a < nalph; ++a)
            if (counts(j,a)) ++adiff[j];

    // Calculate weights
    for (int i = 0; i < nseqs; ++i)
        for( int j = 0; j < ncols; ++j)
            if( adiff[j] > 0 && alignment(i,j) < any )
                weights[i] += 1.0f/( adiff[j] * counts(j, alignment(i,j)) * residues[i] );

    normalize_to_one(&weights[0], weights.size());

    return weights;
}

std::pair< Matrix<float>, std::vector<float> > position_dependent_weights_and_neff(const SequenceAlignment& alignment)
{
    const float kMaxEndgapFraction = 0.1;  // Maximal fraction of sequences with an endgap
    const int kMinNcols = 10;              // Minimum number of columns in subalignments
    const float kZero = 1E-10;             // Zero for calculation of entropy

    const int nseqs  = alignment.nseqs();
    const int ncols  = alignment.ncols();
    const int nalph  = alignment.alphabet()->size()-1;  // alphabet size without ANY character
    const int any    = alignment.alphabet()->any();
    const int endgap = alignment.endgap();

    int ncoli = 0;        // number of columns j that contribute to neff[i]
    int nseqi = 0;        // number of sequences in subalignment i
    int ndiff = 0;        // number of different alphabet letters
    bool change = false;  // has the set of sequences in subalignment changed?
    Matrix<int> n(ncols, endgap+1, 0);    // n(j,a) = number of seq's with some residue at column i AND a at position j
    Matrix<float> w(ncols, nseqs, 0.0f);  // w(i,k) weight of sequence k in column i, calculated from subalignment i
    std::vector<float> fj(nalph, 0.0f);   // to calculate entropy
    std::vector<float> neff(ncols, 0.0f); // diversity of subalignment i
    std::vector<float> wi(nseqs, 0.0f);   // weight of sequence k in column i, calculated from subalignment i
    std::vector<float> wg(global_weights(alignment));  // global weight of sequence k

    std::vector<int> nseqi_debug(ncols, 0); // debugging
    std::vector<int> ncoli_debug(ncols, 0); // debugging

    for (int i = 0; i < ncols; ++i) {
        change = false;
        for (int k = 0; k < nseqs; ++k) {
            if ((i==0 || alignment(k,i-1) >= any) && alignment(k,i) < any) {
                change = true;
                ++nseqi;
                for (int j = 0; j < ncols; ++j)
                    ++n(j, alignment(k,j));
            } else if (i>0 && alignment(k,i-1) < any && alignment(k,i) >= any) {
                change = true;
                --nseqi;
                for (int j=0; j < ncols; ++j)
                    --n(j, alignment(k,j));
            }
        }  // for k over nseqs
        nseqi_debug[i] = nseqi;

        if (change) {  // set of sequences in subalignment has changed
            ncoli = 0;
            reset(&wi[0], nseqs);

            for (int j = 0; j < ncols; ++j) {
                if (n(j,endgap) > kMaxEndgapFraction * nseqi) continue;
                ndiff = 0;
                for (int a = 0; a < nalph; ++a) if (n(j,a)) ++ndiff;
                if (ndiff == 0) continue;
                ++ncoli;
                for (int k = 0; k < nseqs; ++k) {
                    if (alignment(k,i) < any && alignment(k,j) < any) {
                        if (kDebug && n(j, alignment(k,j)) == 0) {
                            fprintf(stderr, "Error: Mi=%i: n[%i][seqs[%i][%i]]=0! (seqs[%i][%i]=%c)\n",
                                    i, j, k, j, k, j, alignment(k,j) );
                        }
                        wi[k] += 1.0f / (n(j, alignment(k,j)) * ndiff);
                    }
                }
            }  // for j over ncols
            normalize_to_one(&wi[0], nseqs);

            if (ncoli < kMinNcols)  // number of columns in subalignment insufficient?
                for (int k = 0; k < nseqs; ++k)
                    if (alignment(k,i) < any)
                        wi[k] = wg[k];
                    else
                        wi[k] = 0.0f;

            neff[i] = 0.0f;
            for (int j = 0; j < ncols; ++j) {
                if (n(j, endgap) > kMaxEndgapFraction * nseqi) continue;
                reset(&fj[0], nalph);

                for (int k=0; k < nseqs; ++k)
                    if (alignment(k,i) < any && alignment(k, j) < any)
                        fj[alignment(k,j)] += wi[k];
                normalize_to_one(&fj[0], nalph);

                for (int a = 0; a < nalph; ++a)
                    if (fj[a] > kZero) neff[i] -= fj[a] * log2(fj[a]);
            }  // for j over ncols

            neff[i] = (ncoli > 0 ? pow(2.0, neff[i] / ncoli) : 1.0f);

        } else {  // set of sequences in subalignment has NOT changed
            neff[i] = (i==0 ? 0.0 : neff[i-1]);
        }

        for (int k = 0; k < nseqs; ++k) w(i,k) = wi[k];
        ncoli_debug[i] = ncoli;
    }  // for i over ncols

    if (kDebug) {
        fprintf(stderr,"\nCalculation of position-dependent weights on subalignments:\n");
        fprintf(stderr,"%-5s  %-5s  %-5s  %5s\n", "i", "ncoli", "nseqi", "neffi");
        for (int i = 0; i < ncols; ++i)
            fprintf(stderr,"%-5i  %-5i  %-5i  %-5.2f\n", i, ncoli_debug[i], nseqi_debug[i], neff[i]);
    }

    return make_pair(w, neff);
}

}//cs


