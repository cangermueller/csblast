/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "alignment.h"

#include <cctype>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <valarray>
#include <vector>
#include <string>

#include "matrix.h"
#include "exception.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"
#include "util.h"

namespace
{

// Debug flag
const bool kDebug = false;

// Converts a character to uppercase and '.' to '-'.
inline char match_chr(char c) { return (isalpha(c) ? toupper(c) : (c == '.' ? '-' : c)); }

// Converts a character to lowercase and '-' to '.'.
inline char insert_chr(char c) { return (isalpha(c) ? tolower(c) : (c == '-' ? '.' : c)); }

// Predicate indicating if character belongs to match column.
inline bool is_match_chr(char c) { return isalpha(c) && isupper(c) || c == '-'; }

// Predicate indicating if character belongs to insert column.
inline char is_insert_chr(char c) { return isalpha(c) && islower(c) || c == '.'; }

} // namespace

namespace cs
{

Alignment::Alignment(std::istream& in, InputFormat format, const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ read(in, format); }

void Alignment::init(std::vector<std::string> headers, std::vector<std::string> seqs, bool case_insens)
{
    if (seqs.empty()) throw Exception("Unable to initialize alignment: no aligned sequences found!");
    if (headers.size() != seqs.size())
        throw Exception("Unable to initialize alignment: unequal number of headers and sequences!");
    // Remove whitespace from sequences
    for (std::vector<std::string>::iterator iter = seqs.begin(); iter != seqs.end(); ++iter)
        iter->erase(remove_if(iter->begin(), iter->end(), isspace), iter->end());

    const int nseqs = seqs.size();
    const int ncols = seqs[0].length();
    for (int k = 1; k < nseqs; ++k)
        if (static_cast<int>(seqs[k].length()) != ncols)
            throw Exception("Bad alignment: sequence %i has length %i but should have %i!", k, seqs[k].length(), ncols);

    // Validate characters and convert to integer representation
    headers_.resize(nseqs);
    resize(seqs.size(), seqs[0].length());
    for (int k = 0; k < nseqs; ++k) {
        headers_[k] = headers[k];
        for (int i = 0; i < ncols; ++i) {
            const char c = seqs[k][i];
            column_indexes_[i] = i;
            mask_[i] = case_insens || is_match_chr(c);
            if (alphabet_->valid(c, true))
                seqs_[i][k] = alphabet_->ctoi(c);
            else
                throw Exception("Invalid character %c at position %i of sequence '%s'", c, i, headers_[k].c_str());
        }
    }

    // Replace gap with endgap for all gaps at either end of a sequence
    for (int k = 0; k < nseqs; ++k) {
        for (int i = 0; i < ncols && seqs_[i][k] == alphabet_->gap(); ++i)
            seqs_[i][k] = alphabet_->endgap();
        for (int i = ncols - 1; i >= 0 && seqs_[i][k] == alphabet_->gap(); --i)
            seqs_[i][k] = alphabet_->endgap();
    }

    // Initialize index array for match columns
    set_match_indexes();
}

void Alignment::set_match_indexes()
{
    const int match_cols = std::count(&mask_[0], &mask_[0] + ncols(), true);
    match_indexes.resize(match_cols);
    match_indexes = column_indexes_[mask_];
}

void Alignment::read(std::istream& in, InputFormat format)
{
    switch (format) {
        case FASTA_INPUT:
            read_fasta_or_a2m(in, true);
            break;
        case A2M_INPUT:
            read_fasta_or_a2m(in, false);
            break;
        default:
            throw Exception("Unsupported alignment input format!");
    }
}

void Alignment::read_fasta_or_a2m(std::istream& in, bool fasta)
{
    const std::string format = fasta ? "FASTA" : "A2M";
    const int kBufferSize = 1048576; //1MB
    std::string buffer;
    buffer.reserve(kBufferSize);

    std::vector<std::string> headers; // temporary for headers
    std::vector<std::string> seqs;    // temporary for sequences
    while (in.good()) {
        //read header
        if (getline(in, buffer)) {
            if (buffer.empty() ||  buffer[0] != '>')
                throw Exception("Bad format: first sequence of %s alignment does not start with '>' character!", format.c_str());
            headers.push_back(std::string(buffer.begin()+1, buffer.end()));
        } else {
            throw Exception("Failed to read from %s formatted input stream!", format.c_str());
        }
        //read sequence
        seqs.push_back("");
        while (in.peek() != '>' && getline(in, buffer))
            seqs[seqs.size()-1].append(buffer.begin(), buffer.end());
    }

    init(headers, seqs, fasta);
}

void Alignment::write(std::ostream& out, OutputFormat format, int width) const
{
    switch (format) {
        case FASTA_OUTPUT:
            write_fasta_or_a2m(out, true, width);
            break;
        case A2M_OUTPUT:
            write_fasta_or_a2m(out, false, width);
            break;
        default:
            throw Exception("Unsupported alignment output format!");
    }
}

void Alignment::write_fasta_or_a2m(std::ostream& out, bool fasta, int width) const
{
    for (int k = 0; k < nseqs(); ++k) {
        out << '>' << headers_[k] << std::endl;
        for (int i = 0; i < ncols(); ++i) {
            out << (fasta || mask_[i] ? match_chr(chr(k,i)) : insert_chr(chr(k,i)));
            if ((i+1) % width == 0) out << std::endl;
        }
        if (ncols() % width != 0) out << std::endl;
    }
}

void Alignment::resize(int nseqs, int ncols)
{
    if (nseqs == 0 || ncols == 0)
        throw Exception("Bad dimensions for alignment resizing: nseqs=%i ncols=%i", nseqs, ncols);

    seqs_.resize(ncols, nseqs);
    column_indexes_.resize(ncols);
    match_indexes.resize(ncols);
    mask_.resize(ncols);
    headers_.resize(nseqs);
}

void Alignment::assign_match_columns_by_first_sequence()
{
    for (int i = 0; i < ncols(); ++i)
        mask_[i] = seqs_[i][0] < alphabet_->gap();
    set_match_indexes();
}

void Alignment::assign_match_columns_by_gap_rule(int gap_threshold)
{
    if (kDebug) std::cerr << "Removing collumns with more than " << gap_threshold << "% of gaps:" << std::endl;

    // global weights are sufficient for calculation of gap percentage
    std::pair<std::vector<float>, float> wg_neff = global_weights_and_diversity(*this);
    for (int i = 0; i < ncols(); ++i) {
        float gap = 0.0f;
        float res = 0.0f;

        for (int k = 0; k < nseqs(); ++k)
            if (seqs_[i][k] < alphabet_->any())
                res += wg_neff.first[k];
            else
                gap += wg_neff.first[k];

        float percent_gaps = 100.0f * gap / (res + gap);

        if (kDebug) fprintf(stderr,"percent gaps[%i]=%-4.1f\n", i, percent_gaps);
        mask_[i] = percent_gaps <= static_cast<float>(gap_threshold);
    }
    set_match_indexes();
}

std::pair<std::vector<float>, float> global_weights_and_diversity(const Alignment& alignment)
{
    const float kZero = 1E-10;  // for calculation of entropy
    const int nseqs = alignment.nseqs();
    const int ncols = alignment.nmatch();
    const int nalph = alignment.alphabet()->size();
    const int any   = alignment.alphabet()->any();

    float neff = 0.0f;                        // diversity of alignment
    std::vector<float> weights(nseqs, 0.0f);  // weights of sequences
    std::vector<int> n(nseqs, 0);             // number of residues in sequence i (excl. ANY)
    std::vector<float> fj(nalph, 0.0f);       // to calculate entropy
    std::vector<int> adiff(ncols, 0);         // number of different alphabet letters in each column
    matrix<int> counts(ncols, nalph, 0);      // counts of alphabet letters in each column (excl. ANY)

    if (kDebug) fprintf(stderr,"\nCalculation of global weights and alignment diversity:\n");

    // Count number of residues in each column
    for (int i = 0; i < ncols; ++i)
        for (int k = 0; k < nseqs; ++k)
            if (alignment[i][k] < any) {
                ++counts[i][alignment[i][k]];
                ++n[k];
            }

    // Count number of different residues in each column
    for (int i = 0; i < ncols; ++i)
        for (int a = 0; a < nalph; ++a)
            if (counts[i][a]) ++adiff[i];

    // Calculate weights
    for(int i = 0; i < ncols; ++i)
        for (int k = 0; k < nseqs; ++k)
            if( adiff[i] > 0 && alignment[i][k] < any )
                weights[k] += 1.0f/( adiff[i] * counts[i][alignment[i][k]] * n[k] );
    normalize_to_one(&weights[0], nseqs);

    // Calculate number of effective sequences
    for (int i = 0; i < ncols; ++i) {
        reset(&fj[0], nalph);
        for (int k=0; k < nseqs; ++k)
            if (alignment[i][k] < any) fj[alignment[i][k]] += weights[k];
        normalize_to_one(&fj[0], nalph);
        for (int a = 0; a < nalph; ++a)
            if (fj[a] > kZero) neff -= fj[a] * log2(fj[a]);
    }
    neff = pow(2.0, neff / ncols);

    if (kDebug) fprintf(stderr,"neff=%-5.2f\n", neff);

    return make_pair(weights, neff);
}

std::pair< matrix<float>, std::vector<float> > position_specific_weights_and_diversity(const Alignment& alignment)
{
    const float kMaxEndgapFraction = 0.1;  // maximal fraction of sequences with an endgap
    const int kMinNcols = 10;  // minimum number of columns in subalignments
    const float kZero = 1E-10;  // for calculation of entropy
    const int nseqs  = alignment.nseqs();
    const int ncols  = alignment.nmatch();
    const int nalph  = alignment.alphabet()->size();
    const int any    = alignment.alphabet()->any();
    const int endgap = alignment.alphabet()->endgap();

    int ncoli = 0;        // number of columns j that contribute to neff[i]
    int nseqi = 0;        // number of sequences in subalignment i
    int ndiff = 0;        // number of different alphabet letters
    bool change = false;  // has the set of sequences in subalignment changed?
    matrix<int> n(ncols, endgap + 1, 0);  // n(j,a) = number of seq's with some residue at column i AND a at position j
    matrix<float> w(ncols, nseqs, 0.0f);  // w(i,k) weight of sequence k in column i, calculated from subalignment i
    std::vector<float> fj(nalph, 0.0f);   // to calculate entropy
    std::vector<float> neff(ncols, 0.0f); // diversity of subalignment i
    std::vector<float> wi(nseqs, 0.0f);   // weight of sequence k in column i, calculated from subalignment i
    std::pair<std::vector<float>, float> wg_neff = global_weights_and_diversity(alignment);  // global weights

    std::vector<int> nseqi_debug(ncols, 0); // debugging
    std::vector<int> ncoli_debug(ncols, 0); // debugging

    if (kDebug) {
        fprintf(stderr,"\nCalculation of position-specific weights and alignment diversity on subalignments:\n");
        fprintf(stderr,"%-5s  %-5s  %-5s  %-5s\n", "i", "ncoli", "nseqi", "neff");
    }

    for (int i = 0; i < ncols; ++i) {
        change = false;
        for (int k = 0; k < nseqs; ++k) {
            if ((i==0 || alignment[i-1][k] >= any) && alignment[i][k] < any) {
                change = true;
                ++nseqi;
                for (int j = 0; j < ncols; ++j)
                    ++n[j][alignment[j][k]];
            } else if (i>0 && alignment[i-1][k] < any && alignment[i][k] >= any) {
                change = true;
                --nseqi;
                for (int j = 0; j < ncols; ++j)
                    --n[j][alignment[j][k]];
            }
        }  // for k over nseqs
        if (kDebug) nseqi_debug[i] = nseqi;

        if (change) {  // set of sequences in subalignment has changed
            ncoli = 0;
            reset(&wi[0], nseqs);

            for (int j = 0; j < ncols; ++j) {
                if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
                ndiff = 0;
                for (int a = 0; a < nalph; ++a) if (n[j][a]) ++ndiff;
                if (ndiff == 0) continue;
                ++ncoli;
                for (int k = 0; k < nseqs; ++k) {
                    if (alignment[i][k] < any && alignment[j][k] < any) {
                        if (kDebug && n[j][alignment[j][k]] == 0) {
                            fprintf(stderr, "Error: Mi=%i: n[%i][seqs[%i][%i]]=0! (seqs[%i][%i]=%c)\n",
                                    i, j, k, j, k, j, alignment[j][k] );
                        }
                        wi[k] += 1.0f / static_cast<float>((n[j][alignment[j][k]] * ndiff));
                    }
                }
            }  // for j over ncols
            normalize_to_one(&wi[0], nseqs);

            if (ncoli < kMinNcols)  // number of columns in subalignment insufficient?
                for (int k = 0; k < nseqs; ++k)
                    if (alignment[i][k] < any)
                        wi[k] = wg_neff.first[k];
                    else
                        wi[k] = 0.0f;

            neff[i] = 0.0f;
            for (int j = 0; j < ncols; ++j) {
                if (n[j][endgap] > kMaxEndgapFraction * nseqi) continue;
                reset(&fj[0], nalph);

                for (int k = 0; k < nseqs; ++k)
                    if (alignment[i][k] < any && alignment[j][k] < any)
                        fj[alignment[j][k]] += wi[k];
                normalize_to_one(&fj[0], nalph);
                for (int a = 0; a < nalph; ++a)
                    if (fj[a] > kZero) neff[i] -= fj[a] * log2(fj[a]);
            }  // for j over ncols

            neff[i] = (ncoli > 0 ? pow(2.0, neff[i] / ncoli) : 1.0f);

        } else {  // set of sequences in subalignment has NOT changed
            neff[i] = (i==0 ? 0.0 : neff[i-1]);
        }

        for (int k = 0; k < nseqs; ++k) w[i][k] = wi[k];
        if (kDebug) {
            ncoli_debug[i] = ncoli;
            fprintf(stderr,"%-5i  %-5i  %-5i  %-5.2f\n", i, ncoli_debug[i], nseqi_debug[i], neff[i]);
        }
    }  // for i over ncols

    return make_pair(w, neff);
}

}//cs


