/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "alignment.h"

#include <cctype>
#include <cstdio>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <vector>
#include <string>

#include "log.h"
#include "matrix.h"
#include "exception.h"
#include "sequence.h"
#include "sequence_alphabet.h"
#include "shared_ptr.h"
#include "util.h"

using std::setw;
using std::left;

namespace
{

// Converts a character to uppercase and '.' to '-'.
inline char to_match_chr(char c)
{
    return (isalpha(c) ? toupper(c) : (c == '.' ? '-' : c));
}

// Converts a character to lowercase and '-' to '.'.
inline char to_insert_chr(char c)
{
    return (isalpha(c) ? tolower(c) : (c == '-' ? '.' : c));
}

// Predicate indicating if character belongs to match column.
inline bool match_chr(char c)
{
    return isalpha(c) && isupper(c) || c == '-';
}

// Predicate indicating if character belongs to insert column.
inline char insert_chr(char c)
{
    return isalpha(c) && islower(c) || c == '.';
}

} // namespace


namespace cs
{

Alignment::Alignment(std::istream& in, Format format, const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{
    read(in, format);
}

void Alignment::init(const std::vector<std::string>& headers, const std::vector<std::string>& seqs)
{
    if (seqs.empty()) throw Exception("Unable to initialize alignment: no aligned sequences found!");
    if (headers.size() != seqs.size())
        throw Exception("Unable to initialize alignment: unequal number of headers and sequences!");

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
            match_column_[i] = match_chr(c);
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
    const int match_cols = std::count(&match_column_[0], &match_column_[0] + ncols(), true);
    match_indexes.resize(match_cols);
    match_indexes = column_indexes_[match_column_];
}

void Alignment::read(std::istream& in, Format format)
{
    LOG(DEBUG1) << "Reading alignment from stream ...";

    std::vector<std::string> headers;
    std::vector<std::string> seqs;

    switch (format) {
        case FASTA:
            read_fasta(in, headers, seqs);
            break;
        case A2M:
            read_a2m(in, headers, seqs);
            break;
        case A3M:
            read_a3m(in, headers, seqs);
            break;
        default:
            throw Exception("Unsupported alignment input format %i!", format);
    }

    init(headers, seqs);

    LOG(DEBUG1) << *this;
}

void Alignment::read_fasta_flavors(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    headers.clear();
    seqs.clear();
    std::string buffer;

    while (in.good()) {
        // Read header
        if (getline(in, buffer)) {
            if (buffer.empty() ||  buffer[0] != '>')
                throw Exception("Bad format: header of sequence %i does not start with '>'!", headers.size() + 1);
            headers.push_back(std::string(buffer.begin()+1, buffer.end()));
        } else {
            throw Exception("Failed to read from alignment input stream!");
        }

        // Read sequence
        seqs.push_back("");
        while (in.peek() != '>' && getline(in, buffer)) {
            if (buffer.empty()) break;
            seqs.back().append(buffer.begin(), buffer.end());
        }
        // Remove whitespace from sequence
        seqs.back().erase(remove_if(seqs.back().begin(), seqs.back().end(), isspace), seqs.back().end());
        // Terminate reading when empty line encountered
        if (buffer.empty()) break;
    }
}

void Alignment::read_fasta(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);
    // convert all characters to match characters
    for (std::vector<std::string>::iterator it = seqs.begin(); it != seqs.end(); ++it)
        transform(it->begin(), it->end(), it->begin(), to_match_chr);
}

void Alignment::read_a2m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);
}

void Alignment::read_a3m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);

    // Check number of match states
    const int nseqs = seqs.size();
    const int nmatch = count_if(seqs[0].begin(), seqs[0].end(), match_chr);
    for (int k = 1; k < nseqs; ++k) {
        const int nmatch_k = count_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (nmatch_k != nmatch)
            throw Exception("Bad alignment: sequence %i has %i match columns but should have %i!", k, nmatch_k, nmatch);
        if (count(seqs[k].begin(), seqs[k].end(), '.') > 0)
            throw Exception("Bad alignment: sequence %i in A3M formatted alignment contains gap characters '.'!", k);
    }

    // Insert gaps into A3M alignment
    std::vector<std::string> seqs_a2m(seqs.size(), "");     // sequences in A2M format
    matrix<std::string> inserts(seqs.size(), nmatch, "");   // inserts of sequence k after match state i
    std::vector<int> max_insert_len(nmatch, 0);             // maximal length of inserts after match state i

    // Move inserts BEFORE first match state into seqs_a2m
    for (int k = 0; k < nseqs; ++k) {
        std::string::iterator i = find_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (i != seqs[k].end()) {
            seqs_a2m[k].append(seqs[k].begin(), i);
            seqs[k].erase(seqs[k].begin(), i);
        }
    }

    // Extract all inserts and keep track of longest insert after each match column
    for (int k = 0; k < nseqs; ++k) {
        int i = -1;
        std::string::iterator insert_end = seqs[k].begin();
        std::string::iterator insert_start = find_if(insert_end, seqs[k].end(), insert_chr);
        while (insert_start != seqs[k].end() && insert_end != seqs[k].end()) {
            i += insert_start - insert_end;
            insert_end = find_if(insert_start, seqs[k].end(), match_chr);
            inserts[k][i] = std::string(insert_start, insert_end);
            max_insert_len[i] = std::max(static_cast<int>(inserts[k][i].length()), max_insert_len[i]);
            insert_start = find_if(insert_end, seqs[k].end(), insert_chr);
        }
    }

    // Build new A2M alignment
    for (int k = 0; k < nseqs; ++k) {
        int i = 0;
        std::string::iterator match = seqs[k].begin();
        while (match != seqs[k].end()) {
            seqs_a2m[k].append(1, *match);
            if (max_insert_len[i] > 0) {
                seqs_a2m[k].append(inserts[k][i]);
                seqs_a2m[k].append(max_insert_len[i] - inserts[k][i].length(), '.');
            }
            match = find_if(match+1, seqs[k].end(), match_chr);
            ++i;
        }
    }

    // Overwrite original A3M alignment with new A2M alignment
    seqs = seqs_a2m;
}

void Alignment::write(std::ostream& out, Format format, int width) const
{
    switch (format) {
        case FASTA:
        case A2M:
        case A3M:
            write_fasta_flavors(out, format, width);
            break;
        case CLUSTAL:
        case PSI:
            write_clustal_flavors(out, format, width);
            break;
        default:
            throw Exception("Unsupported alignment output format %i!", format);
    }
}

void Alignment::write_fasta_flavors(std::ostream& out, Format format, int width) const
{
    for (int k = 0; k < nseqs(); ++k) {
        out << '>' << headers_[k] << std::endl;
        int j = 0;  // counts printed characters
        for (int i = 0; i < ncols(); ++i) {
            switch (format) {
                case FASTA:
                    out << to_match_chr(chr(k,i));
                    ++j;
                    break;
                case A2M:
                    out << (match_column_[i] ? to_match_chr(chr(k,i)) : to_insert_chr(chr(k,i)));
                    ++j;
                    break;
                case A3M:
                    if (match_column_[i]) {
                        out << to_match_chr(chr(k,i));
                        ++j;
                    } else if (seq(k,i) != alphabet_->gap() && seq(k,i) != alphabet_->endgap()) {
                        out << to_insert_chr(chr(k,i));
                        ++j;
                    }
                    break;
                default:
                    throw Exception("Unsupported alignment output format %i!", format);
            }
            if (j % width == 0) out << std::endl;
        }
        if (j % width != 0) out << std::endl;
    }
}

void Alignment::write_clustal_flavors(std::ostream& out, Format format, int width) const
{
    const int kHeaderLength = 18;

    // convert alignment to character representation
    std::vector<std::string> seqs(nseqs(), "");
    for (int k = 0; k < nseqs(); ++k) {
        for (int i = 0; i < ncols(); ++i) {
            char c = alphabet_->itoc(seqs_[i][k]);
            if (c != '-' && !match_column_[i] && format == PSI) c = to_insert_chr(c);
            seqs[k].push_back(c);
        }
    }

    // print alignment in blocks
    if (format == CLUSTAL) out << "CLUSTAL\n\n";
    while (!seqs.front().empty()) {
        for (int k = 0; k < nseqs(); ++k) {
            std::string header(headers_[k].substr(0, headers_[k].find_first_of(' ')));
            if (static_cast<int>(header.length()) <= kHeaderLength)
                out << header + std::string(kHeaderLength - header.length(), ' ') << ' ';
            else
                out << header.substr(0, kHeaderLength) << ' ';

            out << seqs[k].substr(0, std::min(width, static_cast<int>(seqs[k].length()))) << std::endl;
            seqs[k].erase(0, std::min(width, static_cast<int>(seqs[k].length())));
        }
        out << std::endl; // blank line after each block
    }
}

void Alignment::resize(int nseqs, int ncols)
{
    if (nseqs == 0 || ncols == 0)
        throw Exception("Bad dimensions for alignment resizing: nseqs=%i ncols=%i", nseqs, ncols);

    seqs_.resize(ncols, nseqs);
    column_indexes_.resize(ncols);
    match_indexes.resize(ncols);
    match_column_.resize(ncols);
    headers_.resize(nseqs);
}

void Alignment::assign_match_columns_by_sequence(int k)
{
    for (int i = 0; i < ncols(); ++i)
        match_column_[i] = seqs_[i][k] < alphabet_->gap();
    set_match_indexes();
}

void Alignment::assign_match_columns_by_gap_rule(int gap_threshold)
{
    LOG(DEBUG) << "Marking columns with more than " << gap_threshold << "% of gaps as insert columns ...";

    // global weights are sufficient for calculation of gap percentage
    std::vector<float> wg;
    global_weights_and_diversity(*this, wg);
    for (int i = 0; i < ncols(); ++i) {
        float gap = 0.0f;
        float res = 0.0f;

        for (int k = 0; k < nseqs(); ++k)
            if (seqs_[i][k] < alphabet_->any())
                res += wg[k];
            else
                gap += wg[k];

        float percent_gaps = 100.0f * gap / (res + gap);

        LOG(DEBUG1) << "percent gaps[" << i << "]=" << percent_gaps;

        match_column_[i] = percent_gaps <= static_cast<float>(gap_threshold);
    }
    set_match_indexes();

    LOG(DEBUG) << "ncols=" << ncols() << "  nmatch=" << nmatch() << "  ninsert=" << ninsert();
}

float global_weights_and_diversity(const Alignment& alignment, std::vector<float>& wg)
{
    const float kZero = 1E-10;  // for calculation of entropy
    const int nseqs = alignment.nseqs();
    const int ncols = alignment.nmatch();
    const int nalph = alignment.alphabet()->size();
    const int any   = alignment.alphabet()->any();

    // Return values
    wg.assign(nseqs, 0.0f);  // global sequence weights
    float neff = 0.0f;       // diversity of alignment

    // Helper variables
    std::vector<int> n(nseqs, 0);         // number of residues in sequence i (excl. ANY)
    std::vector<float> fj(nalph, 0.0f);   // to calculate entropy
    std::vector<int> adiff(ncols, 0);     // number of different alphabet letters in each column
    matrix<int> counts(ncols, nalph, 0);  // counts of alphabet letters in each column (excl. ANY)

    LOG(INFO) << "Calculation of global weights and alignment diversity ...";

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
                wg[k] += 1.0f/( adiff[i] * counts[i][alignment[i][k]] * n[k] );
    normalize_to_one(&wg[0], nseqs);

    // Calculate number of effective sequences
    for (int i = 0; i < ncols; ++i) {
        reset(&fj[0], nalph);
        for (int k=0; k < nseqs; ++k)
            if (alignment[i][k] < any) fj[alignment[i][k]] += wg[k];
        normalize_to_one(&fj[0], nalph);
        for (int a = 0; a < nalph; ++a)
            if (fj[a] > kZero) neff -= fj[a] * log2(fj[a]);
    }
    neff = pow(2.0, neff / ncols);

    LOG(DEBUG) << "neff=" << neff;

    return neff;
}

std::vector<float> position_specific_weights_and_diversity(const Alignment& alignment, matrix<float>& w)
{
    const float kMaxEndgapFraction = 0.1;  // maximal fraction of sequences with an endgap
    const int kMinNcols = 10;  // minimum number of columns in subalignments
    const float kZero = 1E-10;  // for calculation of entropy
    const int nseqs  = alignment.nseqs();
    const int ncols  = alignment.nmatch();
    const int nalph  = alignment.alphabet()->size();
    const int any    = alignment.alphabet()->any();
    const int endgap = alignment.alphabet()->endgap();

    // Return values
    std::vector<float> neff(ncols, 0.0f);  // diversity of subalignment i
    w.resize(ncols, nseqs);                // w(i,k) weight of sequence k in column i, calculated from subalignment i
    w = matrix<float>(ncols, nseqs, 0.0f); // init to zero

    // Helper variables
    int ncoli = 0;        // number of columns j that contribute to neff[i]
    int nseqi = 0;        // number of sequences in subalignment i
    int ndiff = 0;        // number of different alphabet letters
    bool change = false;  // has the set of sequences in subalignment changed?

    matrix<int> n(ncols, endgap + 1, 0);     // n(j,a) = number of seq's with some residue at column i AND a at position j
    std::vector<float> fj(nalph, 0.0f);      // to calculate entropy
    std::vector<float> wi(nseqs, 0.0f);      // weight of sequence k in column i, calculated from subalignment i
    std::vector<float> wg;                   // global weights
    std::vector<int> nseqi_debug(ncols, 0);  // debugging
    std::vector<int> ncoli_debug(ncols, 0);  // debugging

    global_weights_and_diversity(alignment, wg);  // compute global sequence weights

    LOG(INFO) << "Calculation of position-specific weights and alignment diversity on subalignments ...";
    LOG(DEBUG1) << "i      ncol   nseq   neff";

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
        nseqi_debug[i] = nseqi;

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
                        if (n[j][alignment[j][k]] == 0) {
                            LOG(WARNING) << "Mi=" << i << ": n[" << j << "][ali[" << j << "][" << k << "]=0!";
                        }
                        wi[k] += 1.0f / static_cast<float>((n[j][alignment[j][k]] * ndiff));
                    }
                }
            }  // for j over ncols
            normalize_to_one(&wi[0], nseqs);

            if (ncoli < kMinNcols)  // number of columns in subalignment insufficient?
                for (int k = 0; k < nseqs; ++k)
                    if (alignment[i][k] < any)
                        wi[k] = wg[k];
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
        ncoli_debug[i] = ncoli;

        LOG(DEBUG1) << left << setw(7) << i << setw(7) << ncoli_debug[i] << setw(7) << nseqi_debug[i] << setw(7) << neff[i];
    }  // for i over ncols

    return neff;
}

}  // cs
