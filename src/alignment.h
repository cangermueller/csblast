#ifndef CS_ALIGNMENT_H
#define CS_ALIGNMENT_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class for an alignment of sequences over an alphabet.

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
#include "shared_ptr.h"
#include "util.h"

using std::setw;
using std::left;

namespace
{

// Converts a character to uppercase and '.' to '-'.
inline char to_match_chr(char c)
{
    return isalpha(c) ? toupper(c) : (c == '.' ? '-' : c);
}

// Converts a character to lowercase and '-' to '.'.
inline char to_insert_chr(char c)
{
    return isalpha(c) ? tolower(c) : (c == '-' ? '.' : c);
}

// Predicate indicating if character belongs to match column.
inline bool match_chr(char c)
{
    return (isalpha(c) && isupper(c)) || c == '-';
}

// Predicate indicating if character belongs to insert column.
inline char insert_chr(char c)
{
    return (isalpha(c) && islower(c)) || c == '.';
}

} // namespace

namespace cs
{

template<class Alphabet_T>
class Alignment
{
  public:
    // Supported alignment formats for in- and output.
    enum Format {
        FASTA    = 0,
        A2M      = 1,
        A3M      = 2,
        CLUSTAL  = 3,
        PSI      = 4
    };

    // Public typedefs
    typedef matrix<char>::row_type col_type;
    typedef matrix<char>::const_row_type const_col_type;
    typedef matrix<char>::col_type row_type;
    typedef matrix<char>::const_col_type const_row_type;
    typedef matrix<char>::iterator iterator;
    typedef matrix<char>::const_iterator const_iterator;

    // Constructs alignment multi FASTA formatted alignment read from input stream.
    Alignment(std::istream& in, Format format);

    ~Alignment() {}

    // Reads all available alignments from the input stream and returns them in a vector.
    static void readall(std::istream& in, Format format, std::vector< shared_ptr<Alignment> >& v);

    // Access methods to get the integer representation of character in match column i of sequence k.
    col_type operator[](int i) { return seqs_[match_indexes[i]]; }
    const_col_type operator[](int i) const { return seqs_[match_indexes[i]]; }
    // Returns the integer in column i of sequence k.
    char seq(int k, int i) const { return seqs_[i][k]; }
    // Returns the character in column i of sequence k.
    char chr(int k, int i) const { return Alphabet_T::instance().itoc(seqs_[i][k]); }
    // Returns the number of sequences in the alignment.
    int num_seqs() const { return seqs_.num_cols(); }
    // Returns the total number of alignment columns.
    int num_cols() const { return seqs_.num_rows(); }
    // Returns the number of match columns.
    int num_match_cols() const { return match_indexes.size(); }
    // Returns the number of insert columns.
    int num_insert_cols() const { return num_cols() - num_match_cols(); }
    // Returns the total number of characters in the alignment (incl. inserts).
    int size() const { return seqs_.size(); }

    row_type seq_begin(int k) { return seqs_.col_begin(k); }
    // Returns an iterator just past the end of alignment sequence k.
    row_type seq_end(int k) { return seqs_.col_end(k); }
    // // Returns a const iterator to the first element in alignment sequence k.
    const_row_type seq_begin(int k) const { return seqs_.col_begin(k); }
    // // Returns a const iterator just past the end of alignment sequence k.
    const_row_type seq_end(int k) const { return seqs_.col_end(k); }

    // Returns an iterator to the first element in the alignment matrix (incl. inserts).
    iterator begin() { return seqs_.begin(); }
    // Returns an iterator just past the end of the alignment matrix (incl. inserts).
    iterator end() { return seqs_.end(); }
    // Returns a const iterator to the first element in the alignment matrix (incl. inserts).
    const_iterator begin() const { return seqs_.begin(); }
    // Returns a const iterator just past the end of the alignment matrix (incl. inserts).
    const_iterator end() const { return seqs_.end(); }

    // Returns the header of sequence k.
    std::string header(int k) const { return headers_[k]; }
    // Sets the header of sequence k.
    void set_header(int k, const std::string& header) { headers_[k] = header; }
    // Makes all columns with a residue in sequence k match columns.
    void assign_match_columns_by_sequence(int k = 0);
    // Makes all columns with less than X% gaps match columns.
    void assign_match_columns_by_gap_rule(int gap_threshold = 50);
    // Initializes the alignment object with an alignment in FASTA format read from given stream.
    void read(std::istream& in, Format format);
    // Writes the alignment in given format to ouput stream.
    void write(std::ostream& out, Format format, int width = 100) const;
    // Returns true if column i is a match column.
    bool match_column(int i) const { return match_column_[i]; }
    // Removes all insert columns from the alignment.
    void remove_insert_columns();

    // Prints the Alignment in A2M format for debugging.
    friend std::ostream& operator<< (std::ostream& out, const Alignment& alignment)
    {
        alignment.write(out, Alignment::CLUSTAL);
        return out;
    }

  private:
    // Disallow copy and assign
    Alignment(const Alignment&);
    void operator=(const Alignment&);

    // Initializes alignment with given headers and sequences.
    void init(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);
    // Resize the sequence matrix and header vector to given dimensions. Attention: old data is lost!
    void resize(int num_seqs, int num_cols);
    // Fills match_indexes_ with the indexes of all match columns.
    void set_match_indexes();
    // Reads an alignment in FASTA format.
    void read_fasta(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Reads an alignment in A2M format from given stream.
    void read_a2m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Reads an alignment in A3M format from given stream.
    void read_a3m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Helper method that reads a FASTA, A2M, or A3M formatted alignment.
    void read_fasta_flavors(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs);
    // Writes the alignment in FASTA, A2M, or A3M format to output stream.
    void write_fasta_flavors(std::ostream& out, Format format, int width = 100) const;
    // Writes the alignment in CLUSTAL or PSI format to output stream.
    void write_clustal_flavors(std::ostream& out, Format format, int width = 100) const;

    // Row major matrix with sequences in integer representation.
    matrix<char> seqs_;
    // Array with indices of all columns [0,1,2,...,num_cols-1].
    std::valarray<int> column_indexes_;
    // Array with indices of match columns.
    std::valarray<int> match_indexes;
    // Array mask indicating match and insert columns with true and false respectively.
    std::valarray<bool> match_column_;
    // Headers of sequences in the alignment.
    std::vector<std::string> headers_;
};  // Alignment



// Calculates global sequence weights by maximum entropy weighting (Henikoff&Henikoff '94).
template<class Alphabet_T>
float global_weights_and_diversity(const Alignment<Alphabet_T>& alignment, std::vector<float>& wg);

// Calculates position-dependent sequence weights and number of effective sequences on subalignments.
template<class Alphabet_T>
std::vector<float> position_specific_weights_and_diversity(const Alignment<Alphabet_T>& alignment, matrix<float>& w);

// Returns the alignment format corresponding to provided filename extension
template<class Alphabet_T>
typename Alignment<Alphabet_T>::Format alignment_format_from_string(const std::string& s);



template<class Alphabet_T>
inline Alignment<Alphabet_T>::Alignment(std::istream& in, Format format)
{
    read(in, format);
}

template<class Alphabet_T>
inline void Alignment<Alphabet_T>::readall(std::istream& in, Format format, std::vector< shared_ptr<Alignment> >& v)
{
    while (in.peek() && in.good()) {
        shared_ptr<Alignment> p(new Alignment(in, format));
        v.push_back(p);
    }
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::init(const std::vector<std::string>& headers, const std::vector<std::string>& seqs)
{
    if (seqs.empty()) throw Exception("Unable to initialize alignment: no aligned sequences found!");
    if (headers.size() != seqs.size())
        throw Exception("Unable to initialize alignment: unequal number of headers and sequences!");

    const int num_seqs = seqs.size();
    const int num_cols = seqs[0].length();
    for (int k = 1; k < num_seqs; ++k)
        if (static_cast<int>(seqs[k].length()) != num_cols)
            throw Exception("Bad alignment: sequence %i has length %i but should have %i!", k+1, seqs[k].length(), num_cols);

    // Validate characters and convert to integer representation
    const Alphabet_T& alphabet = Alphabet_T::instance();
    headers_.resize(num_seqs);
    resize(seqs.size(), seqs[0].length());
    for (int k = 0; k < num_seqs; ++k) {
        headers_[k] = headers[k];
        for (int i = 0; i < num_cols; ++i) {
            const char c = seqs[k][i];
            column_indexes_[i] = i;
            match_column_[i] = match_chr(c);
            if (alphabet.valid(c, true))
                seqs_[i][k] = alphabet.ctoi(c);
            else
                throw Exception("Invalid character %c at position %i of sequence '%s'", c, i, headers_[k].c_str());
        }
    }

    // Replace gap with endgap for all gaps at either end of a sequence
    for (int k = 0; k < num_seqs; ++k) {
        for (int i = 0; i < num_cols && seqs_[i][k] == alphabet.gap(); ++i)
            seqs_[i][k] = alphabet.endgap();
        for (int i = num_cols - 1; i >= 0 && seqs_[i][k] == alphabet.gap(); --i)
            seqs_[i][k] = alphabet.endgap();
    }

    // Initialize index array for match columns
    set_match_indexes();
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::set_match_indexes()
{
    const int match_cols = std::count(&match_column_[0], &match_column_[0] + num_cols(), true);
    match_indexes.resize(match_cols);
    match_indexes = column_indexes_[match_column_];
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::read(std::istream& in, Format format)
{
    LOG(DEBUG4) << "Reading alignment from stream ...";

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

    LOG(DEBUG4) << *this;
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::read_fasta_flavors(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    headers.clear();
    seqs.clear();
    std::string buffer;

    while (in.peek() != '>' && in.good()) getline(in, buffer);  // advance to first sequence
    while (in.good()) {
        // Terminator encountered
        if (buffer[0] == '#') break;

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
            if (buffer.empty()) continue;
            if (buffer[0] == '#') break;
            seqs.back().append(buffer.begin(), buffer.end());
        }
        // Remove whitespace from sequence
        seqs.back().erase(remove_if(seqs.back().begin(), seqs.back().end(), isspace), seqs.back().end());
        LOG(DEBUG) << headers.back();
        LOG(DEBUG) << "seqlen=" << seqs.back().length();
    }
    if (headers.empty()) throw Exception("Bad alignment input: no alignment data found in stream!");
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::read_fasta(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);
    // convert all characters to match characters
    for (std::vector<std::string>::iterator it = seqs.begin(); it != seqs.end(); ++it)
        transform(it->begin(), it->end(), it->begin(), to_match_chr);
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::read_a2m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::read_a3m(std::istream& in, std::vector<std::string>& headers, std::vector<std::string>& seqs)
{
    read_fasta_flavors(in, headers, seqs);

    // Check number of match states
    const int num_seqs = seqs.size();
    const int num_match_cols = count_if(seqs[0].begin(), seqs[0].end(), match_chr);
    for (int k = 1; k < num_seqs; ++k) {
        const int num_match_cols_k = count_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (num_match_cols_k != num_match_cols)
            throw Exception("Bad alignment: sequence %i has %i match columns but should have %i!", k, num_match_cols_k, num_match_cols);
        if (count(seqs[k].begin(), seqs[k].end(), '.') > 0)
            throw Exception("Bad alignment: sequence %i in A3M formatted alignment contains gap characters '.'!", k);
    }

    // Insert gaps into A3M alignment
    std::vector<std::string> seqs_a2m(seqs.size(), "");     // sequences in A2M format
    matrix<std::string> inserts(seqs.size(), num_match_cols, "");   // inserts of sequence k after match state i
    std::vector<int> max_insert_len(num_match_cols, 0);             // maximal length of inserts after match state i

    // Move inserts BEFORE first match state into seqs_a2m
    for (int k = 0; k < num_seqs; ++k) {
        std::string::iterator i = find_if(seqs[k].begin(), seqs[k].end(), match_chr);
        if (i != seqs[k].end()) {
            seqs_a2m[k].append(seqs[k].begin(), i);
            seqs[k].erase(seqs[k].begin(), i);
        }
    }

    // Extract all inserts and keep track of longest insert after each match column
    for (int k = 0; k < num_seqs; ++k) {
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
    for (int k = 0; k < num_seqs; ++k) {
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

template<class Alphabet_T>
void Alignment<Alphabet_T>::write(std::ostream& out, Format format, int width) const
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

template<class Alphabet_T>
void Alignment<Alphabet_T>::write_fasta_flavors(std::ostream& out, Format format, int width) const
{
    for (int k = 0; k < num_seqs(); ++k) {
        out << '>' << headers_[k] << std::endl;
        int j = 0;  // counts printed characters
        for (int i = 0; i < num_cols(); ++i) {
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
                    } else if (seq(k,i) != Alphabet_T::instance().gap() && seq(k,i) != Alphabet_T::instance().endgap()) {
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

template<class Alphabet_T>
void Alignment<Alphabet_T>::write_clustal_flavors(std::ostream& out, Format format, int width) const
{
    const int HEADER_WIDTH = 18;

    // convert alignment to character representation
    std::vector<std::string> seqs(num_seqs(), "");
    for (int k = 0; k < num_seqs(); ++k) {
        for (int i = 0; i < num_cols(); ++i) {
            char c = Alphabet_T::instance().itoc(seqs_[i][k]);
            if (c != '-' && !match_column_[i] && format == PSI) c = to_insert_chr(c);
            seqs[k].push_back(c);
        }
    }

    // print alignment in blocks
    if (format == CLUSTAL) out << "CLUSTAL\n\n";
    while (!seqs.front().empty()) {
        for (int k = 0; k < num_seqs(); ++k) {
            std::string header(headers_[k].substr(0, headers_[k].find_first_of(' ')));
            if (static_cast<int>(header.length()) <= HEADER_WIDTH)
                out << header + std::string(HEADER_WIDTH - header.length(), ' ') << ' ';
            else
                out << header.substr(0, HEADER_WIDTH) << ' ';

            out << seqs[k].substr(0, std::min(width, static_cast<int>(seqs[k].length()))) << std::endl;
            seqs[k].erase(0, std::min(width, static_cast<int>(seqs[k].length())));
        }
        out << std::endl; // blank line after each block
    }
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::resize(int num_seqs, int num_cols)
{
    if (num_seqs == 0 || num_cols == 0)
        throw Exception("Bad dimensions for alignment resizing: num_seqs=%i num_cols=%i", num_seqs, num_cols);

    seqs_.resize(num_cols, num_seqs);
    column_indexes_.resize(num_cols);
    match_indexes.resize(num_cols);
    match_column_.resize(num_cols);
    headers_.resize(num_seqs);
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::assign_match_columns_by_sequence(int k)
{
    for (int i = 0; i < num_cols(); ++i)
        match_column_[i] = seqs_[i][k] < Alphabet_T::instance().gap();
    set_match_indexes();
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::assign_match_columns_by_gap_rule(int gap_threshold)
{
    LOG(DEBUG) << "Marking columns with more than " << gap_threshold << "% of gaps as insert columns ...";

    // global weights are sufficient for calculation of gap percentage
    std::vector<float> wg;
    global_weights_and_diversity(*this, wg);
    for (int i = 0; i < num_cols(); ++i) {
        float gap = 0.0f;
        float res = 0.0f;

        for (int k = 0; k < num_seqs(); ++k)
            if (seqs_[i][k] < Alphabet_T::instance().any())
                res += wg[k];
            else
                gap += wg[k];

        float percent_gaps = 100.0f * gap / (res + gap);

        LOG(DEBUG3) << "percent gaps[" << i << "]=" << percent_gaps;

        match_column_[i] = percent_gaps <= static_cast<float>(gap_threshold);
    }
    set_match_indexes();

    LOG(DEBUG) << "num_cols=" << num_cols() << "  num_match_cols=" << num_match_cols() << "  num_insert_cols=" << num_insert_cols();
}

template<class Alphabet_T>
void Alignment<Alphabet_T>::remove_insert_columns()
{
    // create new sequence matrix
    const int match_cols = num_match_cols();
    matrix<char> new_seqs(match_cols, num_seqs());
    for (int i = 0; i < match_cols; ++i) {
        for (int k = 0; k < num_seqs(); ++k) {
            new_seqs[i][k] = seqs_[match_indexes[i]][k];
        }
    }
    seqs_.resize(match_cols, num_seqs());
    seqs_ = new_seqs;

    // update match indexes
    column_indexes_.resize(match_cols);
    match_indexes.resize(match_cols);
    match_column_.resize(match_cols);
    for (int i = 0; i < match_cols; ++i) {
        column_indexes_[i] = i;
        match_column_[i]   = true;
    }
    set_match_indexes();
}



template<class Alphabet_T>
float global_weights_and_diversity(const Alignment<Alphabet_T>& alignment, std::vector<float>& wg)
{
    const float ZERO = 1E-10;  // for calculation of entropy
    const int num_seqs = alignment.num_seqs();
    const int num_cols = alignment.num_match_cols();
    const int alphabet_size = Alphabet_T::instance().size();
    const int any   = Alphabet_T::instance().any();

    // Return values
    wg.assign(num_seqs, 0.0f);  // global sequence weights
    float neff = 0.0f;          // diversity of alignment

    // Helper variables
    std::vector<int> n(num_seqs, 0);                 // number of residues in sequence i (excl. ANY)
    std::vector<float> fj(alphabet_size, 0.0f);      // to calculate entropy
    std::vector<int> adiff(num_cols, 0);             // number of different alphabet letters in each column
    matrix<int> counts(num_cols, alphabet_size, 0);  // counts of alphabet letters in each column (excl. ANY)

    LOG(INFO) << "Calculation of global weights and alignment diversity ...";

    // Count number of residues in each column
    for (int i = 0; i < num_cols; ++i)
        for (int k = 0; k < num_seqs; ++k)
            if (alignment[i][k] < any) {
                ++counts[i][alignment[i][k]];
                ++n[k];
            }

    // Count number of different residues in each column
    for (int i = 0; i < num_cols; ++i)
        for (int a = 0; a < alphabet_size; ++a)
            if (counts[i][a]) ++adiff[i];

    // Calculate weights
    for(int i = 0; i < num_cols; ++i)
        for (int k = 0; k < num_seqs; ++k)
            if( adiff[i] > 0 && alignment[i][k] < any )
                wg[k] += 1.0f/( adiff[i] * counts[i][alignment[i][k]] * n[k] );
    normalize_to_one(&wg[0], num_seqs);

    // Calculate number of effective sequences
    for (int i = 0; i < num_cols; ++i) {
        reset(&fj[0], alphabet_size);
        for (int k=0; k < num_seqs; ++k)
            if (alignment[i][k] < any) fj[alignment[i][k]] += wg[k];
        normalize_to_one(&fj[0], alphabet_size);
        for (int a = 0; a < alphabet_size; ++a)
            if (fj[a] > ZERO) neff -= fj[a] * log2(fj[a]);
    }
    neff = pow(2.0, neff / num_cols);

    LOG(DEBUG) << "neff=" << neff;

    return neff;
}

template<class Alphabet_T>
std::vector<float> position_specific_weights_and_diversity(const Alignment<Alphabet_T>& alignment, matrix<float>& w)
{
    const float MAX_ENDGAP_FRACTION = 0.1;  // maximal fraction of sequences with an endgap
    const int MIN_COLS = 10;   // minimum number of columns in subalignments
    const float ZERO = 1E-10;  // for calculation of entropy
    const int num_seqs  = alignment.num_seqs();
    const int num_cols  = alignment.num_match_cols();
    const int alphabet_size  = Alphabet_T::instance().size();
    const int any    = Alphabet_T::instance().any();
    const int endgap = Alphabet_T::instance().endgap();

    // Return values
    std::vector<float> neff(num_cols, 0.0f);     // diversity of subalignment i
    w.resize(num_cols, num_seqs);                // w(i,k) weight of sequence k in column i, calculated from subalignment i
    w = matrix<float>(num_cols, num_seqs, 0.0f); // init to zero

    // Helper variables
    int ncoli = 0;        // number of columns j that contribute to neff[i]
    int nseqi = 0;        // number of sequences in subalignment i
    int ndiff = 0;        // number of different alphabet letters
    bool change = false;  // has the set of sequences in subalignment changed?

    matrix<int> n(num_cols, endgap + 1, 0);      // n(j,a) = number of seq's with some residue at column i AND a at position j
    std::vector<float> fj(alphabet_size, 0.0f);  // to calculate entropy
    std::vector<float> wi(num_seqs, 0.0f);       // weight of sequence k in column i, calculated from subalignment i
    std::vector<float> wg;                       // global weights
    std::vector<int> nseqi_debug(num_cols, 0);   // debugging
    std::vector<int> ncoli_debug(num_cols, 0);   // debugging

    global_weights_and_diversity(alignment, wg);  // compute global sequence weights

    LOG(INFO) << "Calculation of position-specific weights and alignment diversity on subalignments ...";
    LOG(DEBUG1) << "i      ncol   nseq   neff";

    for (int i = 0; i < num_cols; ++i) {
        change = false;
        for (int k = 0; k < num_seqs; ++k) {
            if ((i==0 || alignment[i-1][k] >= any) && alignment[i][k] < any) {
                change = true;
                ++nseqi;
                for (int j = 0; j < num_cols; ++j)
                    ++n[j][alignment[j][k]];
            } else if (i>0 && alignment[i-1][k] < any && alignment[i][k] >= any) {
                change = true;
                --nseqi;
                for (int j = 0; j < num_cols; ++j)
                    --n[j][alignment[j][k]];
            }
        }  // for k over num_seqs
        nseqi_debug[i] = nseqi;

        if (change) {  // set of sequences in subalignment has changed
            ncoli = 0;
            reset(&wi[0], num_seqs);

            for (int j = 0; j < num_cols; ++j) {
                if (n[j][endgap] > MAX_ENDGAP_FRACTION * nseqi) continue;
                ndiff = 0;
                for (int a = 0; a < alphabet_size; ++a) if (n[j][a]) ++ndiff;
                if (ndiff == 0) continue;
                ++ncoli;
                for (int k = 0; k < num_seqs; ++k) {
                    if (alignment[i][k] < any && alignment[j][k] < any) {
                        if (n[j][alignment[j][k]] == 0) {
                            LOG(WARNING) << "Mi=" << i << ": n[" << j << "][ali[" << j << "][" << k << "]=0!";
                        }
                        wi[k] += 1.0f / static_cast<float>((n[j][alignment[j][k]] * ndiff));
                    }
                }
            }  // for j over num_cols
            normalize_to_one(&wi[0], num_seqs);

            if (ncoli < MIN_COLS)  // number of columns in subalignment insufficient?
                for (int k = 0; k < num_seqs; ++k)
                    if (alignment[i][k] < any)
                        wi[k] = wg[k];
                    else
                        wi[k] = 0.0f;

            neff[i] = 0.0f;
            for (int j = 0; j < num_cols; ++j) {
                if (n[j][endgap] > MAX_ENDGAP_FRACTION * nseqi) continue;
                reset(&fj[0], alphabet_size);

                for (int k = 0; k < num_seqs; ++k)
                    if (alignment[i][k] < any && alignment[j][k] < any)
                        fj[alignment[j][k]] += wi[k];
                normalize_to_one(&fj[0], alphabet_size);
                for (int a = 0; a < alphabet_size; ++a)
                    if (fj[a] > ZERO) neff[i] -= fj[a] * log2(fj[a]);
            }  // for j over num_cols

            neff[i] = (ncoli > 0 ? pow(2.0, neff[i] / ncoli) : 1.0f);

        } else {  // set of sequences in subalignment has NOT changed
            neff[i] = (i==0 ? 0.0 : neff[i-1]);
        }

        for (int k = 0; k < num_seqs; ++k) w[i][k] = wi[k];
        ncoli_debug[i] = ncoli;

        LOG(DEBUG1) << left << setw(7) << i << setw(7) << ncoli_debug[i] << setw(7) << nseqi_debug[i] << setw(7) << neff[i];
    }  // for i over num_cols

    return neff;
}

template<class Alphabet_T>
typename Alignment<Alphabet_T>::Format alignment_format_from_string(const std::string& s)
{
    if (s == "fas")
        return Alignment<Alphabet_T>::FASTA;
    else if (s == "a2m")
        return Alignment<Alphabet_T>::A2M;
    else if (s == "a3m")
        return Alignment<Alphabet_T>::A3M;
    else if (s == "clu")
        return Alignment<Alphabet_T>::CLUSTAL;
    else if (s == "psi")
        return Alignment<Alphabet_T>::PSI;
    else
        throw Exception("Unable to infer alignment format from filename extension '%s'!", s.c_str());
}

}  // cs

#endif
