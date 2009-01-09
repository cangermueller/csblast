/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "sequence.h"

namespace cs
{

Sequence::Sequence(int length,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet),
          sequence_(length)
{}

Sequence::Sequence(const std::string& header,
                   const std::string& sequence,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ init(header, sequence); }

Sequence::Sequence(std::istream& in,
                   const SequenceAlphabet* alphabet)
        : alphabet_(alphabet)
{ in >> *this; }

Sequence::~Sequence()
{}

std::vector<Sequence*> Sequence::read(std::istream& in,
                                      const SequenceAlphabet* alphabet)
{
    std::vector<Sequence*> sequences;
    std::string data((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());

    size_t index = 0;
    std::pair<std::string, std::string > header_seq_pair = parse_fasta_sequence(data, index);
    Sequence* seq = new Sequence(header_seq_pair.first,
                                 header_seq_pair.second,
                                 alphabet);
    sequences.push_back(seq);

    return sequences;
}

// Initializes sequence object with given header string and sequence vector.
// Note that the sequence undergoes a validation check: White space is silently
// removed whereas invalid characters (according to the alphabet) trigger an
// execption.
void Sequence::init(const std::string& header, const std::string& sequence)
{
    std::string().swap(header_); //clear and minimize capacity
    header_.insert(header_.begin(), header.begin() + (header[0]=='>' ? 1 : 0), header.end());

    std::vector<char>().swap(sequence_); //clear and minimize capacity
    sequence_.reserve(sequence.size());
    sequence_.insert(sequence_.begin(), sequence.begin(), sequence.end());
    //strip whitespace and newlines from sequence.
    sequence_.erase(remove_if(sequence_.begin(), sequence_.end(), isspace),
                       sequence_.end());

    //validate each character and convert to integer representation
    int len = sequence_.size();
    for (int i = 0; i < len; ++i)
        if (alphabet_->valid(sequence_[i])) {
            sequence_[i] = alphabet_->ctoi(sequence_[i]);
        } else {
            char c = sequence_[i];
            sequence_.clear();
            throw MyException("Invalid character %c at position %i of sequence '%s'", c, i, header_.c_str());
        }
}

std::pair<std::string, std::string > parse_fasta_sequence(const std::string& data, size_t& index)
{
    size_t i = data.find_first_of('>', index); //find start of header
    if (i == std::string::npos) throw MyException("Bad format: FASTA input does not contain '>'!");
    size_t j = data.find_first_of('\n', i); //find end of header
    if (j == std::string::npos) throw MyException("Bad format: FASTA header does not terminate with newline!");
    size_t k = data.find_first_of('>', j); //find start of next sequence if there is one
    index = k == std::string::npos ? 0 : k;

    std::string header(data.begin()+i, data.begin()+j);
    std::string raw_sequence(data.begin()+j, k == std::string::npos ? data.end() : data.begin()+k);

    return std::make_pair(header, raw_sequence);
}

std::istream& operator>> (std::istream& in, Sequence& sequence)
{
    std::string data((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());

    size_t index = 0;
    std::pair<std::string, std::string > header_seq_pair = parse_fasta_sequence(data, index);
    sequence.init(header_seq_pair.first, header_seq_pair.second);

    return in;
}

std::ostream& operator<< (std::ostream& out, const Sequence& sequence)
{
    const int kLineLength = 80;

    out << '>' << sequence.header() << std::endl;
    int len = sequence.length();
    for (int i = 0; i < len; ++i) {
        out << sequence.chr(i);
        if ((i+1) % kLineLength == 0) out << std::endl;
    }

    return out;
}

}//cs
