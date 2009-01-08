#ifndef CS_SEQUENCE_H
#define CS_SEQUENCE_H
/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

// DESCRIPTION:
// A container class representing a sequence consisting of letters over a
// sequence alphabet.

#include <string>
#include <vector>

#include "sequence_alphabet.h"

namespace cs
{

class Sequence
{
public:
    friend std::istream& operator>> (std::istream& i, Sequence& sequence);
    friend std::ostream& operator<< (std::ostream& o, const Sequence& sequence);

    enum Format
    {
        PLAIN     = 0,
        FASTA     = 1,
        CLUSTAL   = 2  //not implemented
    };

    Sequence(int length,
             const SequenceAlphabet* alphabet);
    Sequence(const std::string header,
             const std::string sequence,
             const SequenceAlphabet* alphabet);
    Sequence(const char* header,
             const char* sequence,
             const SequenceAlphabet* alphabet);
    virtual ~Sequence();

    char&       operator() (int i);
    const char& operator() (int i) const;
    char chr(int i) const;
    int length() const;
    const SequenceAlphabet& alphabet() const;

private:
    void init(const char* header, const char* sequence);

    const SequenceAlphabet* const alphabet_;
    std::string header_;
    std::vector<char> sequence_;
};//Sequence



std::istream& operator>> (std::istream& i, Sequence& sequence);

std::ostream& operator<< (std::ostream& o, const Sequence& sequence);


inline char& Sequence::operator() (int i)
{ return sequence_[i]; }

inline const char& Sequence::operator() (int i) const
{ return sequence_[i]; }

inline char Sequence::chr(int i) const
{ return alphabet_->itoc(sequence_[i]); }

inline int Sequence::length() const
{ return sequence_.size(); }

inline const SequenceAlphabet& Sequence::alphabet() const
{ return *alphabet_; }

}//cs

#endif
