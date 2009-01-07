/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "Sequence_alphabet.h"


Sequence_alphabet::Sequence_alphabet()
{ }


Sequence_alphabet::~Sequence_alphabet()
{ }


inline bool Sequence_alphabet::valid(char letter) const
{ return ctoi_[letter] != -1; }


inline size_t Sequence_alphabet::size() const
{ return itoc_.size(); }


inline int Sequence_alphabet::ctoi(char letter) const
{ return ctoi_[letter]; }


inline char Sequence_alphabet::itoc(int letter) const
{ return itoc_[letter]; }


inline int Sequence_alphabet::any() const
{ return ctoi_[itoc_[itoc_.size()-1]]; }


inline Sequence_alphabet::const_iterator Sequence_alphabet::begin() const
{ return itoc_.begin(); }


inline Sequence_alphabet::const_iterator Sequence_alphabet::end() const
{ return itoc_.end(); }


inline void Sequence_alphabet::init()
{
    itoc_ = itoc();
    const size_t char_size = static_cast<int>(pow(2, 8*sizeof(char)));
    for(size_t i=0; i<char_size; ++i) ctoi_.push_back(-1);
    for(int i=0; i<itoc_.size(); ++i) ctoi_[itoc_[i]]=i;
}


std::vector<char> Sequence_alphabet::itoc() const
{
    std::vector<char> dummy_itoc;
    return dummy_itoc;
}



