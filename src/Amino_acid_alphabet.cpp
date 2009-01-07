/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "Amino_acid_alphabet.h"


Amino_acid_alphabet::Amino_acid_alphabet()
{
    init();
    // Conversion of non-standard amino acids to integers
    ctoi_['J'] = any();
    ctoi_['O'] = any();
    ctoi_['U'] = ctoi_['C'];
    ctoi_['B'] = ctoi_['D'];
    ctoi_['Z'] = ctoi_['E'];
}


Amino_acid_alphabet::~Amino_acid_alphabet()
{ }


std::vector<char> Amino_acid_alphabet::itoc() const
{
    char itoc_arr[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X'};
    std::vector<char> itoc(itoc_arr, itoc_arr + sizeof(itoc_arr)/sizeof(*itoc_arr));
    return itoc;
}
