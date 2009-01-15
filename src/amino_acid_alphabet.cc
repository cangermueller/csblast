/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "amino_acid_alphabet.h"

namespace cs
{

AminoAcidAlphabet::AminoAcidAlphabet()
{
    init();
    // Conversion of non-standard amino acids to integers
    ctoi_['J'] = any();
    ctoi_['O'] = any();
    ctoi_['U'] = ctoi_['C'];
    ctoi_['B'] = ctoi_['D'];
    ctoi_['Z'] = ctoi_['E'];
}

AminoAcidAlphabet::~AminoAcidAlphabet()
{}

AminoAcidAlphabet* AminoAcidAlphabet::instance()
{
    static AminoAcidAlphabet inst;
    return &inst;
}

void AminoAcidAlphabet::init_itoc()
{
    char itoc_arr[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X'};
    itoc_.insert(itoc_.begin(), itoc_arr, itoc_arr + sizeof(itoc_arr)/sizeof(*itoc_arr));
}

}//CS