/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "amino_acid_alphabet.h"

namespace cs
{

const char AminoAcidAlphabet::amino_acids_[] = "ARNDCQEGHILKMFPSTWYVBJZX";

AminoAcidAlphabet::AminoAcidAlphabet()
{
    init();
    // Conversion of non-standard amino acids to integers
    set_ctoi('O', any());
    set_ctoi('U', ctoi('C'));
}

AminoAcidAlphabet* AminoAcidAlphabet::instance()
{
    static AminoAcidAlphabet inst;
    return &inst;
}

}//CS
