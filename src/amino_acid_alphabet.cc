/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "amino_acid_alphabet.h"

namespace cs
{

const char AminoAcidAlphabet::amino_acids_[] = "ARNDCQEGHILKMFPSTWYV";

AminoAcidAlphabet::AminoAcidAlphabet()
        : SequenceAlphabet(20, 'X')
{
    init();
    // Conversion of non-standard amino acids to integers
    set_ctoi('O', any());
    set_ctoi('J', any());
    set_ctoi('U', ctoi('C'));
    set_ctoi('B', ctoi('D'));
    set_ctoi('Z', ctoi('E'));
}

AminoAcidAlphabet* AminoAcidAlphabet::instance()
{
    static AminoAcidAlphabet inst;
    return &inst;
}

}//CS
