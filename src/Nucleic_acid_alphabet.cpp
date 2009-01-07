/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include "Nucleic_acid_alphabet.h"


Nucleic_acid_alphabet::Nucleic_acid_alphabet()
{
    init();
}


Nucleic_acid_alphabet::~Nucleic_acid_alphabet()
{ }


Nucleic_acid_alphabet* Nucleic_acid_alphabet::instance()
{
    static Nucleic_acid_alphabet inst;
    return &inst;
}


std::vector<char> Nucleic_acid_alphabet::itoc() const
{
    char itoc_arr[] = {'A','C','G','T','N'};
    std::vector<char> itoc(itoc_arr, itoc_arr + sizeof(itoc_arr)/sizeof(*itoc_arr));
    return itoc;
}
