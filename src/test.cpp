/***************************************************************************
 *   Copyright (C) 2006 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/
#include <iostream>
#include <cstddef>
#include <vector>
#include <cmath>

#include "Amino_acid_alphabet.h"
#include "Nucleic_acid_alphabet.h"

using std::cout;
using std::endl;

int main(int argc, const char **argv) {
    Amino_acid_alphabet* aa = Amino_acid_alphabet::instance();
    Nucleic_acid_alphabet* na = Nucleic_acid_alphabet::instance();

    for(Sequence_alphabet::const_iterator iter=aa->begin(); iter != aa->end(); ++iter)
        std::cout << *iter << std::endl;

    for(Sequence_alphabet::const_iterator iter=na->begin(); iter != na->end(); ++iter)
        std::cout << *iter << std::endl;

    return 0;
}
