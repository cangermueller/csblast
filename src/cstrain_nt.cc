/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include <iostream>

#include "nucleotide.h"
#include "cstrain.h"
#include "log.h"

int main(int argc, char* argv[])
{
    GetOpt_pp options(argc, argv);
    cs::CSTrain<cs::Nucleotide> cstrain;

    if(argc < 2 || options >> OptionPresent('h', "help")) {
        cstrain.usage(std::cout);
        return 2;
    }
    try {
        cstrain.run(options);
    } catch(const std::exception& e) {
        LOG(ERROR) << e.what();
        std::cout << e.what() << std::endl;
        return 3;
    }
    return 0;
}
