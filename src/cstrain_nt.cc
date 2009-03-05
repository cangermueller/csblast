/***************************************************************************
 *   Copyright (C) 2008 by Andreas Biegert                                 *
 *   andreas.biegert@googlemail.com                                        *
 ***************************************************************************/

#include <iostream>

#include "nucleotide.h"
#include "cstrain.h"
#include "getopt_pp.h"
#include "log.h"

using namespace GetOpt;

int main(int argc, char* argv[])
{
    GetOpt_pp options(argc, argv);
    cs::Params params;

    if(argc < 2 || options >> OptionPresent('h', "help")) {
        cs::usage<cs::Nucleotide>(params, std::cout);
        return 2;
    }
    try {
        params.parse_options(options);
        cs::cstrain<cs::Nucleotide>(params);
    } catch(const std::exception& e) {
        LOG(ERROR) << e.what();
        std::cout << e.what() << std::endl;
        return 3;
    }
    return 0;
}
