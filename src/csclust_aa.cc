// Copyright 2009, Andreas Biegert

#include <iostream>

#include "amino_acid.h"
#include "csclust.h"
#include "getopt_pp.h"
#include "log.h"

using namespace GetOpt;

int main(int argc, char* argv[]) {
  GetOpt_pp options(argc, argv);
  cs::Params params;

  if (argc < 2 || options >> OptionPresent('h', "help")) {
    cs::usage<cs::AminoAcid>(params, std::cout);
    return 2;
  }
  try {
    params.parse_options(options);
    cs::csclust<cs::AminoAcid>(params);
  } catch(const std::exception& e) {
    LOG(ERROR) << e.what();
    std::cout << e.what() << std::endl;
    return 3;
  }
  return 0;
}
