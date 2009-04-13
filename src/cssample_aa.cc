// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "amino_acid.h"
#include "cssample_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSSampleApp<cs::AminoAcid>().main(argc, argv, stdout, "cssample");
}

