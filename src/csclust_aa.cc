// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "amino_acid.h"
#include "csclust_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSClustApp<cs::AminoAcid>().main(argc, argv, stdout, "csclust");
}
