// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "amino_acid.h"
#include "csbuild_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSBuildApp<cs::AminoAcid>().main(argc, argv, stdout, "csbuild");
}
