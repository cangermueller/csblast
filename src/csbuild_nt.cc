// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "nucleotide.h"
#include "csbuild_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSBuildApp<cs::Nucleotide>().main(argc, argv, stdout, "csbuild");
}

