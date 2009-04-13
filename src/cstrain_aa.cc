// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "amino_acid.h"
#include "cstrain_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSTrainApp<cs::AminoAcid>().main(argc, argv, stdout, "cstrain");
}
