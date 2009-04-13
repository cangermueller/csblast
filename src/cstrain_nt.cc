// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "nucleotide.h"
#include "cstrain_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSTrainApp<cs::Nucleotide>().main(argc, argv, stdout, "cstrain");
}
