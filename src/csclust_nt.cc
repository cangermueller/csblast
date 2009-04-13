// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "nucleotide.h"
#include "csclust_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSClustApp<cs::Nucleotide>().main(argc, argv, stdout, "csclust");
}
