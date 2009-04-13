// Copyright 2009, Andreas Biegert

#include <cstdio>

#include "nucleotide.h"
#include "cssample_app.cc"

int main(int argc, char* argv[]) {
  return cs::CSSampleApp<cs::Nucleotide>().main(argc, argv, stdout, "cssample");
}
