#!/bin/bash

source $HOME/lib/utils.sh

if [ $# -eq 0 ]; then
  echo "sgdviz.sh MODEL+"
  exit 1;
fi

MODELS=$@
for M in $MODELS; do
  OUTDIR=`printf '%s/plots' $(dirname $M)`
  mkdir -p $OUTDIR
  sgdviz.pl -i $M -o `printf '%s/%s.pdf' $OUTDIR $(filename $M)`
done
