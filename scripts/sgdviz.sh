#!/bin/bash

source $HOME/lib/utils.sh

if [ $# -eq 0 ]; then
  echo "sgdviz.sh SGDOUT+"
  exit 1;
fi

MODELS=$@
PARAMS="ll-train ll-val"
for M in $MODELS; do
  OUTDIR=`printf '%s/plots' $(dirname $M)`
  mkdir -p $OUTDIR
  sgdviz.pl -i $M -o `printf '%s/%s.pdf' $OUTDIR $(filename $M)` -p $PARAMS
done
