#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -lt 1 ]; then
	echo Model missing!
	echo Usage: profileviz.sh MODEL+
	exit 1
fi

MODELS=$@
SEQ=$HOME/seqs/zinc_finger.seq
for M in $MODELS; do
  NAME=`basename $M`
  OUTDIR=`dirname $M`/plots/profile
  mkdir -p $OUTDIR
  OUTBASE=$OUTDIR/$NAME
  csbuild -i $SEQ -o $OUTBASE.prf -D $M
  csviz -i $OUTBASE.prf -o $OUTBASE.pdf
  rm $OUTBASE.prf
done

