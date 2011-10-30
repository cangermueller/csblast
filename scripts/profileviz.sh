#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -lt 1 ]; then
	echo "Model missing!"
	echo "Usage: profileviz.sh MODEL [OUTBASE]"
	exit 1
fi

MODEL=$1
OUTBASE=$2
SEQ=$HOME/seqs/zinc_finger.seq

if [ -z "$OUTBASE" ]; then
  NAME=`basename $M`
  OUTDIR=`dirname $M`/plots/profile
  mkdir -p $OUTDIR
  OUTBASE=$OUTDIR/$NAME
fi
csbuild -i $SEQ -o $OUTBASE.prf -D $MODEL
csviz -i $OUTBASE.prf -o $OUTBASE.pdf
rm $OUTBASE.prf
