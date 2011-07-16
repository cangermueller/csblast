#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -eq 0 ]; then
  echo "crfviz.sh MODEL [ARGS]"
  exit 1
fi

MODEL=$1
shift
ARGS=$@
S="u"
for A in $ARGS; do
  if [ $A == "-s" ]; then 
    S="s"
  fi
done

OUTDIR=`dirname $MODEL`/plots
OUTBASE=$OUTDIR/`basename $MODEL`_$S
mkdir -p $OUTDIR
crfviz -i $MODEL -o $OUTBASE $ARGS
