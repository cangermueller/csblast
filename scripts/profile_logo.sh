#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -lt 2 ]; then
    echo "Missing arguments!"
	echo "profile_logo.sh MODEL SEQUENCE_FILE+"
	exit 1
fi

MODEL=$1
shift
SEQS=$@
OUT_DIR=$CSD/profile_logos

for SEQ in $SEQS; do
	OUT_BASE=`basename $SEQ`
	OUT_BASE=$OUT_DIR/${OUT_BASE%.*}
	csbuild -i $SEQ -D $MODEL -o $OUT_BASE.prf
	csviz -i $OUT_BASE.prf -o $OUT_BASE.pdf
done

