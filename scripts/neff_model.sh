#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -ne 2 ]; then
	echo "neff_model.sh MODEL ADMIX"
	exit 1
fi

MODEL=$1
ADMIX=$2
DB="$DBS/scop20_1.75_opt"
OUT_DIR=$CSD/neff/`basename $MODEL`
OUT_FILE=`printf '%s/x%.2f.dat' $OUT_DIR $ADMIX`

mkdir -p $OUT_DIR
rm -f $OUT_FILE
for SEQ in $DB/*seq; do
	NEFF=`neff_seq.sh $MODEL $ADMIX $SEQ`
  if [[ ! $NEFF =~ ^[0-9]+\.[0-9]+$ ]]; then
    echo $NEFF
    echo $SEQ
    exit 0
  else
    echo $NEFF >> $OUT_FILE
  fi
done
