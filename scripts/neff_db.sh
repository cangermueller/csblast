#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -ne 3 ]; then 
  echo "Missing arguments!"
  echo "neff_db.sh SEQ_GLOB MODEL OUT-FILE"
  exit 1
fi

GLOB=$1
MODEL=$2
OUT_FILE=$3

rsub --quiet --no-sync --mult 500 -g "$GLOB" -c "neff_seq.sh FILENAME $MODEL >> $OUT_FILE"
