#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -ne 3 ]; then
  echo "Missing arguments!"
  echo "db_build_a3m.sh INPUT-A3M OUTPUT-DIR NEFF"
  exit 1;
fi

IN_A3M=$1
OUT_DIR=$2
NEFF_MIN=$3

if [ -e $OUT_DIR/`basename $IN_A3M .a3m`_1.a3m ]; then exit 0; fi

NEFF=`hhmake -i $IN_A3M -o /dev/null | awk '{ if (match($0, /exp\(entropy\) = ([^[:space:]]+)/, m)) print m[1]}'`
if [ -z $NEFF ]; then
  echo "Error executing hhmake"
  exit 1
fi
if [[ $NEFF > $NEFF_MIN ]]; then
  /cluster/bioprogs/hhblits/hhblits -i $IN_A3M -o /dev/null -oalis $OUT_DIR/`basename $IN_A3M .a3m` -n 8 -neffmax 11 -cov 80 -mact 0.1 -cpu 4
fi
