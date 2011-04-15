#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -ne 2 ]; then
  echo "Missing arguments!"
  echo "neff_seq.sh SEQ_FILE MODEL"
  exit 1
fi

SEQ_FILE=$1
MODEL=$2
CP=`basename $SEQ_FILE`
CP=${CP%.*}.prf

csbuild -i $SEQ_FILE -D $MODEL -o $CP > /dev/null
if [ $? -ne 0 ]; then exit 1; fi
NEFF=`count_profile_neff -i $CP`
if [ $? -ne 0 ]; then exit 1; fi
echo $NEFF | awk '{ match($0, /^Neff[[:space:]]*=[[:space:]]*([^[:space:]]+)/, a); print a[1] }'
rm -f $CP
