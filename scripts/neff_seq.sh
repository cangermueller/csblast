#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -ne 3 ]; then
  echo "neff_seq.sh MODEL ADMIX SEQ"
  exit 1
fi


MODEL=$1
ADMIX=$2
SEQ=$3
CP=`basename $SEQ`
CP=/tmp/neff$$${CP%.*}.prf

csbuild -i $SEQ -D $MODEL -x $ADMIX -o $CP > /dev/null
if [ $? -ne 0 ]; then exit 1; fi
NEFF_LINE=`cscp_neff -i $CP`
if [ $? -ne 0 ]; then exit 1; fi
rm -f $CP
NEFF=`echo $NEFF_LINE | awk '{ if (match($0, /^Neff[[:space:]]*=[[:space:]]*([[:digit:]]+(\.[[:digit:]]+)?)/, a)) print a[1]; else exit 1 }'`
if [ $? -ne 0 ]; then echo $NEFF_LINE; exit 1; fi
echo $NEFF
