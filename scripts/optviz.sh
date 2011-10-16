#!/bin/bash

WORKDIR=$PWD
CRFVIZ=0
WCENTER=
WDECAY=

function help {
  cat <<END 
optviz.sh [OPTIONS]
  OPTIONS:
    -i WORKDIR         Working directory
    -c                 Visualize crf via crfviz.
    -w WEIGHT-CENTER   Weight-center for crf visualization.
    -d WEIGHT-DECAY    Exponential decay for crf visualization.
    -h                 Print this help message.
END
}

while getopts ":i:cw:d:h" OPT; do
  case $OPT in
    'i') WORKDIR=$OPTARG;;
    'c') CRFVIZ=1;;
    'w') WCENTER=$OPTARG;;
    'd') WDECAY=$OPTARG;;
    'h') 
      help
      exit 0;;
    *)
      echo "Invalid option!" 1>&2
      help 1>&2
      exit 1
  esac
done
if [ ! -z "$WCENTER" -a  ! -z "$WDECAY" ]; then CRFVIZ=1; fi
NAME=`basename $(dirname $WORKDIR)`_`basename $WORKDIR`
PLOTDIR=$WORKDIR/plots
mkdir -p $PLOTDIR
BASENAME=$PLOTDIR/$NAME
MODEL=$WORKDIR/cssgd/*v.crf

if [ -f $MODEL ]; then
  crfvizdist -i $MODEL -o ${BASENAME}_dist -p bw,cw,cws,pc -c 0,6
fi
SGD=$WORKDIR/cssgd/out
if [ -f $SGD ]; then
  sgdviz.pl -i $SGD -o ${BASENAME}_sgd.pdf
fi
CSBLAST=$WORKDIR/csblast
if [ -d $CSBLAST ]; then
  benchviz.pl \
    -d  $REF/nr20_1hhblits_g1.00_n4.0_m4.0_y7.0_N3.0M_K4000_b10.0_c10.0_p10.0_q0_mK4000.lib_v.crf \
        $REF/K4000.lib \
        $REF/nr30_neff2.5_hhblast_1round_W13_N5M_neff6.0_K200.crf \
        $REF/blast \
        $CSBLAST \
    -l  "K4000.crf" \
        "K4000.lib" \
        "andreas" \
        "blast" \
        "${NAME//_/ }" \
    -o  ${BASENAME}

fi
if [ $CRFVIZ -eq 1 ]; then
  for M in 1 2; do
    crfviz -i $MODEL -o ${BASENAME}_crf_m$M -m $M
    if [ ! -z "$WCENTER" -a ! -z "$WDECAY" ]; then
      crfviz -i $MODEL -o ${BASENAME}_crf_m${M}_w${WCENTER}_d${WDECAY} -m $M --weight-center $WCENTER --weight-decay $WDECAY
    fi
  done
fi
