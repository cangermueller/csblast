#!/bin/bash

source $CS/.cs.sh

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
NAME=`basename $WORKDIR`
PLOTDIR=$WORKDIR/plots
mkdir -p $PLOTDIR
BASENAME=$PLOTDIR/$NAME
shopt -s extglob
MODEL=$WORKDIR/cssgd/@(out|*v).crf

if [ -f $MODEL ]; then
  crfvizdist -i $MODEL -o ${BASENAME}_dist -p bw,cw,cws,pc -c 0,6
fi
SGD=$WORKDIR/cssgd/out
if [ -f $SGD ]; then
  sgdviz.pl -i $SGD -o ${BASENAME}_sgd.pdf
fi
CSBLAST=$WORKDIR/csblast
if [ ! -d $CSBLAST ]; then
  CSBLAST=$WORKDIR/csopt/csblast
fi
if [ -d $CSBLAST ]; then
  benchviz.pl \
    -d  $CSBENCH/scop20_1.73_opt/K4000.lib_zX/csblast \
        $CSBENCH/scop20_1.73_opt/K4000.lib_xX/csblast \
        $CSBENCH/scop20_1.73_opt/blast \
        $CSBLAST \
    -l  "K4000.lib z" \
        "K4000.lib x" \
        "blast" \
        "${NAME//_/ }" \
    -o  ${BASENAME}
fi
if [ $CRFVIZ -eq 1 ]; then
  crfviz -i $MODEL -o ${BASENAME}_crf -m 1
  if [ ! -z "$WCENTER" -a ! -z "$WDECAY" ]; then
    crfviz -i $MODEL -o ${BASENAME}_crf_wcenter${WCENTER}_wdecay${WDECAY} -m 1 --weight-center $WCENTER --weight-decay $WDECAY
  fi
fi
