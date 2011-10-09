#!/bin/bash

WORKDIR=${1:-$PWD}
NAME=`basename $(dirname $WORKDIR)`_`basename $WORKDIR`
PLOTDIR=$WORKDIR/plots
mkdir -p $PLOTDIR
BASENAME=$PLOTDIR/$NAME

MODEL=$WORKDIR/cssgd/*v.crf
if [ -f $MODEL ]; then
  crfvizdist -i $MODEL -o ${BASENAME}_dens -p bw,cw,cws,pc -c 0,6
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
