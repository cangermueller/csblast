#!/bin/bash

source $HOME/src/cs/.cs.sh
source $HOME/lib/utils.sh

if [ $# -eq 0 ]; then
	echo "neffviz.sh MODEL [TRAINSET]"
	exit 1
fi

MODEL=$1
TRAINSET=$2
OUT_DIR=$CSC/plots/neff/`basename $MODEL`
ADMIX=`seq 0.85 0.05 1.0`
QJOBS=

mkdir -p $OUT_DIR
for A in $ADMIX; do
  qsubmit $CSS/neff_model.sh $MODEL $A
done
if [ ! -z $TRAINSET ]; then
	qsubmit -b y cstrainset_neff -i $TRAINSET -x $OUT_DIR/trainset_x.dat -y $OUT_DIR/trainset_y.dat
fi
qwait
R --vanilla --slave -f $CSS/neff_plot.R --args $OUT_DIR/neff.pdf `basename $MODEL` $OUT_DIR/*.dat
