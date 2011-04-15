#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 4
#$ -q mpi

source $HOME/src/cs/.cs.sh

if [ $# -eq 0 ]; then
    echo "Missing arguments!"
	echo "neff_model.sh MODEL TRAINSET"
	exit 1
fi

MODEL=$1
TRAINSET=$2
DB="$DBS/scop20_1.75_bench"
ROOT="$CSD/neff"
OUT_DIR=$ROOT/`basename $MODEL`

mkdir -p $OUT_DIR
for SEQ in $DB/*seq; do
	neff_seq.sh $SEQ $MODEL >> $OUT_DIR/prediction.dat
done

if [ ${#TRAINSET} -gt 0 ]; then
	trainset_neff -i $TRAINSET -x $OUT_DIR/trainset_x.dat -y $OUT_DIR/trainset_y.dat
fi

R --vanilla --slave -f $CSS/neff_plot.R --args $OUT_DIR/neff.pdf `basename $MODEL` $OUT_DIR/*.dat

