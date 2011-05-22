#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 1
#$ -q mpi

source $HOME/src/cs/.cs.sh

if [ $# -lt 1 ]; then
	echo Model missing!
	echo Usage: profile_logo.sh MODEL_FILE+
	exit 1
fi
MODEL_FILES=$@

ROOT=$CSD/profile_logos
for MODEL_FILE in $MODEL_FILES; do
	MODEL_NAME=`basename $MODEL_FILE`
	mkdir -p $ROOT/$MODEL_NAME

	for SEQ in $ROOT/seqs/*seq; do
		OUT_BASE=`basename $SEQ`
		OUT_BASE=$ROOT/$MODEL_NAME/${MODEL_NAME}_${OUT_BASE%\.???}
		csbuild -i $SEQ -o $OUT_BASE.prf -D $MODEL_FILE
		csviz -i $OUT_BASE.prf -o $OUT_BASE.pdf
	done
done


