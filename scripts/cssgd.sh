#!/bin/sh

#$ -cwd
#$ -pe pe_mpi 8
#$ -q mpi

source $HOME/src/cs/.cs.sh

if [ $# -eq 0 ]; then 
    echo "Missing arguments!"
    echo "cssgd.sh [BIN] TRAINSET"
    exit 1
fi

BIN=cssgd
if [ $1 == "DEBUG" ]; then 
    BIN=${BIN}_debug 
    shift
fi
TRAINSET=$1
VALSET=$TRAINSET
K=75       # states
P=2        # prior
SC=1.0     # sigma context
SD=1.0     # sigma decay
SB=1.0     # sigma bias 

BASENAME=`basename $TRAINSET`
BASENAME=${BASENAME%\.*}`printf '_K%d_p%d_sb%1.2f_sc%1.2f_sd%1.2f' $K $P $SB $SC $SD`
if [ $2 ]; then BASENAME=${BASENAME}_$2; fi
$BIN	  -i $TRAINSET \
		    -j $VALSET \
		    -o $CM/$BASENAME.crf \
        -e 0.0001 \
        -B 500 \
        -K $K \
        -P $P \
        --sigma-bias $SB \
        --sigma-context $SC \
        --sigma-decay $SD \
		    &> $CM/$BASENAME.log
