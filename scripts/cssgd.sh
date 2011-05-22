#!/bin/sh

#$ -q normal
#$ -pe threads.pe 1 

source $HOME/src/cs/.cs.sh

if [ $# -eq 0 ]; then 
    echo "Missing arguments!"
    echo "cssgd.sh TRAINSET [STATES] [SB] [SC] [SP] [SPE] [MODEL]"
    exit 1
fi

BIN=cssgd
TRAINSET=$1
VALSET=$VSET
K=${2:-200}       # states
P=2
SB=${3:-10.0}     # sigma bias 
SC=${4:-10.0}     # sigma context
SD=1.0            # sigma decay
SP=${5:-10.0}
SPE=${6:-0}
SPD=0.001
MODEL=$7
ETA=0.001
BLOCKS=1000

BASENAME=`basename $TRAINSET`
BASENAME=${BASENAME%\.*}`printf '_K%d_b%.1f_c%.1f_p%.1f' $K $SB $SC $SP`
if [ $SPE -gt 0 ]; then
  SPD=2.0
  BASENAME=`printf '%s_q%d_r%.1f' $BASENAME $SPE $SPD`
fi


if [ ! -z $MODEL ]; then
  BASENAME=${BASENAME}_m`basename $MODEL`
  M="-m $MODEL"
fi

$BIN  -i $TRAINSET \
      -j $VALSET \
		  -o $CSC/$BASENAME.crf \
      -e $ETA \
      -B $BLOCKS \
      -K $K \
      -P $P \
      -b $SB \
      -c $SC \
      -d $SD \
      -p $SP \
      -q $SPE \
      -Q $SPD \
      $M \
      --save \
		  &> $CSC/$BASENAME.log
