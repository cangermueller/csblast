#!/bin/bash


source $HOME/src/cs/.cs.sh

if [ $# -lt 2 ]; then
	echo "bench.sh NEFF|ADMIX MODEL|BENCH_SCRIPT [TYPE] [DIR]" 1>&2
	exit 1
fi

NEFF=$1
if [[ ! $NEFF =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
  echo "NEFF|ADMIX invalid!"
  exit 1
fi
MODEL=$2
if [ ! -e $MODEL ]; then
  echo "MODEL|BENCH_SCRIPT does not exist!"
  exit 1
fi
TYPE=${3:-0}
DIR=$4
if [ -z $DIR ]; then
  if [ ${MODEL:0:${#CSC}} == $CSC ]; then
    DIR=`dirname $MODEL`
    DIR=${DIR#$CSC/}
  else
    DIR=share
  fi
fi


case $TYPE in
  0) DB=$DBS/scop20_1.75_bench;;
  1) DB=$DBS/scop20_1.75_opt;;
  2) DB=$DBS/scop20_1.73_bench;;
  *) DB=$DBS/scop20_1.73_opt;;
esac
DB_FILE=${DB}_db
PLOTS="tpfp,wtpfp,rocx,evalue,pvalue"

echo test
if [ ${MODEL##*.} == "sh" ]; then
  rsub --mult 100 --quiet -g "$DB/*.seq" -s $MODEL
else
  if [ ${NEFF%.*} -gt 1 ]; then
    PC="z"
  else 
    PC="x"
  fi
  NAME=`printf '%s_%s%.2f' $(basename $MODEL) $PC $NEFF`
  case $TYPE in
    0) OUTDIR="$CSD/bench/scop20_1.75/bench";;
    1) OUTDIR="$CSD/bench/scop20_1.75/opt";;
    2) OUTDIR="$CSD/bench/scop20_1.73/bench";;
    *) OUTDIR="$CSD/bench/scop20_1.73/opt";;
  esac
  OUTDIR=$OUTDIR/$DIR/$NAME
  if [ -e $OUTDIR ]; then exit 0; fi
  mkdir -p $OUTDIR
  rsub --logfile $HOME/jobs/bench$$.log --mult 100 --quiet -g "$DB/*.seq" -c "csblast -i FILENAME -D $MODEL --blast-path /cluster/bioprogs/blast/bin -o $OUTDIR/BASENAME.bla -d $DB_FILE -e 1e5 -v 10000 -b 0 -$PC $NEFF"
  qsub -b y csbin -i \"$OUTDIR/\*.bla\" -d $OUTDIR -s $DB_FILE -p $PLOTS
fi
