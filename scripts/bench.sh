#!/bin/bash

source $HOME/src/cs/.cs.sh

function synopsis {
  cat <<END
bench.sh - Benchmarks a given model

bench.sh [OPTIONS] -m MODEL

  OPTIONS:
    -m MODEL                 The input model
    -x TAU                   Admixture coefficient tau
    -z GAMMA                 Admixture coefficient gamma
    -c CATEGORY              Output catecory
    -d DB                    Database name
    -t                       Use the test set instead of the optimization set
    -h                       Show this text
END
}

DB="scop20_1.75"
MODEL=
X=
T=
CAT=share
TSET=0


### Initialization ###


while getopts "m:x:z:c:d:th" OPT; do
  case $OPT in
    m) MODEL=$OPTARG;;
    x) X=$OPTARG;;    
    z) Z=$OPTARG;;    
    c) CAT=$OPTARG;;    
    d) DB=$OPTARG;;
    t) TSET=1 ;;
    h) synopsis; exit 0 ;;
    *) exit 1
  esac
done
if [ -z "$MODEL" ]; then 
  echo "No model provided!"
  synopsis
  exit 1
fi
if [ -z "$X" -a -z "$Z" ]; then Z=12.5; fi

if [ $TSET -eq 0 ]; then
  TYPE="opt"
else
  TYPE="test"
fi
DB_DIR="$DBS/${DB}_$TYPE"
if [ ! -d $DB_DIR ]; then
  echo "'$DB_DIR' does not exist!"
  exit 1
fi
DB_FILE=${DB_DIR}_db
if [ ! -f $DB_FILE ]; then 
  echo "'$DB_FILE' does not exits!"
  exit 1
fi
if [ -z "$X" ]; then 
  ADMIX="-z $Z"
else
  ADMIX="-x $X"
fi
OUTDIR=$CSD/bench/$DB/$TYPE/share
if [ ! -d "$OUTDIR" ]; then 
  echo "'$OUTDIR' does not exist!"
  exit 1
fi
if [ $CAT ]; then OUTDIR=$OUTDIR/$CAT; fi
OUTDIR=$OUTDIR/`basename $MODEL`


### Run ##


CMD=
if [ $MODEL == "blast" ]; then
  CMD="$BLAST_PATH/blastpgp -i FILENAME -o $OUTDIR/BASENAME.bla -d $DB_FILE -e 1e5 -v 10000 -b 0"
else
  CMD="csblast -D $MODEL --blast-path $BLAST_PATH -i FILENAME -o $OUTDIR/BASENAME.bla -d $DB_FILE -e 1e5 -v 10000 -b 0 $ADMIX"
fi
mkdir -p $OUTDIR
rsub --logfile $HOME/jobs/bench$$.log --mult 100 --quiet -g "$DB_DIR/*.seq" -c "$CMD"
qsub -b y csbin -i \"$OUTDIR/\*.bla\" -d $OUTDIR -s $DB_FILE -p tpfp,wtpfp,rocx,evalue,pvalue
