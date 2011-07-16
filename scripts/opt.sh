#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -lt 2 ]; then
  echo "opt.sh OPT MODEL [TYPE] [DIR]"
  exit 1
fi

OPT=$1
MODEL=$2
TYPE=${3:-0}
DIR=${4:-share}
if [ ! -f $YML ]; then "YML does not exist!"; exit 1; fi
if [ ! -f $MODEL ]; then "MODEL does not exist!"; exit 1; fi

DB=
ROOT=
case $TYPE in
  0) DB=$DBS/scop20_1.75_bench; ROOT=$CSBENCH;;
  1) DB=$DBS/scop20_1.75_opt; ROOT=$CSOPT;;
  2) DB=$DBS/scop20_1.73_bench; ROOT=$CSD/bench/scop20_1.73/bench;;
  *) DB=$DBS/scop20_1.73_opt; ROOT=$CSD/bench/scop20_1.73/opt;;
esac
DB_FILE=${DB}_db

YML=$CSE/opt$OPT.yml
BASENAME=`printf '%s/%s/%s' $ROOT $DIR $(basename $MODEL)`
WORKDIR=${BASENAME}_${OPT}opt
OUTFILE=${BASENAME}_${OPT}opt.log
mkdir -p $WORKDIR
CSBLAST=$WORKDIR/csblast.sh
cat > $CSBLAST <<END
csblast -i FILENAME -D $MODEL --blast-path /cluster/bioprogs/blast/bin -o <%= "#{workdir}/#{basename}.bla" %> -d $DB_FILE -e 1e5 -v 10000 -b 0 `opt_args.pl $YML`
END
csopt -p $YML -b $CSBLAST -d $WORKDIR -o $OUTFILE -g "$DB/*.seq" -s $DB_FILE -r 1 -k
