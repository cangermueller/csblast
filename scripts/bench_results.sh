#!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -lt 1 ]; then
  echo "bench_results.sh DIR [TYPE]"
  exit 1
fi

DIR=$1
TYPE=${2:-0}
case $TYPE in
  0) DB=$DBS/scop20_1.75_bench;;
  1) DB=$DBS/scop20_1.75_opt;;
  2) DB=$DBS/scop20_1.73_bench;;
  *) DB=$DBS/scop20_1.73_opt;;
esac
DB_FILE=${DB}_db
PLOTS="tpfp,wtpfp,rocx,evalue,pvalue"

function is_result_dir {
  D=$1
  if [ -d $D -a `ls $D/*.bla 2> /dev/null | wc -w` -gt 0 -a `ls $D/*.dat 2> /dev/null | wc -w` -eq 0 ]; then
    echo 1
  else 
    echo 0
  fi
}

function submit_job {
  D=$1
  qsub -b y csbin -i \"$D/\*.bla\" -d $D -s $DB_FILE -p $PLOTS
}

if [ `is_result_dir $DIR` -eq 1 ]; then
  submit_job $DIR
else
  for D in $DIR/*; do
    if [ `is_result_dir $D` -eq 1 ]; then
      submit_job $D
    fi
  done
fi
