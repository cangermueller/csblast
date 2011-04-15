#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 1
#$ -q mpi

source $HOME/src/cs/.cs.sh

DB=/data/fs03/scratch_ssd_raid/hhblits_dbs/uniprot20_16Nov10/second_round/alignments
NEFF_MIN=2.5
OUT_DIR=$DBS/uniprot20_neff$NEFF_MIN

function build_a3m {
    SEQ=$1
    OUT_BASE=$OUT_DIR/`basename $SEQ .a3m`
    if [ ! -e "${OUT_BASE}_1.a3m" ]; then
        NEFF=`hhmake -i $SEQ -o /dev/null | awk '{ if (match($0, /exp\(entropy\) = ([^[:space:]]+)/, m)) print m[1]}'`
        if [ -z $NEFF ]; then
          echo "Error executing hhmake"
          exit 1
        fi
        if [[ $NEFF > $NEFF_MIN ]]; then
          cp $SEQ $OUT_BASE.seq  
          rsub --quiet --no-sync -c "/cluster/bioprogs/hhblits/hhblits -i $SEQ -o /dev/null -oalis $OUT_BASE -n 8 -neffmax 11 -cov 80 -mact 0.1 -cpu 4"
        fi
    fi
}

mkdir -p $OUT_DIR
find $DB -name "*a3m" | while read SEQ; do build_a3m $SEQ; done;
