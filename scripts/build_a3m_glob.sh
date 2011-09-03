#!/bin/bash

source $HOME/src/cs/.cs.sh


NEFF_MIN=2.5
OUT_DIR=$DBS/uni20v2
if [ $1 ]; then 
  DB=$1
else 
  DB=/data/fs03/scratch_ssd_raid/hhblits_dbs/uniprot20_29Mar11/second_round/alignments
fi

function build_a3m {
    SEQ=$1
    OUT_BASE=$OUT_DIR/`basename $SEQ .a3m`
    if [ ! -e "${OUT_BASE}_1.a3m" ]; then
        NEFF=`hhmake -i $SEQ -o /dev/null | awk '{ if (match($0, /exp\(entropy\) = ([^[:space:]]+)/, m)) print m[1]}'`
        if [ -z "$NEFF" ]; then
          echo "Error executing hhmake"
          exit 1
        fi
        if [[ $NEFF > $NEFF_MIN ]]; then
          rsub --quiet --no-sync -c "/cluster/bioprogs/hhblits/hhblits -i $SEQ -o /dev/null -oalis $OUT_BASE -n 8 -neffmax 12 -mact 0.35 -cpu 1"
          cp $SEQ $OUT_BASE.seq  
        fi
    fi
}

mkdir -p $OUT_DIR
find $DB -name "*a3m" | while read SEQ; do build_a3m $SEQ; done;
