#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 1
#$ -q mpi

source $HOME/src/cs/.cs.sh

DB=$DBS/uniprot20_neff2.5
DB_FROM=/data/fs03/scratch_ssd_raid/hhblits_dbs/uniprot20_16Nov10/second_round/alignments


PREV=1
for F in $DB/*a3m; do
    F=`basename $F`
    F=${F%_?\.*}
    if [ $F != $PREV -a ! -e $DB/$F.seq ]; then
        cp $DB_FROM/${F:0:2}/$F.a3m $DB/$F.seq
        PREV=$F
    fi
done

