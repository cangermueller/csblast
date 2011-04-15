#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 8
#$ -q mpi

source $HOME/src/cs/.cs.sh

DB="$DBS/nr20_sampled_clusters_hhblits_3rounds"
NEFF=8
WLEN=13
N=2000000
S=4

HOME=/cluster/user/christof
for DB_ in $DB; do
	for NEFF_ in $NEFF; do
		for WLEN_ in $WLEN; do
            for S_ in $S; do
                OUT_BASE=`printf '%s/%s_n%i_W%i_N%iM_s%i' $CST $(basename $DB_) $NEFF_ $WLEN ${N:0:2} $S`
                EXT="tpr"
                if [ $S_ == 1 ]; then EXT="tsq"; fi
                cstrainset -d $DB_ -o $OUT_BASE.$EXT -s $S_ --neff $NEFF_ --wlen $WLEN_ --size $N &> $OUT_BASE.log
            done
		done
	done
done



