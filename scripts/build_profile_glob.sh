#!/bin/bash

#$ -cwd
#$ -pe pe_mpi 1
#$ -q mpi

source $HOME/src/cs/.cs.sh

MIN_NEFF=1.0
MAX_NEFF=7.0
DB_ALI=$DBS/uniprot20_neff2.5
DB_PRF=${DB_ALI}_neff$MIN_NEFF-$MAX_NEFF


mkdir -p $DB_PRF
rsub --quiet --no-sync --mult 10 -g "$DB_ALI/*.seq" -c "build_profile.pl -i DIRBASENAME*.a3m -o $DB_PRF/BASENAME.prf -min $MIN_NEFF -max $MAX_NEFF"

