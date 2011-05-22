#!/bin/bash

#$ -q normal
#$ -pe threads.pe 2

source $HOME/src/cs/.cs.sh

DB=$DBS/nr20_neff3.0_1hhblits
#DB=$DBS/nr20_3hhblits
#DB=$DBS/scop20_1.75_opt
#DB=$DBS/uniprot20_neff2.5_neff2.5-6.5
#DB_PC=$DBS/uniprot20_neff2.5_neff7.0-12.0
#BASENAME="uniprot20_neff2.5-6.0_7.0-12.0"

S=3
N=3000000
G=0.5
NEFF_X=6.0
NEFF_Y=6.0
NEFF_Z=0.0
WLEN=13

if [ ! -z $DB_PC ]; then E="-e $DB_PC"; fi
if [ $N -lt 1000000 ]; then
    N_SHORT=$N
else 
    N_SHORT=`printf '%d.%dM' ${N:0:1} ${N:1:1}`
fi
for DB_ in $DB; do
  BASENAME=${BASENAME:-`basename $DB_`}
  for S_ in $S; do
    if [ $S_ -eq 1 ]; then
      OUT_BASE=`printf '%s/%s_s%i_n%.1f_z%.1f_W%i_N%s' $CST $BASENAME $S_ $NEFF_X $NEFF_Z $WLEN $N_SHORT`
      EXT="tsq"
    else 
      OUT_BASE=`printf '%s/%s_s%i_g%.2f_n%.1f_m%.1f_z%.1f_W%i_N%s' $CST $BASENAME $S_ $G $NEFF_X $NEFF_Y $NEFF_Z $WLEN $N_SHORT`
      EXT="tpr"
    fi
    cstrainset -d $DB_ $E -o $OUT_BASE.$EXT -g $G -n $NEFF_X -m $NEFF_Y -z $NEFF_Z -W $WLEN -N $N -s $S_ -x 0.0 -D $K4000 &> $OUT_BASE.log
  done
done



