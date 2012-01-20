shopt -s extglob

export CS=${CS:-$HOME/src/cs}

export CSVERSION="2.2.0"

export CSB=$CS/bin
export CSS=$CS/scripts
export CS_DATA=$CS/data

export CSD=$HOME/data/cs
export CSM=$CSD/models
export CSC=$CSM/crf
export CST=$CSD/trainsets

export CSBENCH=$CSD/bench
export CSIBENCH=$CSD/bench_csiblast

export K4000=$CSM/lib/K4000.lib
export K4000CRF=$CSC/uni20v2_t_N6.0M_g1.0_u4.00_U10.00_M_K4000_e0.10_D2.5_P2_b1.0_c1.6_d0.85_p1.0_me0.13.crf
export K4000CRFL=`basename $K4000CRF`

export DBS=$HOME/databases
export SEQS=$HOME/seqs
export ZINC=$SEQS/zinc_finger.seq

export BLAST_PATH=/cluster/bioprogs/blast-2.2.19/bin
if [ ! -d $BLAST_PATH ]; then
  BLAST_PATH=$HOME/bin/blast/bin
fi
export PATH=$PATH:$CSB:$CSS

# Temporary variables

export O5=$CSBENCH/scop20_1.75_opt
export O5M=$O5/models
export O5Mu=$O5M/uni20v2_t_N3.0M_g1.0_u4.00_U10.00_M
export O5MU=$O5M/uni20v2_t_N6.0M_g1.0_u4.00_U10.00_M
export T5=$CSBENCH/scop20_1.75_test
export T5M=$T5/models

export O=$CSBENCH/scop20_1.73_opt
export Ou=$O/uni20v2_t_N3.0M_g1.0_u4.00_U10.00_M
export OU=$O/uni20v2_t_N6.0M_g1.0_u4.00_U10.00_M
export T=$CSBENCH/scop20_1.73_test

export OI=$CSIBENCH/scop20_1.73_opt
export TI=$CSIBENCH/scop20_1.73_test
export OIP=$CSIBENCH/scop20_1.73_opt_1psi_neff1+


### Functions ###


function sgdviz {
  LOGS=$@
  mkdir -p $PWD/plots; 
  for LOG in $LOGS; do
    OUT=`basename $LOG`;
    OUT=$PWD/plots/${OUT%log}pdf
    sgdviz.pl -i $LOG -o $OUT -p ll-val ll-train neff prior
  done
}

function cpop {
  cp ${1%/}/RUNME ${1%/}/*@(sh|yml) $2
}



alias csupdate="source $CS/.cs.sh"
alias visgd='vi -p `find -path "*/0?/cssgd/out"`'
alias viopt='vi -p `find -path "*/0?/csopt/out"`'
