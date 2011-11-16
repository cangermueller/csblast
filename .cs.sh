export CS=${CS:-$HOME/src/cs}

export CSB=$CS/bin
export CSS=$CS/scripts
export CS_DATA=$CS/data

export CSD=$HOME/data/cs
export CSM=$CSD/models
export CSC=$CSM/crf
export CST=$CSD/trainsets

export CSBENCH=$CSD/bench/scop20_1.75
export CSOPT=$CSBENCH/opt
export CSTEST=$CSBENCH/test
export CSCROSS=$CSBENCH/cross

export CSE=$HOME/etc/cs
export CSOZ=$CSE/optz.yml
export CSOX=$CSE/optx.yml

export K4000=$CSM/lib/K4000.lib

export DBS=$HOME/databases
export SEQS=$HOME/seqs
export ZINC=$SEQS/zinc_finger.seq

export TOPT=$CST/nr20_1hhblits_t_g1.00_n4.0_m4.0_y7.0_N3.0M.tsq
export K4000CRF=$CSC/00_opt/nr20_1hhblits_g1.00_n4.0_m4.0_y7.0_N3.0M_K4000_b10.0_c10.0_p10.0_q0_mK4000.lib_v.crf

export BLAST_PATH=/cluster/bioprogs/blast/bin
export PATH=$PATH:$CSB:$CSS

export OS=$CSOPT/cssgd
export OST=$OS/test
export OSS=$OS/N300000
export OSL=$OS/N3.0M
export OSM=$OS/N1.0M

export OT=$CSOPT/cstrainset
export OTT=$OT/test
export OTN=$OT/nr20_1hhblits
export OTU=$OT/uni20v2

export OB=$CSOPT/csblast
export OBU=$OB/uni20v2_t_N3.0M_g1.0_u4.00_U10.00_M

export REF=$CSOPT/share/ref


function csseriesz {
  MODEL=$1  
  TYPE=${2-0}
  for Z in `seq 9.0 0.5 14`; do
    qsub $CSS/bench.sh $Z $MODEL $TYPE series
  done
}

function csseriesx {
  MODEL=$1  
  TYPE=${2-0}
  for X in `seq 0.8 0.02 1.0`; do
    qsub $CSS/bench.sh $X $MODEL $TYPE series
  done
}

function cpbench {
  TO=${1%/}
  shift
  FROM=$@
  for F in $FROM; do
    if [ -d $F ]; then
      T=$TO/`basename $F`
      mkdir -p $T
      cp $F/*dat $T
    fi
  done
}

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
alias sgdresults="sgdresults.pl -i *log > results"
