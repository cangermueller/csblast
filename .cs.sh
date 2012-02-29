shopt -s extglob

export CS=${CS:-$HOME/src/cs}

export CSVERSION="2.2.1"

export CSB=$CS/bin
export CSS=$CS/scripts


export CSD=$HOME/data/cs
export CSM=$CSD/models
export CSC=$CSM/crf
export CST=$CSD/trainsets

export CSBENCH=$CSD/bench
export CSIBENCH=$CSD/bench_csiblast
export ALIBENCH=$CSD/bench_ali
export CSBENCHO=$CSD/bench.2012-02-04
export CSIBENCHO=$CSD/bench_csiblast.2012-02-04

export L4000=$CSM/lib/K4000.lib
export C4000=$CSC/uni20v2_t_N6.0M_g1.0_u4.00_U10.00_M_K4000_e0.10_D2.5_P2_b1.0_c1.6_d0.85_p1.0_me0.13.crf

export BLAST_PATH=/cluster/bioprogs/blast-2.2.19/bin
export HHREPID_PATH=/cluster/bioprogs/hhrepid

if [ ! $DEFCSPATH ]; then
  export PATH=$PATH:$CSB:$CSS
  export DEFCSPATH=1
fi


### Temporary variables ###


export O=$CSBENCH/scop20_1.73_opt
export T=$CSBENCH/scop20_1.73_test
export OI=$CSIBENCH/scop20_1.73_opt
export TI=$CSIBENCH/scop20_1.73_test

export OO=$CSBENCHO/scop20_1.73_opt
export TO=$CSBENCHO/scop20_1.73_test
export OIO=$CSIBENCHO/scop20_1.73_opt
export TIO=$CSIBENCHO/scop20_1.73_test

export ALI=$ALIBENCH/scop30_1.75_ali_seqid0.3_tmscore0.6

export DO=$DBS/scop20_1.73_opt
export DOS=$DBS/scop20_1.73_opt_db
export DT=$DBS/scop20_1.73_test
export DTS=$DBS/scop20_1.73_test_db

export DA=$DBS/scop30_1.75_ali_seqid0.3_tmscore0.6


### Functions ###


alias csupdate="source $CS/.cs.sh"

function benchdir {
  MODEL=$1
  J=${2:-1}
  shift; shift
  APP=$@

  DIR=v$CSVERSION
  if [ $J -gt 1 ]; then
    DIR=`printf "%s_j%d_h%s" $DIR $J "1e-3"`
  fi
  DIR=${DIR}_`basename $MODEL`
  if [ $APP ]; then
    DIR=${DIR}_$APP
  fi
  echo $DIR
}

function alibenchdir {
  MODEL=`basename $1`
  shift
  APP=$@
  DIR=`printf v%s_%s $CSVERSION $MODEL`
  if [ $APP ]; then
    DIR=${DIR}_$APP
  fi
  echo $DIR
}

function count {
  EXT=$1
  shift
  DIRS=$@
  for DIR in $DIRS; do
    N=`ls ${DIR%\/}/*$EXT 2> /dev/null | wc -l`
    echo $N
  done
}

function rmfast {
  EXT=$1
  shift 
  DIRS=$@
  if [ -z $EXT ]; then
    echo "No file extenstion provided!"
    return 1
  fi
  for DIR in $DIRS; do
    if [ ! -d $DIR ]; then
      echo "'$DIR' does not exist!"
      return 1
    fi
    find $DIR -maxdepth 1 -name "*.$EXT" -delete
  done
}