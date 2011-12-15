#!/bin/bash

source $HOME/src/cs/.cs.sh

function synopsis {
  cat <<END
bench.sh - Benchmarks a given model

bench.sh [OPTIONS] -m MODEL -o OUTDIR

  OPTIONS:
    -m MODEL        The input model
    -o OUTDIR       Output directory
    -d DB           Database name
    -p PARAMS       Command line parameters
    -s              Do not submit the job
    -h              Show this text
END
}

function check_arg {
  if [ -z "$2" ]; then 
    echo "No $1 provided!"
    synopsis
    exit 1
  fi
}

MODEL=
OUTDIR=
DB="scop20_1.73_opt"
PARAMS=
SUBMIT=1


### Initialization ###


while getopts "m:o:d:p:sh" OPT; do
  case $OPT in
    m) MODEL=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    d) DB=$OPTARG;;
    p) PARAMS=$OPTARG;;
    s) unset SUBMIT;;
    h) synopsis; exit 0 ;;
    *) exit 1
  esac
done
check_arg model $MODEL
check_arg outdir $OUTDIR

DB_DIR="$DBS/$DB"
DB_FILE="${DB_DIR}_db"
if [ ! -d $DB_DIR ]; then
  echo "'$DB_DIR' does not exist!"
  exit 1
fi
if [ ! -f $DB_FILE ]; then 
  echo "'$DB_FILE' does not exits!"
  exit 1
fi


### Create and submit job scripts ###


# command file
mkdir -p $OUTDIR
if [ $MODEL == "blast" ]; then
  CMDFILE="$OUTDIR/blast.sh"
  cat > $CMDFILE <<END 
#!/bin/bash
$BLAST_PATH/blastpgp \\
  -i FILENAME \\
  -o $OUTDIR/BASENAME.bla \\
  -d $DB_FILE \\
  -e 1e5 -v 10000 -b 0 $PARAMS
END
else
  CMDFILE="$OUTDIR/csblast.sh"
  cat > $CMDFILE <<END
#!/bin/bash
csblast \\
  -D $MODEL \\
  -i FILENAME \\
  -o $OUTDIR/BASENAME.bla \\
  -d $DB_FILE \\
  --blast-path $BLAST_PATH \\
  -e 1e5 -v 10000 -b 0 $PARAMS
END
fi

# runner file
RUNFILE="$OUTDIR/RUNME"
cat > $RUNFILE <<END
#!/bin/bash
rsub --logfile $HOME/jobs/bench$$.log --mult 100 --quiet -g "$DB_DIR/*.seq" -s $CMDFILE
SHOULD=\`ls $DB_DIR/*.seq 2> /dev/null | wc -l\`
IS=\`ls $OUTDIR/*bla 2> /dev/null\`
if [ ! \$IS -eq \$SHOULD ]; then
  echo "There are only \$IS instead of \$SHOULD blast result files in '$OUTDIR'!"
  exit 1
fi
csbin -i "$OUTDIR/*.bla" -d $OUTDIR -s $DB_FILE -p tpfp,wtpfp,rocx,evalue,pvalue
if [ ! \$? -eq 0 ]; then
  echo "Error calling csbin!"
  exit 1
fi
find $OUTDIR -maxdepth 1 -type f -name "*.bla" -delete
END
chmod u+x $RUNFILE
if [ $SUBMIT ]; then
  qsub -N csbench -o "$OUTDIR/log" $RUNFILE
fi
