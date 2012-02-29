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
    -j ITER         Number of iterations
    -p PARAMS       Command line parameters
    -k KEEP         Keep blast files
    -s              Submit the job
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
J=1
PARAMS=
KEEP=1
SUBMIT=

CSBLAST="csblast-$CSVERSION"
NR="$DBS/nr_2011-12-09/nr"
E="1e5"
V=1000
B=0
H="1e-3"
MULT=(50 10 5)


### Initialization ###


while getopts "m:o:d:j:p:k:sh" OPT; do
  case $OPT in
    m) MODEL=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    d) DB=$OPTARG;;
    j) J=$OPTARG;;
    p) PARAMS=$OPTARG;;
    s) SUBMIT=1;;
    k) KEEP=$OPTARG;;
    h) synopsis; exit 0 ;;
    *) exit 1
  esac
done
check_arg model $MODEL
check_arg outdir $OUTDIR

eval $CSBLAST &> /dev/null
if [ $? -ne 0 ]; then
  echo "'$CSBLAST' binary not found!"
  exit 1
fi
DB_DIR=$DB
if [ ! -d $DB_DIR ]; then
  DB_DIR="$DBS/$DB"
fi
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
  if [ $J -eq 1 ]; then
    # blastpgp 1 round
    cat > $CMDFILE <<END 
#!/bin/bash

<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>

$BLAST_PATH/blastpgp \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -d $DB_FILE \\
  -e $E -v $V -b $B $PARAMS
END
  else
    # blastpgp >1 round
    cat > $CMDFILE <<END 
#!/bin/bash

<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>
<% chkfile = "#{dirbase}.chk" %>

$BLAST_PATH/blastpgp \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -C <%= chkfile %> \\
  -d $NR \\
  -e $E -h $H -j $J $PARAMS

$BLAST_PATH/blastpgp \\
  -i FILENAME \\
  -R <%= chkfile %> \\
  -o <%= outfile %> \\
  -d $DB_FILE \\
  -e $E -v $V -b $B $PARAMS
END
  fi
else
  CMDFILE="$OUTDIR/csblast.sh"
  if [ $J -eq 1 ]; then
    # csblast 1 round
  cat > $CMDFILE <<END
#!/bin/bash

<% model = "$MODEL" %>
<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>

$CSBLAST \\
  -D <%= model %> \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -d $DB_FILE \\
  --blast-path $BLAST_PATH \\
  -e $E -v $V -b $B $PARAMS
END
  else
    # csblast >1 round
    cat > $CMDFILE <<END 
#!/bin/bash

<% model = "$MODEL" %>
<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>
<% chkfile = "#{dirbase}.chk" %>

$CSBLAST \\
  -D <%= model %> \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -C <%= chkfile %> \\
  -d $NR \\
  --blast-path $BLAST_PATH \\
  -e $E -h $H -j $J $PARAMS

$CSBLAST \\
  -D <%= model %> \\
  -i FILENAME \\
  -R <%= chkfile %> \\
  -o <%= outfile %> \\
  -d $DB_FILE \\
  --blast-path $BLAST_PATH \\
  -e $E -v $V -b $B $PARAMS
END
  fi
fi

# runner file
RUNFILE="$OUTDIR/RUNME"
cat > $RUNFILE <<END
#!/bin/bash

OUTDIR=$OUTDIR
rsub --logfile \$OUTDIR/log --mult ${MULT[$((J-1))]} --quiet -g "$DB_DIR/*.seq" -s $CMDFILE
SHOULD=\`ls $DB_DIR/*.seq 2> /dev/null | wc -l\`
IS=\`ls \$OUTDIR/*.bla 2> /dev/null | wc -l\`
if [ ! \$IS -eq \$SHOULD ]; then
  echo "There are only \$IS instead of \$SHOULD blast result files in '\$OUTDIR'!"
  exit 1
fi

csbin -i "\$OUTDIR/*.bla" -d \$OUTDIR -o \$OUTDIR/csbin.dat -s $DB_FILE -p tpfp,wtpfp,fdr,rocx,evalue,pvalue --max-evalue 300
if [ ! \$? -eq 0 ]; then
  echo "Error calling csbin!"
  exit 1
fi
END
if [ $KEEP -eq 0 ]; then
  cat >> $RUNFILE <<END

find \$OUTDIR -maxdepth 1 -type f -name "*.bla" -delete
END
  if [ $J -gt 1 ]; then
    cat >> $RUNFILE <<END
find \$OUTDIR -maxdepth 1 -type f -name "*.chk" -delete
END
  fi
fi
chmod u+x $RUNFILE
if [ $SUBMIT ]; then
  qsub -N csbench -o "$OUTDIR/log" $RUNFILE
fi
