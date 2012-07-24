#!/bin/bash

function synopsis {
  cat <<END
csalibench.sh - Performs alignment benchmark for a given model

csalibench.sh [OPTIONS] -m MODEL -o OUTDIR

  OPTIONS:
    -m MODEL        The input model
    -o OUTDIR       Output directory
    -d DB           Database name
    -p PARAMS       Command line parameters
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
DB="scop30_1.75_ali_seqid0.3_tmscore0.6"
PARAMS=
SUBMIT=

CSVERSION=2.2.1
CSBLAST="csblast-$CSVERSION"
E="1e4"
V=100
B=100
MULT=100


### Initialization ###


while getopts "m:o:d:j:p:k:sh" OPT; do
  case $OPT in
    m) MODEL=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    d) DB=$OPTARG;;
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
  DB_DIR=$DBS/$DB
fi
if [ ! -d $DB_DIR ]; then
  echo "'$DB_DIR' does not exist!"
  exit 1
fi


### Create and submit job scripts ###


# command file
mkdir -p $OUTDIR
if [ $MODEL == "blast" ]; then
  CMDFILE="blast.sh"
  cat > $OUTDIR/$CMDFILE <<END 
#!/bin/bash

<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>
<% dbfile = "#{dirname}_db" %>

$BLAST_PATH/blastpgp \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -d <%= dbfile %> \\
  -e $E -v $V -b $B $PARAMS
END

else
  CMDFILE="csblast.sh"
  cat > $OUTDIR/$CMDFILE <<END
#!/bin/bash

<% model = "$MODEL" %>
<% dirbase = "$OUTDIR/BASENAME" %>
<% outfile = "#{dirbase}.bla" %>
<% dbfile = "#{dirname}_db" %>

$CSBLAST \\
  -D <%= model %> \\
  -i FILENAME \\
  -o <%= outfile %> \\
  -d <%= dbfile %> \\
  --blast-path $BLAST_PATH \\
  -e $E -v $V -b $B $PARAMS
END
fi

# runner file
RUNFILE="$OUTDIR/RUNME"
cat > $RUNFILE <<END
#!/bin/bash

DB=$DB_DIR
OUTDIR=$OUTDIR
rsub -g "\$DB/fam/*/*.seq" -s \$OUTDIR/$CMDFILE --logfile \$OUTDIR/log --mult $MULT --quiet --seed $RANDOM 
SHOULD=\`ls \$DB/fam/*/*.seq 2> /dev/null | wc -l\`
IS=\`ls \$OUTDIR/*.bla 2> /dev/null | wc -l\`
if [ ! \$IS -eq \$SHOULD ]; then
  echo "There are only \$IS instead of \$SHOULD blast result files in '\$OUTDIR'!"
  exit 1
fi
echo "Done!"
END

chmod u+x $RUNFILE
if [ $SUBMIT ]; then
  qsub -N csalibench $RUNFILE
fi
