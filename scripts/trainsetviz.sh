  #!/bin/bash

source $HOME/src/cs/.cs.sh

if [ $# -eq 0 ]; then
  echo "trainsetviz.sh TRAINSET+"
  exit 1
fi

N=100000
TSETS=$@
for TSET in $TSETS; do
  OUTDIR=`dirname $TSET`/plots
  mkdir -p $OUTDIR
  trainsetviz.pl -i $TSET -o $OUTDIR/`basename $TSET`.pdf -N $N
done
exit 0



# neff_plot.R
for TSET in $TSETS; do
  DIR=`mktemp -d`
  NAME=`basename $TSET`
  echo "$NAME..."
  cstrainset_neff -i $TSET -x $DIR/x.dat -y $DIR/y.dat -N $N > $DIR/log
  if [ $? -ne 0 ]; then exit 1; fi
  R --vanilla --slave -f $CSS/neff_plot.R --args $DIR/$NAME.pdf $NAME $DIR/x.dat $DIR/y.dat
  if [ $? -ne 0 ]; then exit 1; fi
  OUTDIR=`dirname $TSET`/plots
  mkdir -p $OUTDIR
  mv $DIR/$NAME.pdf $OUTDIR
  rm -r $DIR
  echo "$NAME.pdf created!"
done

