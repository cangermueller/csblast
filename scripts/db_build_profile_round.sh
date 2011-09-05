#!/bin/bash

source $HOME/src/cs/.cs.sh

ROUND=1
DB_ALI=$DBS/uni20_v
DB_PRF=${DB_ALI}_R$ROUND

mkdir -p $DB_PRF
rsub --quiet --no-sync --mult 100 -g "$DB_ALI/*.seq" -c "if [ -f DIRBASENAME_$ROUND.a3m ]; then csbuild -x 0 -i DIRBASENAME_$ROUND.a3m -o $DB_PRF/BASENAME.prf; cp DIRBASENAME_$ROUND.a3m $DB_PRF/BASENAME.a3m; fi"

