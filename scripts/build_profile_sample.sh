#!/bin/bash

source $HOME/src/cs/.cs.sh

MIN_NEFF=2.5
MAX_NEFF=6.5
MODE_NEFF=4.5
DB_ALI=$DBS/uni20_t
DB_PRF=${DB_ALI}_neff$MIN_NEFF-$MAX_NEFF


mkdir -p $DB_PRF
rsub --quiet --no-sync --mult 10 -g "$DB_ALI/*.seq" -c "build_profile.pl -i DIRBASENAME*.a3m -o $DB_PRF/BASENAME -min $MIN_NEFF -max $MAX_NEFF -mode $MODE_NEFF -filter"

: <<'END'
for F in $DB_ALI/*.seq; do
  echo test
  NAME=`basename $F`
  NAME=${NAME%.*}
  if [ ! -e "$DB_PRF/$NAME.prf" ]; then
    rsub --quiet --no-sync -c "build_profile.pl -i ${F%.*}*.a3m -o $DB_PRF/$NAME -min $MIN_NEFF -max $MAX_NEFF"
  fi
done
END

