#!/bin/bash

if [ $# -eq 0 ]; then
	echo "sample_db.sh DB+"
	exit 1
fi

DBS=$@
for DB in $DBS; do
  echo "Sampling from $DB"
  sample_files.pl -i $DB -o ${DB}_v -p ${DB}_t -n 0.2
done
