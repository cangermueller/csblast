#!/bin/bash

ARGS=$@
if [[ ! $ARGS =~ -i ]]; then
  ARGS="$ARGS -i $ZINC"
fi
if [[ ! $ARGS =~ -d ]]; then
  ARGS="$ARGS -d $DBS/scop20_1.73_opt_db"
fi
if [[ ! $ARGS =~ -D ]]; then
  ARGS="$ARGS -D $K4000"
fi
if [[ ! $ARGS =~ --blast-path ]]; then
  ARGS="$ARGS --blast-path $BLAST_PATH"
fi

csblast $ARGS
