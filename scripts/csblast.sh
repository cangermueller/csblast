#!/bin/bash

MODEL=${1:-$K4000}
SEQ=${2:-$ZINC}

csblast -i $SEQ -D $MODEL -d $DBS/scop20_1.75_opt_db --blast-path /cluster/bioprogs/blast/bin
