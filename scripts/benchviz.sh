#!/bin/bash

function help {
  cat <<END
benchviz.sh -d "DIRECTORY+" -l "LABEL+"
END
}

DIRS=
LABELS=
while getopts ":d:l:" OPT; do
  case $OPT in
    d) 
      DIRS="$DIRS $OPTARG"
      shift $((OPTIND-1))
      OPTIND=1
      ;;
    l) 
      LABELS="$LABELS $OPTARG"
      shift $((OPTIND-1))
      OPTIND=1
      ;;
  esac
done
OPTS=$@

benchviz.pl \
  -d  $REF/nr20_1hhblits_g1.00_n4.0_m4.0_y7.0_N3.0M_K4000_b10.0_c10.0_p10.0_q0_mK4000.lib_v.crf \
      $REF/K4000.lib \
      $REF/nr30_neff2.5_hhblast_1round_W13_N5M_neff6.0_K200.crf \
      $REF/blast \
      $DIRS \
  -l  "K4000.crf" \
      "K4000.lib" \
      "andreas" \
      "blast" \
      $LABELS \
  $OPTS
