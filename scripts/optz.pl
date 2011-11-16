#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Template;
use My::Utils qw(filename);
use File::Spec::Functions;
use Cwd qw(abs_path);


=pod
=head1 NAME
  
  optz.pl - Optimizes the target diversity z for a given model.

=head1 SYNOPSIS

  optz.pl [OPTIONS] -m MODEL

  OPTIONS:
    -m, --model MODEL   The model used for optimization.
    -o, --out-dir DIR   The output directory.
    -d, --db DB         The name of the database [def: scop20_1.75].
    -t, --tset          Use the test set instead of the optimization set.
    -c, --cat CAT       Name of the category.
    -h, --help          Print this message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $model;
my $outdir;
my $db = "scop20_1.75";
my $db_type;
my $test;
my $cat;

my $csopt_sh = q(#!/bin/bash

source $CS/.cs.sh

DB=$DBS/[% db_type %]
BASEDIR=[% outdir %]

csopt \
  -p $BASEDIR/csopt.yml \
  -b $BASEDIR/csblast.sh \
  -o $BASEDIR/csopt.out \
  -w $BASEDIR/workdir \
  -d $BASEDIR \
  -e ${DB}_db \
  -g "$DB/*.seq" \
  -n 2 \
  -r 1 \
  --mult 100 \
	--seed 0);

my $csopt_yml = q(---
z:  
  order:    1
  value:    12.0
  add:      0.5
  min:      10
  max:      15);

my $csblast_sh = q(#!/bin/bash

source $CS/.cs.sh;

csblast \
  -i FILENAME \
  -d <%= seqfile %> \
  -D [% model %] \
  -z <%= z %> \
  --blast-path $BLAST_PATH \
  -o <%= "#{csblastdir}/#{basename}.bla" %> \
  -e 1e5 \
  -v 10000 \
  -b 0 \
  2> /dev/null);


### Initialization ###


GetOptions(
  "m|model=s"   => \$model,
  "o|out-dir=s" => \$outdir,
  "d|db=s"      => \$db,
  "t|tset!"     => \$test,
  "c|cat=s"     => \$cat,
  "h|help"      => sub { pod2usage(2); }
) or die pod2usage(1);

unless ($model) { die pod2usage("No model provided!"); }
$model = abs_path($model);
my $type = $test ? "test" : "opt";
$db_type = sprintf("%s_%s", $db, $type);
unless ($outdir) {
  $outdir = sprintf("%s/bench/%s/%s/share", $ENV{"CSD"}, $db, $type);
  unless (-d $outdir) { die "'$outdir' does not exist!"; }
  if ($cat) { $outdir .= "/$cat"; }
  $outdir .= sprintf("/%s_optz", filename("$model"));
} else {
  $outdir = abs_path($outdir);
}


### Create script files and submit the optimization job ###


my $tplvars = {
  outdir  => $outdir,
  model   => $model,
  db_type => $db_type
};
my $tpl = Template->new;
mkdir $outdir;
$tpl->process(\$csopt_sh, $tplvars, catfile($outdir, "csopt.sh"));
$tpl->process(\$csopt_yml, $tplvars, catfile($outdir, "csopt.yml"));
$tpl->process(\$csblast_sh, $tplvars, catfile($outdir, "csblast.sh"));
system(sprintf("qsub %s", catfile($outdir, "csopt.sh")));
