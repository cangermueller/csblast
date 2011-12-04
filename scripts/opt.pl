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
  
  opt.pl - Optimizes the admixture coefficient tau for a given model.

=head1 SYNOPSIS

  opt.pl [OPTIONS] -m MODEL

  OPTIONS:
    -m, --model MODEL   The model used for optimization.
    -o, --out-dir DIR   The output directory.
    -x, --pc-admix      Optimize the admixture rate x instead of z.
    -d, --db DB         The name of the database [def: scop20_1.75_opt].
    -t, --tset          Use the test set instead of the optimization set.
    -c, --cat CAT       Name of the category.
    -s, --seed INT      Seed for csopt.
    -h, --help          Print this message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $model;
my $outdir;
my $param = "z";
my $db = "scop20_1.75_opt";
my $test;
my $cat;
my $seed = int(rand(1e6));

my $csopt_sh = q(#!/bin/bash

source $CS/.cs.sh

DB=$DBS/[% db %]
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
  --seed [% seed %]);

my %csopt_yml = (
  "z" => q(---
z:  
  order:    1
  value:    12.0
  add:      0.5
  min:      10
  max:      15),

  "x" => q(---
x:  
  order:    1
  value:    0.9
  add:      0.05
  min:      0.0
  max:      1.0)
);

my $csblast_sh = q(#!/bin/bash

source $CS/.cs.sh;

csblast \
  -i FILENAME \
  -d <%= seqfile %> \
  -D [% model %] \
  -[% param %] <%= [% param %] %> \
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
  "x|pc-admix"  => sub { $param = "x"; },
  "d|db=s"      => \$db,
  "c|cat=s"     => \$cat,
  "s|seed=i"    => \$seed,
  "h|help"      => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($model) { die pod2usage("No model provided!"); }
$model = abs_path($model);

unless ($outdir) {
  my $db_base = $db;
  my $type = "all";
  if ($db_base =~ /_(opt|test)$/) {
    $type = $1;
    $db_base =~ s/_$1$//;
  }  
  $outdir = sprintf("%s/bench/%s/%s/share", $ENV{"CSD"}, $db_base, $type);
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
  param   => $param,
  db      => $db,
  seed    => $seed
};
my $tpl = Template->new;
mkdir $outdir;
$tpl->process(\$csopt_yml{$param}, $tplvars, catfile($outdir, "csopt.yml"));
$tpl->process(\$csblast_sh, $tplvars, catfile($outdir, "csblast.sh"));
my $csopt = catfile($outdir, "csopt.sh");
$tpl->process(\$csopt_sh, $tplvars, $csopt);
system("chmod u+x $csopt");
system("qsub $csopt");
