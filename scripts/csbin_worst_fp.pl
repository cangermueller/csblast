#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Spec::Functions qw(catdir);
use File::Basename qw(dirname);


=pod
=head1 NAME
  
  csbin_worst_fp.pl - Investigates the worst false positives produced by csblast

=head1 SYNOPSIS

  csbin_worst_fp.pl [OPTIONS] -i CSBIN-FILE

  OPTIONS:
    -i, --infile FILE     csbin output file with statistics
    -o, --outdir DIR      Output directory [def: .]
    -b, --csblast FILE    Script for repeating the csblast search
    -e, --evalue [0;inf[  Evalue cut-off
    -d, --db STRING       Database name [def: scop20_1.73_opt]
    -v, --verbose INT     Verbosity [def: 1]
    -h, --help            Print this help text

=head1 AUTHOR
  
  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $infile;
my $outdir = "." ;
my $csblast_file;
my $evalue_cut = 1e-3;
my $db = "scop20_1.73_opt";
my $verbose = 1;
my @hits;


### Initialization ###


GetOptions (
  "i|infile=s"  => \$infile,
  "o|outdir=s"  => \$outdir,
  "b|csblast=s" => \$csblast_file,
  "e|e-value=f" => \$evalue_cut,
  "d|db=s"      => \$db,
  "v|verbose=i" => \$verbose,  
  "h|help"      => sub { pod2usage(2); }
) or pod2usage(1);
unless ($infile) { pod2usage("No input file provided!"); }


### Main ###


# Worst FP hits sorted by E-value:
# No    Query     Hit       Queryfam    Hitfam      Hit-No  E-Value  Weight
# 1     d1hvca_   d1vrta2   b.50.1.1    e.8.1.2     3       3e-07    0.14
open FIN, "< $infile" or die "Can't read from '$infile'!";
my $is_worst_fp = 0;
while (<FIN>) {
  if ($is_worst_fp) {
    if (/^\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\d+\s+(\S+)/) {
      my $hit = {
        query    => $1,
        hit      => $2,
        queryfam => $3,
        hitfam   => $4,
        evalue   => $5 * 1.0
      };
      if ($hit->{evalue} > $evalue_cut) { last; }
      push(@hits, $hit);
    } else {
      last;
    }
  } else {
    if (/^Worst FP hits/) {
      <FIN>;
      $is_worst_fp = 1;
    }
  }
}
close FIN;

@hits = sort({ $a->{evalue} <=> $b->{evalue} } @hits);
my $num = 1;
foreach my $hit (@hits) {
  my $name = sprintf("%02d_%.1e_%s_%s", $num++, $hit->{evalue}, $hit->{query}, $hit->{hit});
  print "Creating '$name' ...\n" if $verbose;
  my $dir = catdir($outdir, $name);
  mkdir $dir;
  &copy_seq($hit->{query}, $dir);
  &copy_seq($hit->{hit}, $dir);
  &hhrepid($hit->{query}, $dir);
  &hhrepid($hit->{hit}, $dir);
  if ($csblast_file) { &call_csblast($hit->{query}); }
  &align($hit, $dir);
}

print "Done!\n" if $verbose;


### Functions ###
  

sub execute {
  my ($cmd) = @_;
  print "$cmd\n" if $verbose > 1;
  `$cmd`;
  if ($?) { die "Error executing '$cmd'!"; }
}

sub seq_file {
  my ($id, $ext) = @_;
  unless ($ext) { $ext = "seq"; }
  return catdir($ENV{"DBS"}, $db, "$id.$ext");
}

sub copy_seq {
  my ($id, $dir) = @_;
  my $seq_file = &seq_file($id);
  &execute("cp $seq_file $dir");
}

sub hhrepid {
  my ($id, $dir) = @_;
  my $seq_file = &seq_file($id, "a3m");
  &execute("$ENV{'HHREPID_PATH'}/bin/hhrepid -i $seq_file -d $ENV{'HHREPID_PATH'}/cal.hhm -o $dir/$id.hhrepid.out -pdir $dir > $dir/$id.hhrepid.log 2> /dev/null");
}

sub call_csblast {
  my ($id) = @_;
  my $seq_file = &seq_file($id);
  &execute("rsub -g $seq_file -s $csblast_file --quiet");
}

sub align {
  my ($hit, $dir) = @_;
  my $query_file = &seq_file($hit->{query});
  my $blast_file = catdir(dirname($infile), "$hit->{query}.bla");
  my $ali_file = catdir($dir, "$hit->{query}.fasta");
  &execute("cp $blast_file $dir");
  &execute("alignhits $blast_file $ali_file -Q $query_file -e $hit->{evalue} -fas");
}
