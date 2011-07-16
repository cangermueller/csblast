#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Temp qw(tempfile);
use My::Utils qw(filename);


=pod
=head1 NAME

  neffseq.pl - Calculates the Neff in the position specific profile of a sequence or profile.

=head1 SYNOPSIS

  neffseq.pl [OPTIONS] -i INFILE-GLOB -D MODEL

  OPTIONS:
    -i, --infile INFILE-GLOB+       Input file glob with alignment or sequence
    -D, --context-data MODEL        Model to be used for calculating pseudocounts
    -x, --pc-admix [0;1]            Pseudocount admixture for context-specific pseudocounts [default: 0.9]
    -h, --help                      Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @inglobs;
my @infiles;
my $D;
my $x = 0.9;


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "i|infile=s{1,}" => \@inglobs,
  "D|model=s" => \$D,
  "x|pc-admix=f" => \$x,
  "h|help" => sub { pod2usage(2); }
) or pod2usage(1);
unless ($D) { pod2usage("No model provided!"); }
unless (defined($x)) { pod2usage("No pc-admix provided!"); }
foreach my $inglob (@inglobs) {
  foreach my $f (glob($inglob)) {
    if (-f $f) { push(@infiles, $f); }
  }
}
unless (@infiles) { pod2usage("No input files provided!"); }


### Compute the average Neff ###


my $neff = 0;
my $neff_n = 0;
foreach my $file (@infiles) {
  ($_, my $cp) = tempfile(filename($file) . "XXXX", SUFFIX => ".tmp", DIR => "/tmp");
  &exec("csbuild -i $file -D $D -x $x -o $cp");
  $_ = &exec("cscp_neff -i $cp");
  if (/^Neff\s*=\s*(\d+(\.\d+)?)/m) {
    $neff += $1;
    $neff_n++;
    &exec("rm -f $cp");
  } else { die "Error executing 'cscp_neff -i $cp'!"; }
}
if ($neff_n) { $neff /= $neff_n; }
printf("%.2f\n", $neff);


### Misc ###


sub exec {
  my ($cmd) = @_;
  my $rv = `$cmd 2> /dev/null`;
  if ($?) { die "Error executing '$cmd'!"; }
  return $rv;
}
