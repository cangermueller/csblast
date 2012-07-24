#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Spec::Functions qw(catfile);
use File::Basename qw(dirname basename);

=pod1
=head1 NAME

  csfdr.pl - Returns the number of FP/TP fiven a FDR.

=head1 SYNOPSIS

  csfdr.pl [OPTIONS] -i DIR+

  OPTIONS:
    -i, --indir DIR+    Input directory of data file containing FP/TP.
    -f, --fdr [0;1]+    FDR for which the number of FP/TP is to be returned [def: 0.01 0.1 0.2].
    -h, --help          Show this help text.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @indirs;
my @infiles;
my @fdrs;


### Initialization ###


GetOptions(
  "i|indir=s{1,}" => \@indirs,
  "f|fdr=f{1,}"   => \@fdrs,
  "h|help"        => sub { pod2usage(2); }
) or die pod2usage(1);
unless (@indirs) { pod2usage("No input directories provided!"); }
unless(@fdrs) { @fdrs = (0.2, 0.1, 0.01); }


### Main ###


foreach my $dir (@indirs) {
  my $infile = -d $dir ? catfile($dir, "wtpfp.dat") : $dir;
  if (-f $infile) { push(@infiles, $infile); }
}

foreach my $fdr (@fdrs) {
  my @entries;
  foreach my $infile (@infiles) {
    my $cur_fptp = &get_fptp($infile, $fdr);
    if ($cur_fptp) { 
      my $e = {
        FP   => $cur_fptp->[0],
        TP   => $cur_fptp->[1],
        NAME => basename(dirname($infile))
      };
      push(@entries, $e);
    }
  }
  unless (@entries) { last; }
  @entries = sort({ $b->{TP} <=> $a->{TP} } @entries);
  my $max = $entries[0]->{TP};
  printf("== FDR %.2f\n", $fdr);
  printf("%3s   %8s   %8s   %8s   %8s   %s\n", "NR", "TP", "-%", "FP", "FDR", "NAME");
  my $i = 0;
  foreach my $e (@entries) {
    printf("%3d   %8.2f   %8.2f   %8.2f   %8.2f   %s\n", ++$i, 
      $e->{TP}, ($max/$e->{TP}-1)*100, $e->{FP}, $e->{FP}/($e->{FP}+$e->{TP}), $e->{NAME});
  }
}
  

sub get_fptp {
  my ($infile, $fdr) = @_;
  my @lo;
  my @up;

  open FIN, "< $infile" or die "Can't read from '$infile'!";
  while (<FIN>) {
    if (/^#/) {
      next;
    } elsif (my ($fp, $tp) = /^\s*(\d\S*)\s+(\d\S*)\s*$/) {
      my $cur_fdr = $fp / ($fp + $tp);
      if ($cur_fdr <= $fdr) {
        if (!defined($lo[0]) || $fdr - $cur_fdr < $fdr - $lo[0]) {
          @lo = ($cur_fdr, $fp, $tp);
        }
      } elsif (!defined($up[0]) || $cur_fdr - $fdr < $up[0] - $fdr) {
        @up = ($cur_fdr, $fp, $tp);
      }
    } else {
      return undef;
    }
  }
  close FIN;
  if (@lo && @up) {
    my $lo_dev = $fdr - $lo[0];
    my $up_dev = $up[0] - $fdr;
    my $w = $lo_dev / ($lo_dev + $up_dev);
    return [$lo[1] * (1 - $w) + $up[1] * $w, $lo[2] * (1 - $w) + $up[2] * $w];
  }
  return undef;
}
