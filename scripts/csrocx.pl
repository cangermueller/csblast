#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(catdir);
use File::Basename qw(dirname basename);


=pod
=head1 NAME

  csrocx.pl - Computes the mean rocx score of given rocx curves

=head1 SYNOPSIS

  csrocx.pl [OPTIONS] -i DIR+

  OPTIONS:
    -i, --indir DIR+          Input directory of data files containing ROCX data.
    -h, --help                Print this help text.
 
=head1 AUTHOR

 Angermueller Christof
 angermue@in.tum.de

=cut


### Variables ###


my @indirs;
my @rocx;


### Initialization ###


GetOptions(
  "i|indir=s{1,}" => \@indirs,
  "h|help"        => sub { pod2usage(2); }
  ) or pod2usage(1);
unless (scalar(@indirs)) { pod2usage("No input directories provided!"); }


### Main ###


foreach my $i (@indirs) {
  my $infile = -d $i ? catdir($i, "rocx.dat") : $i;
  if (-f $infile) { 
    my $cur_rocx = &get_rocx($infile);
    if ($cur_rocx) { push(@rocx, [$cur_rocx, basename(dirname($infile))]); }
  }
}
@rocx = sort({ $b->[0] <=> $a->[0] } @rocx);
printf("%3s   %10s   %6s   %s\n", "NR", "ROCX", "-%", "NAME"); 
unless (@rocx) { exit 0; }
my $max = $rocx[0]->[0];
for my $i (0 .. $#rocx) {
  printf("%3d   %10.5f   %6.2f  %s\n", $i + 1, 
    $rocx[$i]->[0], ($max/$rocx[$i]->[0]-1)*100, $rocx[$i]->[1]);
}


sub get_rocx {
  my ($file) = @_;
  my $rocx = 0;
  my $x_last = 0;
  my $y_last = 1;
  open(FIN, "< $file") or die "$file: $!";
  while (<FIN>) {
    if (/^#/) { next; }
    elsif (my ($x, $y) = /^\s*(\d\S*)\s+(\d\S*)\s*$/) {
      $rocx += ($x - $x_last) * 0.5 * ($y + $y_last);
      ($x_last, $y_last) = ($x, $y);
    } else { return undef; }
  }
  return $rocx;
}
