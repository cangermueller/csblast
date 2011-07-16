#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


=pod
=head1 NAME

 rocx.pl - Computes the mean rocx score of given rocx curves

=head1 SYNOPSIS

 rocx.pl [OPTIONS] --dir DIR+

 OPTIONS:
 -d, --dir DIR+       Input directory list
 -c, --curve CURVE    Name of the rocx curve file in the given directory [default: rocx.dat]
 
=head1 AUTHOR

 Angermueller Christof
 angermue@in.tum.de

=cut


### Variables ###


my @dirs;
my $curve = "rocx.dat";
my @rocx;


### Initialization ###


GetOptions(
  "dir|d=s{1,}" => \@dirs,
  "curve|c=s"   => \$curve,
  "help|h"      => sub { pod2usage(2); }
) or pod2usage(1);
unless (scalar(@dirs)) { pod2usage("Input directories missing!"); }


### Compute mean rocx scores ###


foreach my $dir (@dirs) {
  if (-e "$dir/$curve") { push(@rocx, [&rocx("$dir/$curve"), $dir]); }
}
@rocx = sort({ $b->[0] <=> $a->[0] } @rocx);
printf("%3s\t%10s\t%s\n", "NR", "ROCX", "NAME"); 
for my $i (0 .. $#rocx) {
  printf("%3d\t%10.5f\t%s\n", $i + 1, $rocx[$i]->[0], $rocx[$i]->[1]);
}


sub rocx {
  my ($file) = @_;
  my $rocx = 0;
  my $x_last = 0;
  my $y_last = 1;
  open(FIN, "< $file") or die "$file: $!";
  while (<FIN>) {
    if (my ($x, $y) = /^([^\#]\S*)\s+(\S+)\s*$/) {
      $rocx += ($x - $x_last) * 0.5 * ($y + $y_last);
      ($x_last, $y_last) = ($x, $y);
    }
  }
  return $rocx;
}

