#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(catdir);


=pod
=head1 NAME

 rocx.pl - Computes the mean rocx score of given rocx curves

=head1 SYNOPSIS

 rocx.pl [OPTIONS] --dir DIR+

 OPTIONS:
 -d, --dir DIR+         Input directory list
 -b, --basename NAME    Basename of rocx files [def: ]
 
=head1 AUTHOR

 Angermueller Christof
 angermue@in.tum.de

=cut


### Variables ###


my @dirs;
my $basename;
my @rocx;


### Initialization ###


GetOptions(
  "d|dir=s{1,}"      => \@dirs,
  "b|basename=s"     => \$basename,
  "h|help"           => sub { pod2usage(2); }
  ) or pod2usage(1);
unless (scalar(@dirs)) { pod2usage("Input directories missing!"); }


### Compute mean rocx scores ###


my $rocx_file = sprintf("%srocx.dat", $basename ? "${basename}_" : "");
foreach my $dir (@dirs) {
  my $path = catdir($dir, $rocx_file);
  if (-f $path) { push(@rocx, [&rocx($path), $dir]); }
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

