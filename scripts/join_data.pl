#!/usr/bin/perl -w

use strict;
use warnings;
use My::Utils qw(max);

my %data;
my @files = @ARGV;
my @max = (0, 0);

foreach my $f (@files) {
	open FIN, "< $f" or die "Can't read from '$f'!";
	while (<FIN>) {
		if (/^\s*([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)/) {
			$max[0] = max($max[0], length($1));
			$max[1] = max($max[1], length($2));
			push(@{$data{$1}}, $2);
		}
	}
	close FIN;
}

my @xx = sort({ $a <=> $b } keys(%data));
foreach my $x (@xx) {
  my @flds;
  push(@flds, sprintf("%$max[0]s", $x));
  foreach my $v (@{$data{$x}}) {
	  push(@flds, sprintf("%$max[1]s", $v));
	}
  print join("   ", @flds), "\n";
}
