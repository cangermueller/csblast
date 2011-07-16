#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename qw(dirname);

if (scalar(@ARGV) < 2) {
	print STDERR "error_rate.pl ERR DIRECTORY+";
	exit 1;
}

my $err = shift(@ARGV);
my @rocx = @ARGV;
my $file = "wtpfp.dat";
my @fptp;

foreach my $r (@rocx) {
	$r =~ s/\/$//;
	$r .= "/$file";
	my @best;
	my $best_err;
	if (! -e $r) { next; }
	open FIN, "< $r" or die "Can't read from '$r'!";
	while (<FIN>) {
		if (my ($fp, $tp) = /^\s*([^#]\S+)\s*(\S+)\s*$/) {
			my $e = $fp / ($fp + $tp);
			if (!defined($best_err) || abs($err - $best_err) > abs($err - $e)) {
				@best = ($fp, $tp);
				$best_err = $e;
			}
		}
	}
	close FIN;
	push(@fptp, \@best);
}

printf("%10s   %10s   %6s   %s\n", "FP", "TP", "GAIN", "DIR");
if (@fptp) {
	my $tp = $fptp[0]->[1];
	for my $i (0 .. $#fptp) {
		my $e = $fptp[$i];
		printf("%10.1f   %10.1f   %6.2f   %s\n", $e->[0], $e->[1], ($e->[1] - $tp) / $tp * 100, $rocx[$i]);
	}
}
				
