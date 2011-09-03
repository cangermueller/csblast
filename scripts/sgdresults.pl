#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(dirname basename);
use My::Utils qw(is_numeric);


=pod
=head1 NAME

  sgdresults.pl - Creates results file for the given SGD output files.

=head1 SYNOPSIS

  sgdresults.pl [OPTIONS] --infile INFILE

  OPTIONS:
  -i, --infile INFILE+    SGD output files.
  -b, --bench DIR         Directory with bench results.
  -p, --param PARAM+      Parameters to be listed
  -s, --sort FIELD        Field to be used for sorting [default: ll-val]
                          FIELD = name|param|ll-train|ll-val|neff|rocx
  -c, --complete          Only list complete SGD runs [default: false]
  -h, --help              Show this help message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @infiles;
my $benchdir;
my @params;
my $sort;
my $complete;
my @results;


### Initialization ###


GetOptions(
  "i|infile=s{1,}" => \@infiles,
  "b|bench=s"      => \$benchdir,
  "p|param=s{1,}"  => \@params,
  "s|sort=s"       => \$sort,
  "c|complete!"    => \$complete,
  "h|help"         => sub { pod2usage(2); }
) or die pos2usage(1);
unless (@infiles) { pod2usage("Not input files provided!"); }
unless ($benchdir) {
  $benchdir = sprintf("%s/%s", $ENV{CSBENCH}, basename(dirname(abs_path($infiles[0]))));
  unless (-d $benchdir) { $benchdir = undef; }
}
unless ($sort) { $sort = ($benchdir ? "rocx" : "ll-val"); }

foreach my $infile (@infiles) {
	if (-e $infile) { 
		my $r = &get_results($infile);
		if ($r) { push(@results, $r); }
	}
}
if ($sort eq "name") { @results = sort({ $a->{NAME} cmp $b->{NAME} } @results); }
elsif ($sort eq "param" && scalar(@params)) { @results = sort({ 
			for my $i (0 .. $#{$a->{PARAMS}}) {
				if ($a->{PARAMS}->[$i] != $b->{PARAMS}->[$i]) { return $a->{PARAMS}->[$i] cmp $b->{PARAMS}->[$i]; }
			}
			0; } @results); }
elsif ($sort eq "ll-train") { @results = sort({ &comp_num($b->{SGD}->{LLT}, $a->{SGD}->{LLT}) } @results); }
elsif ($sort eq "neff") { @results = sort({ $b->{SGD}->{NEFF} <=> $a->{SGD}->{NEFF} } @results); }
elsif ($sort eq "rocx") { @results = sort({ $b->{ROCX}->{ROCX} <=> $a->{ROCX}->{ROCX} } @results); }
else { @results = sort({ &comp_num($b->{SGD}->{LLV}, $a->{SGD}->{LLV}) } @results); }

my @len = ((0) x (scalar(@params) + 1), 1, 3, 3, 8, 8, 8, 8, 8);
foreach my $r (@results) {
	if (length($r->{NAME}) > $len[0]) { $len[0] = length($r->{NAME}); }
	for my $i (0 .. $#{$r->{PARAMS}}) {
		if (length($r->{PARAMS}->[$i]) > $len[$i + 1]) { $len[$i + 1] = length($r->{PARAMS}->[$i]); }
	}
}
my $sep = " " x 2;
my $format = "%-$len[0]s";
foreach my $l (@len[1 .. $#len]) {
	$format .= sprintf("%s%%%ds", $sep, $l);
}
$format .= "\n";
printf($format, "NAME", @params, "C", "EL", "ES", "LL-TR", "LL-VAL", "NEFF", "ROCX", "ROCX_P");
my $len_tot = 0;
foreach my $l (@len) { $len_tot += $l + length($sep); }
printf("%s\n", "-" x $len_tot);
foreach my $r (@results) {
	printf("%-$len[0]s", $r->{NAME});
	for my $i (0 .. $#params) {
		printf("$sep%$len[$i + 1]s", $r->{PARAMS}->[$i]);
	}
	printf("$sep%1s$sep%3d$sep%3d$sep%8.4f$sep%8.4f$sep%8.4f$sep%s$sep%8s\n", $r->{SGD}->{COMPLETE} ? "T" : "F", 
		$r->{SGD}->{EL}, $r->{SGD}->{ES}, $r->{SGD}->{LLT}, $r->{SGD}->{LLV}, $r->{SGD}->{NEFF}, 
		$r->{ROCX}->{ROCX} ? sprintf("%8.4f", $r->{ROCX}->{ROCX}) : "",
    $r->{ROCX}->{P});
}


sub get_results {
	my ($log) = @_;
	if (-e $log) {
		my %r;
		$r{NAME}       = &get_name($log);
		$r{PARAMS}     = &get_params($log);
		$r{SGD}        = &get_sgd($log);
		if ($r{SGD}) {
			$r{ROCX}     = &get_rocx($log);
			return \%r;
		}
	}
	return undef;
}

sub get_name {
	my ($log) = @_;
	$log =~ s/\.log*//;
	return $log;
}

sub get_params {
	my ($log) = @_;
	my @v;
	foreach my $p (@params) {
		if ($log =~ /_$p(\d+(?:\.\d+)?)/) { push(@v, $1); }
		else { push(@v, ""); }
	}
	return \@v;
}

sub get_sgd {
	my ($log) = @_;
	my $c = 0;
	my @last;
	my @best;
	open FIN, "< $log" or die "Can't read from '$log'!";
	while (<FIN>) {
		if (my @line = /^(\d+)\s+\[=+\]\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s*$/) { 
			@last = @line;
			if (scalar(@best) == 0 || $line[3] >= $best[3]) { @best = @line; }
		} elsif (scalar(@last)) { 
			if (/^-----/) { $c = 1; }
			last; 
		}
	}
	close FIN;
	if (scalar(@last) == 0 || $c == 0 && defined($complete)) { return 0; }
	my %sgd;
	$sgd{EL}       = $last[0];
	$sgd{ES}       = $best[0];
	$sgd{LLT}      = $last[1];
	$sgd{LLV}      = $best[3];
	$sgd{NEFF}     = $best[4];
	$sgd{COMPLETE} = $c;
	return \%sgd;
}

sub get_rocx {
	my ($log) = @_;
  my $rocx = {ROCX => 0, P => ""};
	unless ($benchdir) { return $rocx; }
	my $name = &get_name($log);
	$_ = `rocx.pl -d $benchdir/$name*`;
	if ($?) { die "Error calling rocx.pl!"; }
	foreach (split(/\n/)) {
    if (/^\s*(\d+)\s+(\S+)\s+(\S+)\s*$/) {
			if ($2 > $rocx->{ROCX}) { 
        $rocx->{ROCX} = $2;
        my $name = $3;
        my $p = "";
        if ($name =~ /_t\.crf/) { $p = "t"; }
        elsif ($name =~ /_v\.crf/) { $p = "v"; }
        if ($name =~ /_([^_]+)$/) { $p .= $1; }
        $rocx->{P} = $p;
      }
		}
	}
	return $rocx;
}

sub comp_num {
  my ($a, $b) = @_;
  unless (is_numeric($a)) { return -1; }
  unless (is_numeric($b)) { return 1; }
  return $a <=> $b;
}
