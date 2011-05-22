#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use List::Util qw(shuffle);
use Time::HiRes qw(gettimeofday tv_interval);
use My::Config;
use File::Basename qw(basename);

=pod

=head1 NAME
 
  timer.pl - Measure time for CSBLAST pseudocounts generation before a BLAST call.

=head1 SYNOPSIS

  timer.pl [OPTIONS] -D CONTEXT-DATA+

  OPTIONS:
 
    -D, --context-data FILE+  Model with context-data.
    -d, --dir FILE            Directory with sequence files.
    -x, --pc-admix ]0,1]      Pseudocount admix for context-specific pseudocounts [def: 0.90]
    -z, --pc-neff [1,inf[     Target Neff for pseudocounts admixture [def: 0.00]
    -e, --ext EXT             Sequence files extension [def: seq].
    -n, --num [1,inf[         Number of sequence to be used for time measurement.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut

 
### Variables ###


my %cfg           = new My::Config;
my @context_data;
my $pc_admix      = 0.9;
my $pc_neff       = 0.0;
my $dir           = "$cfg{'dbs_path'}/scop20_1.75_opt";
my $ext           = "seq";
my $num           = 10;
my @seqs;
my @time;


### Initialization ###


GetOptions(
  "D|context-data=s{1,}" => \@context_data,
  "d|dir=s"              => \$dir,
  "e|ext=s"              => \$ext,
  "x|pc-admix=f"         => \$pc_admix,
  "z|pc-neff=f"          => \$pc_neff,
  "n|num=i"              => \$num,
  "h|help"               => sub { pod2usage(2); }
) or pod2usage(1);
unless (@context_data) { pod2usage("No context-data provided!"); }
foreach my $cd (@context_data) {
  unless (-f $cd) { pod2usage("'$cd' does not exist!"); }
}
unless (-d $dir) { pod2usage("Directory with sequence files does not exist!"); }


### Time measurement ###


@seqs = glob("$dir/*.$ext");
@seqs = shuffle(@seqs);
if ($num > scalar(@seqs)) { $num = scalar(@seqs); }

printf("%-30s : %d\n", "context-data", scalar(@context_data));
printf("%-30s : %.2f\n", "pc-admix", $pc_admix);
printf("%-30s : %.2f\n", "pc-neff", $pc_neff);
printf("%-30s : %d\n", "samples", $num);

foreach my $cd (@context_data) {
  my $t0 = [gettimeofday];
  for my $i (1 .. $num) {
    system(sprintf("csblast -D %s -i %s -x %.2f -z %.2f --blast-path %s --emulate &> /dev/null",
        $cd, $seqs[$i - 1], $pc_admix, $pc_neff, $cfg{"blast_path"}));
    if ($? != 0) { die "Error calling csblast: $!"; }
  }
  my $t1 = [gettimeofday];
  my $td = tv_interval($t0, $t1) * 1000;
  push(@time, [$td, basename($cd)]);
}

printf("%s\n", "-" x 60);
@time = sort({ $a->[0] <=> $b->[0] } @time);
printf("%12s\t%12s\t%s\n", "TIME", "TIME/SAMPLE", "CONTEXT-DATA");
foreach my $t (@time) {
  printf("%12d\t%12d\t%s\n", $t->[0], $num > 0 ? $t->[0] / $num : 0.0, $t->[1]);
}
