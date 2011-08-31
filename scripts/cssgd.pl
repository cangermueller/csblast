#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(basename);
use My::Utils qw(filename);


=pod
=head1 NAME

  cssgd.pl - Submits cssgd jobs

=head1 SYNOPSIS

  cssgd.pl [OPTIONS] CSSGD-OPTIONS

  OPTIONS:
  -S, --suffix SUFFIX     Suffix to be appended.
  -C, --cat CAT	          Catecory of the model.
  -h, --help              Print this help message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $i;
my $j;
my $K       = 4000;
my $P       = 2;
my $b       = 10.0;
my $c       = 10.0;
my $d       = 1.0;
my $p       = 2.5;
my $q       = 0.002;
my $Q       = 30;
my $T       = 0.001;
my $e       = 0.001;
my $B       = 1000;
my $m;

my $suffix;
my $cat     = "share";
my @args;

my $pe      = "threads.pe 4";
my $queue   = undef;
# my $queue   = "normal,opteron2354,opteron8380,quadcore96gb,x2270,x4100,small.q";

my $HOME    = $ENV{"HOME"};
my $CSD     = "$HOME/data/cs";
my $CSC     = "$CSD/models/crf";


### Initialization ###


Getopt::Long::Configure("pass_through", "no_ignore_case");
GetOptions(
  "i=s"              => \$i,
  "j=s"              => \$j,
  "K=i"              => \$K,
  "P=i"              => \$P,
  "b=f"              => \$b,
  "c=f"              => \$c,
  "d=f"              => \$d,
  "p=f"              => \$p,
  "q=f"              => \$q,
  "Q=i"              => \$Q,
  "T=f"              => \$T,
  "e=f"              => \$e,
  "B=i"              => \$B,
  "m=s"              => \$m,
  "S|suffix=s"       => \$suffix,
  "C|cat=s"          => \$cat,
  "h|help"           => sub { pod2usage(2); }
) or pod2usage(1);
@args = @ARGV;
unless ($i) { pod2usage("No trainset provided!"); }
unless (-f $i) { pod2usage("Trainset does not exist!"); }
unless ($j) { $j = &get_vset($i); }
unless (-f $j) { pod2usage("Validation set does not exits!"); }


### Command composition ###


my $out = sprintf("%s/%s/%s_K%d_b%.1f_c%.1f_p%.1f", $CSC, $cat, &get_tset($i), $K, $b, $c, $p, $q);
if ($m) { $out .= sprintf("_m%s", basename($m)); }
if ($suffix) { $out .= sprintf("_%s", $suffix); }
my $crf_tset = sprintf("%s_t.crf", $out);
my $crf_vset = sprintf("%s_v.crf", $out);
my $crf_log = sprintf("%s.log", $out);
if (-e $crf_tset || -e $crf_vset) { die "CRF already exists!"; }

my $cmd = sprintf(
"qsub -pe $pe %s -o $crf_log -e $crf_log -b y " .
"cssgd -i $i -j $j -K $K -P $P -b $b -c $c -d $d -p $p -q $q -Q $Q " . 
"-T $T -e $e -B $B -o $crf_vset -O $crf_tset %s %s",
  $queue ? "-q '$queue'": "", 
  $m ? "-m $m" : "",
  join(" ", @args));

#print "$cmd\n"; exit 0;
system("mkdir -p $CSC/$cat");
system("$cmd");
if ($?) { die "Error calling cssgd!"; }



sub get_vset {
  my ($t) = @_;
  if ($t =~ s/_t_/_v_/ && $t =~ s/_N\d+(\.\d+M)?/_N*/) {
    @_ = glob($t);
    if (scalar(@_)) { return $_[0]; }
  }
  return "";
}

sub get_tset {
  my ($tset) = @_;
  $tset = filename($tset);
  $tset =~ s/_t_/_/;
  return $tset;
}
