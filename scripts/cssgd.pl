#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(basename);
use File::Temp qw(tempfile);
use File::Spec::Functions;
use My::Utils qw(filename);


=pod
=head1 NAME

  cssgd.pl - Submits cssgd jobs

=head1 SYNOPSIS

  cssgd.pl [OPTIONS] CSSGD-OPTIONS

  OPTIONS:
    --outdir OUTDIR         Output directory.
    --outcat OUTCAT	        Output catecory.
    --outsuffix OUTSUFFIX   Suffix to be appended.
    --repeats REPEATS       Number of parallel CSSGD runs.
    --rounds ROUNDS         Number of CSSGD rounds.
    --seed SEED             CSSGD seed.
    --pe PE                 SGE parallel environment.
    --cpu CPU               Number of CPUs to be used.
    --submit                Submit job script.
    --help                  Print this help message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $i;
my $j;
my $K       = 200;
my $P       = 3;
my $b       = 1.0;
my $c       = 1.0;
my $d       = 0.8;
my $p       = 1.0;
my $q       = 0;
my $t       = 0.001;
my $T       = 0.001;
my $E       = 2;
my $e       = 0.03;
my $D       = 1.5;
my $B       = 1000;
my $m;

my $outdir;
my $outcat  = "share";
my $outsuffix;
my $rounds  = 1;
my $repeats = 1;
my $seed    = 0;
my $submit  = 1;
my @args;

my $cpu     = 1;
my $pe;

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
  "q=i"              => \$q,
  "E=i"              => \$E,
  "e=f"              => \$e,
  "D=f"              => \$D,
  "B=i"              => \$B,
  "t=f"              => \$t,
  "T=f"              => \$T,
  "m=s"              => \$m,
  "outdir=s"         => \$outdir,
  "outcat=s"         => \$outcat,
  "outsuffix=s"      => \$outsuffix,
  "rounds=i"         => \$rounds,
  "repeats=i"        => \$repeats,
  "seed=i"           => \$seed,
  "submit!"          => \$submit,
  "cpu=i"            => \$cpu,
  "pe=s"             => \$pe,
  "h|help"           => sub { pod2usage(2); }
) or pod2usage(1);
@args = @ARGV;
unless ($i) { pod2usage("No trainset provided!"); }
unless (-f $i) { pod2usage("Trainset does not exist!"); }
unless ($j) { $j = &get_vset($i); }
unless (-f $j) { pod2usage("Validation set does not exits!"); }
unless ($pe) {
  my @pl = `qconf -spl`;
  unless (@pl) { die "No parallel environment available!"; }
  $pe = $pl[0];
  foreach my $p (@pl) {
    chomp($p);
    if ($p eq "default" || $p eq "threads.pe") {
      $pe = $p;
      last;
    }
  }
}
unless ($outdir) {
  $outdir = catdir($CSC, $outcat ? $outcat : "share");
}


### Do the job ###


if ($repeats > 1) {
  for my $r (1 .. $repeats) {
    srand($seed + $r);
    &submit(catdir($outdir, sprintf("%02d", $r), int(rand(1e6))));
  }
} else {
  &submit($outdir, $seed);
}


### Utilities ###


sub submit {
  my ($dir, $seed) = @_;
  my $outbase = catdir($dir, &get_tset($i));
  if ($outsuffix) { $outbase .= "_$outsuffix"; }

  system("mkdir -p $dir");
  my $scriptfile = catfile($dir, "cssgd.sh");
  open FOUT, "> $scriptfile" or die "Can't write to '$scriptfile'!";

  print FOUT "export OMP_NUM_THREADS=\$NSLOTS\n";
  my $cmdbase = sprintf(qq/cssgd \\
    -i $i \\
    -j $j \\
    -o CRF-VSET \\
    -O CRF-TSET \\
    -K $K \\
    -P $P \\
    -b $b \\
    -c $c \\
    -d $d \\
    -p $p \\
    -q $q \\
    -E $E \\
    -e $e \\
    -D $D \\
    -B $B \\
    -t $t \\
    -T $T \\
    --weight-center $c \\
    --weight-decay $d \\
    -m MODEL \\
    --seed $seed \\
    %s &> LOGFILE/,
    join(" ", @args));

  for my $r (1 .. $rounds) {
    my $out = $outbase;
    if ($rounds > 1) { $out = sprintf("%s_%02d", $out, $r); }
    my $cmd = $cmdbase;
    $cmd =~ s/CRF-VSET/${out}_v.crf/;
    $cmd =~ s/CRF-TSET/${out}_t.crf/;
    $cmd =~ s/LOGFILE/${out}.log/;
    if ($r > 1) { 
      my $m = sprintf("%s_%02d_v.crf", $outbase, $r - 1);
      $cmd =~ s/MODEL/$m/;
    } elsif ($m) {
      $cmd =~ s/MODEL/$m/;
    } else { 
      $cmd =~ s/\s*-m MODEL.*$//m; 
    }
    print FOUT "$cmd\n";
    print FOUT 'if [ $? -ne 0 ]; then exit $?; fi', "\n";
  }
  close(FOUT);

  system("chmod u+x $scriptfile");
  print "Job script: $scriptfile\n";
  if ($submit) { system("qsub -pe $pe $cpu $scriptfile"); } 
}

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
