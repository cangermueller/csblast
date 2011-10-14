#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(basename);
use File::Temp qw(tempfile);
use My::Utils qw(filename);


=pod
=head1 NAME

  cssgd.pl - Submits cssgd jobs

=head1 SYNOPSIS

  cssgd.pl [OPTIONS] CSSGD-OPTIONS

  OPTIONS:
  -S, --suffix SUFFIX     Suffix to be appended
  -C, --cat CAT	          Catecory of the model
  -R, --rounds ROUNDS     Number of CSSGD rounds
  -h, --help              Print this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $i;
my $j;
my $K       = 50;
my $P       = 2;
my $b       = 1.0;
my $c       = 1.0;
my $d       = 1.0;
my $p       = 1.0;
my $q       = 0;
my $t       = 0.001;
my $T       = 0.001;
my $E       = 2;
my $e       = 0.05;
my $D       = 2.0;
my $B       = 1000;
my $m;

my $suffix;
my $cat     = "share";
my $rounds  = 1;
my $repeats = 1;
my $seed    = 0;
my $submit  = 1;
my @args;

my $cpu     = 1;
my $pe;
my $queue   = undef;

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
  "suffix=s"       => \$suffix,
  "cat=s"          => \$cat,
  "rounds=i"       => \$rounds,
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


### Do the job ###


if ($repeats > 1) {
  for my $r (1 .. $repeats) {
    srand($seed + $r);
    &submit(sprintf("%s/%s/%02d", $CSC, $cat, $r), int(rand(1e6)));
  }
} else {
  &submit(sprintf("%s/%s", $CSC, $cat), $seed);
}


### Utilities ###


sub submit {
  my ($dir, $seed) = @_;
  my $outbase = sprintf("$dir/%s_K%d", &get_tset($i), $K);
  if ($m) { $outbase .= sprintf("_m%s", basename($m)); }
  if ($suffix) { $outbase .= sprintf("_%s", $suffix); }
  my ($fout, $scriptfile) = tempfile("cssgdXXXX", DIR => "/tmp", SUFFIX => ".sh");

  print $fout "export OMP_NUM_THREADS=\$NSLOTS\n";
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
    -m MODEL \\
    --seed $seed \\
    %s &> LOGFILE/,
    join(" ", @args));

  print $fout "mkdir -p $dir\n";
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
    print $fout "$cmd\n";
    print $fout 'if [ $? -ne 0 ]; then exit $?; fi', "\n";
  }
  close($fout);

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
