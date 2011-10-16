#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename qw(basename);
use File::Temp qw(tempfile);


=pod
=head1 NAME

  cstrainset.pl - Submits cstrainset jobs

=head1 SYNOPSIS

  cstrainset.pl [OPTIONS] CSTRAINSET-OPTIONS

  OPTIONS:
  -B, --basename    Basename of the training set.
  -S, --suffix      Suffix to be appended.
  -V, --vset        Create validation set.
  -h, --help        Print this help message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $HOME      = $ENV{"HOME"};
my $DBS       = "$HOME/databases";
my $CSD       = "$HOME/data/cs";
my $K4000     = "$CSD/models/lib/K4000.lib";
my $CST       = "$CSD/trainsets";

my $d;
my $e;
my $N         = 3000000;
my $s;
my $g         = 1.0;
my $n         = 4.0;
my $m         = 4.0;
my $M         = 20.0;
my $y         = 0.0;
my $x         = 0.0;
my $W         = 13;
my $D         = undef;

my $bn;
my $vset      = 1;
my $suffix;
my @args;
my $cpu       = 1;
my $pe;
my $submit    = 1;


### Initialization ###


Getopt::Long::Configure("pass_through", "no_ignore_case");
GetOptions(
  "d=s"          => \$d,
  "e=s"          => \$e,
  "N=i"          => \$N,
  "s=i"          => \$s,
  "g=f"          => \$g,
  "n=f"          => \$n,
  "m=f"          => \$m,
  "M=f"          => \$M,
  "y=f"          => \$y,
  "x=f"          => \$x,
  "D=s"          => \$D,
  "W=i"          => \$W,
  "pe=s"         => \$pe,
  "cpu=i"        => \$cpu,
  "basename=s"   => \$bn,
  "suffix=s"     => \$suffix,
  "vset!"        => \$vset,
  "submit!"      => \$submit,
  "h|help"       => sub { pod2usage(2); }
) or pod2usage(1);
@args = @ARGV;
unless ($d) { pod2usage("No database provided!"); }
unless (-d $d) { pod2usage("Database does not exist!"); }
if (defined($e) && ! -d $e) { pod2usage("Database for pseudocounts column does not exist!"); }
$d =~ s/\/$//;
if (defined($e)) { $e =~ s/\/$//; }
else { $n = $m; }
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
  

&submit;
if ($vset) { &submit(1); }


sub submit {
  my ($vs) = @_;
  my $dd = $d;
  my $ee = $e;
  my $bb = &get_basename;
  my $NN = $N;
  my $ss = $s ? $s : ($g == 1.0 ? 1 : 3);

  if ($vs) {
    $dd = &get_vset($dd);
    unless (-d $dd) { return; }
    if ($ee) {
      $ee = &get_vset($ee);
      unless (-d $ee) { return; }
    }
    $bb = &get_vset($bb);
    unless ($bb) { return; }
    $NN = 1500000;
  }

  my $ext = $g == 1.0 ? "tsq" : "tpr";
  my $out = sprintf("%s/%s_g%.2f%s_m%.1f_M%.1f_y%.1f_N%s%s", $CST, $bb, $g, 
    $ee ? sprintf("_n%.1f", $n) : "", 
    $m, $M, $y, &get_N_short($NN),
    $D ? sprintf("_D%s", basename($D)) : "",
  );
  if ($suffix) { $out .= sprintf("_%s", $suffix); }
  if ($vs && -e "$out.$ext") { print "Validation set already exists!\n"; return; }
  my $cmd = sprintf(qq/cstrainset \\
    -d $dd %s \\
    -s $ss \\
    -g $g \\
    -n $n \\
    -m $m \\
    -M $M \\
    -y $y \\
    -N $NN \\
    -x $x \\
    -W $W \\
    -o $out.$ext %s %s \\
    &> $out.log/, 
    $ee ? "-e $ee" : "", 
    $D ? "-D $D" : "", 
    join(" ", @args));
  $cmd =~ s/^\s+//mg;
  #print "$cmd\n"; exit 0;

  my ($fout, $scriptfile) = tempfile("cstrainsetXXXX", DIR => "/tmp", SUFFIX => ".sh");
  print $fout "#!/bin/bash\n";
  print $fout "export OMP_NUM_THREADS=\$NSLOTS\n";
  print $fout "$cmd\n";
  close($fout);

  system("chmod u+x $scriptfile");
  print "Job script: $scriptfile\n";
  if ($submit) { system("qsub -pe $pe $cpu $scriptfile"); } 
}

sub get_basename {
  if ($bn) { return $bn; }
  else {
    my $bb = basename($d);
    if (defined($e) && $e =~ /neff(.+)$/) { $bb = "${bb}_$1"; }
    return $bb;
  }
}

sub get_vset {
  my ($t) = @_;
  if ($t =~ s/_t(\b|_)/_v$1/) { return $t; }
  else { return ""; }
}

sub get_N_short {
  my ($NN) = @_;
  if ($NN < 1e6) { return "$NN"; }
  else { return sprintf("%.1fM", $NN / 1e6); }
}