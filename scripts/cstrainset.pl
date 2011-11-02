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
my $W         = 13;
my $N         = 3000000;
my $s;
my $g;
my $D;
my $u         = 3.5;
my $U         = 20.0;
my $v         = $u;
my $V         = $U;
my $j         = 0.0;
my $J         = 0.0;
my $y         = 0.0;
my $R         = 0;

my $basename;
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
  "W=i"          => \$W,
  "N=i"          => \$N,
  "s=i"          => \$s,
  "g=f"          => \$g,
  "D=s"          => \$D,
  "u=f"          => \$u,
  "U=f"          => \$U,
  "v=f"          => \$v,
  "V=f"          => \$V,
  "j=f"          => \$j,
  "J=f"          => \$J,
  "y=f"          => \$y,
  "R=i"          => \$R,
  "pe=s"         => \$pe,
  "cpu=i"        => \$cpu,
  "basename=s"   => \$basename,
  "suffix=s"     => \$suffix,
  "vset!"        => \$vset,
  "submit!"      => \$submit,
  "h|help"       => sub { pod2usage(2); }
) or pod2usage(1);
@args = @ARGV;
unless ($d) { pod2usage("No database provided!"); }
unless (-d $d) { pod2usage("Database does not exist!"); }
$d =~ s/\/$//;
if (!defined($s) && !defined($g)) {
  $s = 1.0;
  $g = 1.0;
} elsif (!defined($s)) {
  $s = $g == 1.0 ? 1 : 3;
} else {
  $g = $s == 1 ? 1.0 : 0.0;
}
if ($j == 0.0) {
  $v = $u;
  $V = $U;
}
if ($j > 0.0 && $J == 0.0) { $J = 3.0; }

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
  my $dd;
  my $bn;
  my $NN;

  if ($vs) {
    $dd = &get_vset($d);
    unless (-d $dd) { return; }
    $bn = &get_basename($dd);
    $NN = 1500000;
  } else {
    $dd = $d;
    $bn = &get_basename($dd);
    $NN = $N;
  }

  my $out = sprintf("%s/%s_N%s_g%.2f_u%.1f_U%.1f", 
    $CST, $bn, &get_N_short($NN), $g, $u, $U);
  if ($R) { $out .= sprintf("_R%d", $R); }
  if ($j > 0) { $out .= sprintf("_v%.1f_V%.1f_j%.1f_J%.1f", $v, $V, $j, $J); }
  if ($y) { $out .= sprintf("_y%.1f", $y); }
  if ($D) { $out .= sprintf("_D%s", basename($D)); }
  if ($suffix) { $out .= sprintf("_%s", $suffix); }
  my $ext = $s == 1.0 ? "tsq" : "tpr";
  if ($vs && -e "$out.$ext") { print "Validation set already exists!\n"; return; }

  my $cmd = sprintf(qq/cstrainset \\
    -d $dd \\
    -W $W \\
    -N $NN \\
    -s $s \\
    -g $g \\
    -u $u \\
    -U $U \\
    -v $v \\
    -V $V \\
    -j $j \\
    -J $J \\
    -y $y \\
    -R $R \\
    -o $out.$ext %s %s \\
    &> $out.log/, 
    $D ? "-D $D" : "", 
    join(" ", @args));
  $cmd =~ s/^\s+//mg;
  # print "$cmd\n"; return;

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
  my ($db) = @_;
  if ($basename) { return $basename; }
  else { return basename($db); }
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
