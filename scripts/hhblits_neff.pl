#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use List::Util qw(shuffle);
use File::Basename qw(basename);
use My::Utils qw(min max remove_suffix);


=pod
=head1 NAME
  
  hhblits_neff.pl - Visualizes the average Neff sequences per HHblits round.

=head1 SYNOPSIS

  hhblits_neff.pl [OPTIONS] -d DB

  OPTIONS:

    -d, --db DIR             Database directory
    -o, --outfile            Output file [default: DIR_neff.pdf]
    -n, --num [1;inf[        Number of samples to be drawn [default: all]
    -k, --keep               Keep plot files [default: false]
    -h, --help               Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $db;
my $outfile;
my $nsamples;
my $keep = 1;

my $basename;
my @samples;
my @neff;
my @count;
my $max_round = 0;

my %opts = (
  WIDTH           => 12,
  HEIGHT          => 7,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 4,
  CURVETYPE       => "lines",
  LINETYPE        => 1,
  COLORS          => [
    "#FF0000",  # red
    "#32CD32",  # limegreen
    "#1E90FF",  # dodgerblue
    "#FFA500",  # orange
    "#C71585",  # purple
    "#A0522D",  # sienna
    "#00CED1",  # darkturquoise 
    "#8470ff",  # slate blue
    "#008080",  # teal
    "#FF00FF",  # fuchsia
  ]
);


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "d|db=s"      => \$db,
  "o|outfile=s" => \$outfile,
  "n|num=i"     => \$nsamples,
  "k|keep!"     => \$keep,
  "h|help"      => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($db) { pod2usage("No database provided!"); }
$db =~ s/\/$//;
unless ($outfile) { $outfile = basename($db) . "_neff.pdf"; }
$basename = &remove_suffix($outfile);


### Compute statistics ###


print "Globbing sequence files ...\n";
@samples = glob("$db/*\.seq");
@samples = shuffle(@samples);
if ($nsamples) { $nsamples = min($nsamples, scalar(@samples)); }
else { $nsamples = scalar(@samples); }
printf("Computing the Neff of %d samples ...\n", $nsamples);
for my $i (0 .. $nsamples - 1) {
  my $s = $samples[$i];
  my $n = &neff($s);
  unless ($n) { next; }
  $neff[0] += $n;
  $count[0]++;
  $s =~ s/\.seq$//;
  foreach my $a3m (glob("${s}_*\.a3m")) {
    if ($a3m =~ /_(\d{1,2})\.a3m$/) {
      my $r = $1;
      my $n = &neff($a3m);
      if ($n) {
        $max_round = max($max_round, $r);
        $neff[$r] += $n;
        $count[$r]++;
      }
    }
  }
  if (($i + 1) % 100 == 0) { printf("   %d samples\n", $i + 1); }
}


### Create the plot ###


open FOUT, "> $basename.dat" or die "Can't write to '$basename.dat'!";
for my $round (0 .. $max_round) {
  printf(FOUT "%3d   %8.2f   %8d\n", $round, $neff[$round] / $count[$round], $count[$round]);
}
close FOUT;
my $cmd = &plot;
open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) {
  &exec("rm -f $basename.gp");
  &exec("rm -f $basename.dat");
}
&exec("ps2pdf $basename.ps $basename.pdf");
&exec("rm -f $basename.ps");
&exec("pdfcrop $basename.pdf $outfile");
print "Done!\n";


### Functions ###


sub neff {
  my ($file) = @_;
  $_ = `hhmake -i $file -o /dev/null 2> /dev/null`;
  if (/exp\(entropy\) = (\d\S*)/) { return $1; }
}

sub plot {
  my ($outfile) = @_;
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$basename.ps"
    set key right top
    set title "Average Neff per round: @{[&basename($db)]}"
    set grid
    set xlabel "Round"
    set xtics add 0, 1
    set ylabel "Neff"
    set ytics nomirror
    set y2label "Number of alignments"
    set y2range [0:${\(1.4*$count[0])}]
    set y2tics
    plot @{[&plot_cmd(["Neff", "#Alignments"])]}/;

  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub plot_cmd {
  my ($labels) = @_;
  my @cmds;
  for (my $i = $#{$labels}; $i >= 0; $i--) {
    push(@cmds, sprintf(qq/"$basename.dat" using 1:%d axes x1y%d title "$labels->[$i]" with $opts{CURVETYPE} / . 
        qq/lw $opts{LINEWIDTH} lt $opts{LINETYPE} lc rgb "$opts{COLORS}->[$i]"/,
        $i + 2, $i ? 2 : 1));
  }
  return join(", ", @cmds);
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { print STDERR "Error executing '$cmd'!\n"; }
  return $?;
}
