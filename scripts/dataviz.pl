#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename dirname);
use My::Utils qw(remove_suffix min max);
use POSIX;


=pod
=head1 NAME

  dataviz.pl - Visualizes a datafile: first column = x axis, second column = y axis, third column = y2 axis

=head1 SYNOPSIS

  dataviz.pl [OPTIONS] -i INFILE+

  OPTIONS:
    -i, --infile INFILE+      Input data files
    -o, --outfile OUTFILE     Output file [default: plot.pdf]
    -t, --title TITLE         Title of the plot [default: ]
    -l, --label LABEL+        Labels to be used instead of directory names
    -L, --legend              Show legend [default: false]
    -g, --legend-pos POS      Position of the legende [default: center right]
    -x, --xlabel XLABEL       Label of the x-axis [default: x]
    -y, --ylabel YLABEL       Label of the y-axis [default: y]
    -Y, --y2label Y2LABEL     Label of the y2-axis [default: ]
    -u, --xmin XMIN           Minimum value of the x axis
    -U, --xmax XMAX           Maximum value of the x axis
    -v, --ymin YMIN           Minimum value of the y axis
    -V, --ymax YMAX           Maximum value of the y axis
    -w, --y2min Y2MIN         Minimum value of the y2 axis
    -W, --y2max Y2MAX         Maximum value of the y2 axis
    -c, --xtics XTICS         Gnuplot xtics format
    -C, --ytics YTICS         Gnuplot ytics format
    -s, --xtic-step ]0;inf[   Step size of xtics
    -j, --xlog                Use log scale for the xaxis
    -J, --ylog                Use log scale for the yaxis
    -m, --max                 Highlight maximum [default: true]
    -k, --keep                Keep data and plot files [default: false]
    -G, --group [1;inf[       Member of the same group
    -h, --help                Show this help message


=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @infiles;
my $outfile      = "./rocx.pdf";
my $title        = "";
my @labels;
my $show_legend;
my $legend_pos   = "top right";
my $xlabel       = "x";
my $ylabel       = "y";
my $y2label;
my @xrange;
my @yrange;
my @y2range;
my $xtics;
my $ytics;
my $xtics_step;
my $group_size   = 1;
my $max          = 1;
my $xlog;
my $ylog;
my $keep;

my @entities;
my @range;
my $max_ncols;
my @xtics;
my $y2axis;

my %opts = (
  WIDTH           => 10,
  HEIGHT          => 10,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 4,
  LINEWIDTHT      => 3,
  LINETYPEY       => 1,
  LINETYPEY2      => 3,
  POINTSIZE       => 2,
  POINTSIZET      => 1,
  POINTTYPE       => 7,
  MAXCOLOR        => "#FF0000",
  COLORS          => [
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


Getopt::Long::Configure("pass_through", "no_ignore_case");
GetOptions(
  "i|infile=s{1,}" => \@infiles,
  "o|outfile=s"    => \$outfile,
  "t|title=s"      => \$title,
  "l|label=s{1,}"  => \@labels,
  "L|legend!"      => \$show_legend,
  "g|legend-pos=s" => \$legend_pos,
  "x|xlabel=s"     => \$xlabel,
  "y|ylabel=s"     => \$ylabel,
  "Y|y2label=s"    => \$y2label,
  "u|xmin=f"       => \$xrange[0],
  "U|xmax=f"       => \$xrange[1],
  "v|ymin=f"       => \$yrange[0],
  "V|ymax=f"       => \$yrange[1],
  "w|y2min=f"      => \$y2range[0],
  "W|y2max=f"      => \$y2range[1],
  "c|xtics=s"      => \$xtics,
  "C|ytics=s"      => \$ytics,
  "s|xtics-step=f" => \$xtics_step,
  "m|max!"         => \$max,
  "G|group=i"      => \$group_size,
  "j|xlog!"        => \$xlog,
  "J|ylog!"        => \$ylog,
  "k|keep!"        => \$keep,
  "h|help"        => sub { pod2usage(2); }
) or die pos2usage(1);
unless (@infiles) { pod2usage("No input files provided!"); }
if (scalar(@labels) || scalar(@infiles) > 1) { $show_legend = 1; }


### Create the plot ###


&get_entities;
$y2axis = defined($y2label) && $max_ncols > 2;
print "Creating '$outfile'...\n";

for (1 .. $max_ncols) { push(@range, [&INT_MAX, -&INT_MAX]); }
foreach my $e (@entities) {
  for my $i (0 .. $#{$e->{RANGE}}) {
    $range[$i]->[0] = min($range[$i]->[0], $e->{RANGE}->[$i]->[0]);
    $range[$i]->[1] = max($range[$i]->[1], $e->{RANGE}->[$i]->[1]);    
  }
}
unless (defined($xrange[0])) { $xrange[0] = $range[0]->[0]; }
unless (defined($xrange[1])) { $xrange[1] = $range[0]->[1]; }
unless (defined($yrange[0])) { $yrange[0] = $range[1]->[0]; }
unless (defined($yrange[1])) { $yrange[1] = $range[1]->[1] + 0.1 * ($range[1]->[1] - $range[1]->[0]); }
if ($y2axis) {
  unless (defined($y2range[0])) { $y2range[0] = $range[2]->[0]; }
  unless (defined($y2range[1])) { $y2range[1] = $range[2]->[1]; }
}
unless ($xtics) {
  unless ($xtics_step) { $xtics_step = $entities[0]->{VALUES}->[1]->[0] - $entities[0]->{VALUES}->[0]->[0]; }
  $xtics = "$xrange[0], $xtics_step";
}
unless ($ytics) {
  $ytics = "";
}

for my $i (0 .. $#{$entities[0]->{VALUES}}) {
  my $v = $entities[0]->{VALUES}->[$i]->[0];
  if ($v >= $xrange[0] && $v <= $xrange[1]) { push(@xtics, $v); }
}

my $basename = &remove_suffix($outfile);
my $cmd = &plot("$basename.ps");
open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) { &exec("rm -f $basename.gp"); }
&exec("ps2pdf $basename.ps $basename.pdf");
&exec("rm -f $basename.ps");
&exec("pdfcrop $basename.pdf $outfile");
print "Done!\n";


### Plotting ###


sub plot {
  my ($outfile) = @_;
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$outfile"

    set title "$title"
    set key $legend_pos
    set border 15 back
    set grid

    set xlabel "$xlabel"
    set xrange [$xrange[0]:$xrange[1]]
    set xtics $xtics/;
  if ($xlog) { $cmd .= qq/
    set log x/;
  }

  $cmd .= qq/
    set ylabel "$ylabel"
    set yrange[$yrange[0]:$yrange[1]]
    set ytics $ytics nomirror/;

  if ($xlog) { $cmd .= qq/
    set log y/;
  }

  if ($y2axis) { $cmd .= qq/
    set y2label "$y2label"
    set y2range[$y2range[0]:$y2range[1]]
    set y2tics/;
  }

  $cmd .= qq/
    plot @{[&cmd_curves]}
    /;

  $cmd =~ s/^\s+//mg;
  # print "$cmd\n"; exit 0;
  return "$cmd\n";
}

sub cmd_curves {
  my @curves;
  if ($y2axis) {
    for my $i (0 .. $#entities) {
      my $e = $entities[$i];
      if ($e->{NCOLS} > 2) {
        push(@curves, qq/"$e->{FILE}" using 1:3 axes x1y2 with linespoints lt $opts{LINETYPEY2} lw $opts{LINEWIDTHT} / . 
                      qq/lc rgb "$opts{COLORS}->[$i]" pt $opts{POINTTYPE} ps $opts{POINTSIZET} notitle/);
      }
    }
  }
  for (my $i = 0, my $j = 0; $i <= $#entities; $i += $group_size, $j++) {
    for (my $g = $group_size - 1; $g >= 0; $g--) {      
      if ($i + $g > $#entities) { next; }
      my $e = $entities[$i + $g];
      my $lt = $opts{LINETYPEY} + $g;      
      push(@curves, qq/"$e->{FILE}" using 1:2 with linespoints lt $lt lw $opts{LINEWIDTH} / . 
                    qq/lc rgb "$opts{COLORS}->[$j]" pt $opts{POINTTYPE} ps $opts{POINTSIZE} / .
                    ($show_legend ? qq/title "$e->{LABEL}"/ : qq/notitle/));
      if ($max) {
        my $m = 0;
        my $v = $e->{VALUES};
        for (my $j = 1; $j <= $#{$v}; $j++) {
          if ($v->[$j]->[1] > $v->[$m]->[1]) { $m = $j; }
        }
        push(@curves, qq/"< echo '$v->[$m]->[0] $v->[$m]->[1]'" with linespoints lt $lt lw $opts{LINEWIDTH} / . 
                      qq/lc rgb "$opts{MAXCOLOR}" pt $opts{POINTTYPE} ps $opts{POINTSIZE} notitle/);
      }
    }
  }
  return join(", ", @curves);
}


### Miscellanious ###


sub get_entities {
  $max_ncols = 0;
  my $ncols;
  foreach my $infile (@infiles) {
    open FIN, "< $infile" or die "Can't read from '$infile'!";
    my $msg = "Invalid format: '$infile'!";
    my @values;
    while (<FIN>) {
      if (/^#/) { next; }
      chomp;
      s/^\s+//;
      s/\s+$//;
      my @val = split(/\s+/);
      unless ($ncols) { 
        $ncols = scalar(@val); 
        $max_ncols = max($max_ncols, $ncols);
      }
      if (scalar(@val) < 2 || scalar(@val) != $ncols) { die $msg; }
      foreach my $v (@val) {
        unless ($v =~ /^[+-]?\d+(\.\d+)?$/) { die $msg; }
      }
      push(@values,\@val);
    }
    close FIN;
    unless (@values) { die $msg; }

    my @range;
    for (1 .. $ncols) { push(@range, [&INT_MAX, &INT_MIN]); }
    foreach my $val (@values) {
      for my $i (0 .. $ncols - 1) {
        $range[$i]->[0] = min($range[$i]->[0], $val->[$i]);
        $range[$i]->[1] = max($range[$i]->[1], $val->[$i]);
      }
    }

    my %e = (
      FILE   => $infile,
      NCOLS  => $ncols,
      VALUES => \@values,
      RANGE  => \@range,
      LABEL  => $labels[scalar(@entities)] ? $labels[scalar(@entities)] : basename($infile)
    );
    push(@entities, \%e);
  }
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { die "Error executing '$cmd'!\n"; }
  return $?;
}
