#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use My::Utils qw(remove_suffix min max);
use File::Basename qw(basename);
use POSIX;


=pod
=head1 NAME

  sgdviz.pl - Visualizes parameters of an SGD optimization.

=head1 SYNOPSIS

  sgdviz.pl [OPTIONS] --infile INFILE

  OPTIONS:
    -i, --infile INFILE+      List of SGD log files
    -o, --outfile OUTFILE     Output file [default: INFILE.pdf]
    -p, --param PARAM+        Parameters to be plotted [default: ll-train,ll-val]
    -t, --title TITLE         Title of the plots [default: ]
    -l, --labels LABEL+       List of labels to be used instead of directory names
    -L, --labels              Show labels [default: false]
    -g, --legend POS          Position of the legende [default: center right]
    -y, --ymin YMIN           Minimum value of the y1 axis
    -Y, --ymax YMAX           Maximum value of the y1 axis
    -k, --keep                Keep plot files [default: false]
    -h, --help                Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @infiles;
my $outfile;
my %params;
my @aparams;
my @labels;
my $show_labels;
my $title   = "";
my $legpos  = "center right";
my $keep;

my @sgd;
my $ncols   = 5;
my @range;
my @y1range;
my @y2range;
my $single;

my %opts = (
  WIDTH           => 12,
  HEIGHT          => 7,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 4,
  LINEWIDTHS      => 3,
  CURVETYPE       => "lines",
  COLORS          => [
    "#FF0000",  # red
    "#1E90FF",  # dodgerblue
    "#32CD32",  # limegreen
    "#FFA500",  # orange
    "#C71585",  # purple
    "#A0522D",  # sienna
    "#00CED1",  # darkturquoise 
    "#8470ff",  # slate blue
    "#008080",  # teal
    "#FF00FF",  # fuchsia
  ],
  COLORMAX        => "#696969" # dimgray
);
my %col = (
  IT      => 1,
  LLTRAIN => 2,
  PRIOR   => 3,
  LLVAL   => 4,
  NEFF    => 5
);

### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "i|infile=s{1,}"  => \@infiles,
  "o|outfile=s"     => \$outfile,
  "p|params=s{1,}"  => \@aparams,
  "t|title=s"       => \$title,
  "l|label=s{1,}"   => \@labels,
  "L|labels!"       => \$show_labels,
  "g|legend=s"      => \$legpos,
  "y|ymin=f"        => \$y1range[0],
  "Y|ymax=f"        => \$y1range[1],
  "k|keep!"         => \$keep,
  "h|help"          => sub { pod2usage(2); }
) or die pod2usage(1);
unless (@infiles) { pod2usage("Not input file provided!"); }
unless (@aparams) { @aparams = qw(ll-train ll-val); }
foreach my $p (@aparams) { 
  if ($p eq "ll-train") { $params{LLTRAIN} = 1; }
  elsif ($p eq "ll-val") { $params{LLVAL} = 1; }
  elsif ($p eq "prior") { $params{PRIOR} = 1; }
  elsif ($p eq "neff") { $params{NEFF} = 1; }
}
unless (keys(%params)) { pod2usage("No valid parameters provided!"); }
$params{LLTRAIN} = 1;
if (!defined($show_labels) && (scalar(@labels) || scalar(@infiles) > 1)) {
  $show_labels = 1;
}
unless ($outfile) {
  if (scalar(@infiles) > 1) { $outfile = "sgd.pdf"; }
  else { $outfile = &remove_suffix($infiles[0]) . ".pdf"; }
}
$single = scalar(@infiles) == 1;


### Create the plot ###


print "Creating '$outfile'...\n";
for (1 .. $ncols) { push(@range, [&INT_MAX, -&INT_MAX]); }
for my $i (0 .. $#infiles) {
  my $s = &read_sgd($infiles[$i], $labels[$i]);
  for my $i (0 .. $ncols - 1) {
    $range[$i]->[0] = min($range[$i]->[0], $s->{RANGE}->[$i]->[0]);
    $range[$i]->[1] = max($range[$i]->[1], $s->{RANGE}->[$i]->[1]);
  }
  push(@sgd, $s);
}

if (!defined($y1range[0])) {
  $y1range[0] = min($range[$col{LLTRAIN}-1]->[0], $range[$col{LLVAL}-1]->[0]);
}
if (!defined($y1range[1])) {
  $y1range[1] = max($range[$col{LLTRAIN}-1]->[1], $range[$col{LLVAL}-1]->[1]);
}
$y1range[1] += 0.1 * ($y1range[1]-$y1range[0]);
@y2range = @{$range[$col{NEFF}-1]};

foreach my $s (@sgd) {
  my ($min, $max) = (\$s->{RANGE}->[$col{PRIOR}-1]->[0], \$s->{RANGE}->[$col{PRIOR}-1]->[1]);
  for my $i (0 .. $#{$s->{VALUES}}) {
    my $v = \$s->{VALUES}->[$i]->[$col{PRIOR}-1];
    $$v = sprintf("%.4f", $y1range[0] + ($$v - $$min) / ($$max - $$min) * ($y1range[1] - $y1range[0]));
  }
  $$min = $y1range[0];
  $$max = $y1range[1];
  &write_sgd($s); 
}

my $basename = &remove_suffix($outfile);
my $cmd = &plot("$basename.ps");
open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) {
  &exec("rm -f $basename.gp");
  foreach my $s (@sgd) {
    &exec("rm -f @{[&get_datafile($s)]}");
  }
}
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

    set key $legpos
    set title "$title"
    set grid
    set xlabel "Iteration"
    set ylabel "Likelihood"/;
  if ($params{NEFF}) { $cmd .= qq/
    set y2label "N'(p_{cs}) or {\/Symbol t}=1.0"
    set y2range [$y2range[0]:$y2range[1]]
    set y2tics
    set ytics nomirror/;
  }
  if ($single && $params{LLTRAIN}) { $cmd .=
    &add_label($col{LLTRAIN} - 1, $opts{COLORS}->[0]);
  }
  if ($single && $params{LLVAL}) { $cmd .=
    &add_label($col{LLVAL} - 1, $opts{COLORS}->[1]);
  }
  $cmd .= qq/
    set ytics
    set yrange[$y1range[0]:$y1range[1]]
    plot @{[&cmd_curves]}/;

  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub add_label {
  my ($col, $color) = @_;
  my $max = &get_max($sgd[0], $col);
  return sprintf('
    set label "%.3f" at %d,%f center font "%s, %d" tc rgb "%s" front point lt rgb "%s" pt 2 ps 2 offset 0,1', 
    $max->[1], @{$max}, $opts{FONT}, $opts{FONTSIZE}, $color, $color);
}

sub cmd_curves {
  my @curves;
  for my $i (0 .. $#sgd) {
    my $s = $sgd[$i];
    my $datafile = &get_datafile($s);
    if ($params{PRIOR}) {
      push(@curves, sprintf("\"$datafile\" using 1:$col{PRIOR} axes x1y1 with $opts{CURVETYPE} lt 8 lc rgb \"%s\" lw $opts{LINEWIDTHS} %s",
          $opts{COLORS}->[$single ? 3 : $i], $single ? "title \"Prior\"" : "notitle"));
    }
    if ($params{NEFF}) {
      push(@curves, sprintf("\"$datafile\" using 1:$col{NEFF} axes x1y2 with $opts{CURVETYPE} lt 5 lc rgb \"%s\" lw $opts{LINEWIDTHS} %s",
          $opts{COLORS}->[$single ? 2 : $i], $single ? "title \"N'(p_{cs})\"" : "notitle"));
    }
    if ($params{LLVAL}) {
      push(@curves, sprintf("\"$datafile\" using 1:$col{LLVAL} axes x1y1 with $opts{CURVETYPE} lt 2 lc rgb \"%s\" lw $opts{LINEWIDTH} %s",
          $opts{COLORS}->[$single ? 1 : $i], $single ? "title \"LL-Val\"" : "notitle"));
    }
    if ($params{LLTRAIN}) {
      push(@curves, sprintf("\"$datafile\" using 1:$col{LLTRAIN} axes x1y1 with $opts{CURVETYPE} lt 1 lc rgb \"%s\" lw $opts{LINEWIDTH} %s",
          $opts{COLORS}->[$single ? 0 : $i], $single ? "title \"LL-Train\"" : "title \"$s->{LABEL}\""));
    }
  }
  return join(", ", @curves);
}


### Miscallenous ###


sub read_sgd {
  my ($file, $label) = @_;
  my @values;
  my @range;
  for (1 .. $ncols) { push(@range, [&INT_MAX, -&INT_MAX]); }

  #Epoch                  LL-Train    Prior            LL-Val     Neff
  #1     [==============]   2.3253  +1.0918  -0.0325   2.1410  14.5634
  open FIN, "< $file" or die "Can't read from '$file'!";
  while (<FIN>) {
    if (my @val = /^(\d+)\s+\[=+\]\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) { 
      $val[2] *= -1;
      for my $i (0 .. $ncols - 1) {
        $range[$i]->[0] = min($range[$i]->[0], $val[$i]);
        $range[$i]->[1] = max($range[$i]->[1], $val[$i]);
      }
      push(@values, \@val); 
    } elsif (scalar(@values)) { last; }
  }
  close FIN;
  unless (scalar(@values)) { die "Format of '$file' invalid!"; }

  my %s;
  $s{FILE} = $file;
  if ($label) {
    $s{LABEL} = $label;
  } else {
    $s{LABEL} = basename($file);
    $s{LABEL} =~ s/_/ /g;
  }
  $s{VALUES} = [sort({ $a->[0] <=> $b->[0] } @values)];
  $s{RANGE} = \@range;
  return \%s;
}


sub write_sgd {
  my ($s) = @_;
  my $datafile = &get_datafile($s);
  open FOUT, "> $datafile" or die "Can't write to '$datafile'";
  printf(FOUT "#%9s%10s%10s%10s%10s\n", "Epoch", "LL-Train", "Prior", "LL-Val", "Neff");
  foreach my $v (@{$s->{VALUES}}) {
    foreach my $w (@{$v}) { printf(FOUT "%10s", $w); }
    print(FOUT "\n");
  }
  close FOUT;
}

sub get_datafile {
  my ($s) = @_;
  return "$s->{FILE}.sgd";
}

sub get_max {
  my ($s, $col) = @_;
  my $m = 0;
  for my $i (1 .. $#{$s->{VALUES}}) {
    if ($s->{VALUES}->[$i]->[$col] > $s->{VALUES}->[$m]->[$col]) { $m = $i; }
  }
  return [$s->{VALUES}->[$m]->[0], $s->{VALUES}->[$m]->[$col]];
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { print STDERR "Error executing '$cmd'!\n"; }
  return $?;
}
