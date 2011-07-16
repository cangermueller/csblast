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

  rocxviz.pl - Visualizes the ROCX score of the given directories against a specified parameter.

=head1 SYNOPSIS

  rocxviz.pl [OPTIONS] --dir DIR+ --param PARAM [ROCX.PL OPTIONS]

  OPTIONS:
    -d, --dir DIR           Input directories with rocx curve.
    -p, --param PARAM       Plot parameter which must be part of the directory name (.*_PARAM.*).
    -o, --outfile OUTFILE   Output plot file [default: rocx.pdf]
    -t, --title TITLE       Title of the plot [default: OUTFILE]
    -l, --label LABEL+      List of labels to be used instead of directory names
    -L, --labels            Show labels [default: false]
    -g, --legend POS        Position of the legende [default: center right]
    -s, --sort SORT+        List of directories defining the plot order
    -x, --xlabel XLABEL     Label of this x-axis
    -y, --ylabel YLABEL     Label of this y-axis [default: ROC5]
    -u, --xmin YMIN         Minimum value of the x axis
    -U, --xmax YMAX         Maximum value of the x axis
    -v, --ymin YMIN         Minimum value of the y axis
    -V, --ymax YMAX         Maximum value of the y axis
    -S, --xtic ]0;inf[      Step size of xtics
    -m, --max               Highlight maximum [default: true]
    -k, --keep              Keep data and plot files [default: false]
    -D, --data-only         Only write data files to output directory
    -h, --help              Show this help message


=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my @dirs;
my $param;
my $outfile = "./rocx.pdf";
my $title   = "";
my $xlabel;
my $ylabel  = "ROC_5";
my $show_labels;
my @xrange;
my @yrange;
my @labels;
my @order;
my $max     = 1;
my $keep;
my $data_only;
my $rocx_opts;

my @entities;
my $workdir;
my $basename;
my @range;
my $legpos  = "top right";
my @xtics;

my %opts = (
  WIDTH           => 10,
  HEIGHT          => 10,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 4,
  POINTSIZE       => 2,
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
  "d|dir=s{1,}"   => \@dirs,
  "p|param=s"     => \$param,
  "o|outfile=s"   => \$outfile,
  "t|title=s"     => \$title,
  "x|xlabel=s"    => \$xlabel,
  "y|ylabel=s"    => \$ylabel,
  "l|label=s{1,}" => \@labels,
  "L|labels!"     => \$show_labels,
  "g|legend=s"    => \$legpos,
  "u|xmin=f"      => \$xrange[0],
  "U|xmax=f"      => \$xrange[1],
  "v|ymin=f"      => \$yrange[0],
  "V|ymax=f"      => \$yrange[1],
  "S|xtics=f"     => \$xtics[1],
  "s|sort=s{1,}"  => \@order,
  "m|max!"        => \$max,
  "k|keep!"       => \$keep,
  "D|data-only!"  => \$data_only,
  "h|help"        => sub { pod2usage(2); }
) or die pos2usage(1);
$rocx_opts = join(" ", @ARGV);
unless (@dirs) { pod2usage("No input directories provided!"); }
unless ($param) { pod2usage("No parameter provided!"); }
$workdir = dirname($outfile);
$basename = remove_suffix($outfile);
unless ($xlabel) { $xlabel = $param};
if (@labels) { $show_labels = 1; }


### Create the plot ###


&get_entities;
print "Creating '$outfile'...\n";
for (1 .. 2) { push(@range, [&INT_MAX, -&INT_MAX]); }
foreach my $e (@entities) { 
  for my $i (0 .. 1) {
    $range[$i]->[0] = min($range[$i]->[0], $e->{RANGE}->[$i]->[0]);
    $range[$i]->[1] = max($range[$i]->[1], $e->{RANGE}->[$i]->[1]);
  }
  &write_datafile($e); 
}
if ($data_only) {
  print "Data files written to '$workdir'!\n";
  exit 0;
}
unless (defined($xrange[0])) { $xrange[0] = $range[0]->[0]; }
unless (defined($xrange[1])) { $xrange[1] = $range[0]->[1]; }
unless (defined($yrange[0])) { $yrange[0] = $range[1]->[0]; }
unless (defined($yrange[1])) { $yrange[1] = $range[1]->[1] + 0.1 * ($range[1]->[1] - $range[1]->[0]); }
$xtics[0] = $xrange[0];
unless ($xtics[1]) { $xtics[1] = $entities[0]->{VALUES}->[1]->[0] - $entities[0]->{VALUES}->[0]->[0]; }
my $cmd = &plot("$basename.ps");
open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) {
  &exec("rm -f $basename.gp");
  foreach my $e (@entities) {
    &exec("rm -f @{[&get_datafile($e)]}");
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

    set style line 1 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[0]"
    set style line 2 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[1]"
    set style line 3 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[2]"
    set style line 4 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[3]"
    set style line 5 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[4]"
    set style line 6 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[5]"
    set style line 7 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[6]"
    set style line 8 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[7]"
    set style line 9 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{MAXCOLOR}"

    set key $legpos
    set border 15 back
    set title "$title"
    set grid
    set xrange[$xrange[0]:$xrange[1]]
    set xtics $xtics[0], $xtics[1]
    set xlabel "$xlabel"
    set yrange[$yrange[0]:$yrange[1]]
    set ylabel "$ylabel"

    plot @{[&cmd_curves]}
    /;
  $cmd =~ s/^\s+//mg;
  # print "$cmd\n"; exit 0;
  return "$cmd\n";
}

sub cmd_curves {
  my @curves;
  for my $i (0 .. $#entities) {
    my $e = $entities[$i];
    my $curve = sprintf('"%s" with linespoints ls %d %s',
      &get_datafile($e), $i + 1, $show_labels ? "title \"$e->{LABEL}\"" : "notitle");
    if ($max) {
      my $m = 0;
      my $v = $e->{VALUES};
      for (my $j = 1; $j < $#{$v}; $j++) {
        if ($v->[$j]->[1] > $v->[$m]->[1]) { $m = $j; }
      }
      $curve .= ", \"< echo '$v->[$m]->[0] $v->[$m]->[1]'\" with points ls 9 notitle";
    }
    push(@curves, $curve);
  }
  return join(", ", @curves);
}


### Miscellanious ###


sub get_datafile {
  my ($e) = @_;
  return sprintf("%s/%s.dat", $workdir, $e->{NAME});
}

sub write_datafile {
  my ($e) = @_;
  my $filename = &get_datafile($e);
  my @values = sort({$a->[0] <=> $b->[0]} @{$e->{VALUES}});
  open FOUT, "> $filename" or die "Can't write to '$filename'!";
  printf(FOUT "#%9s%10s\n", $param, "ROCX");
  foreach my $v (@values) {
    printf(FOUT "%10s%10s\n", $v->[0], $v->[1]);
  }
  close FOUT;
}

sub get_entities {

  my %order;
  my $i = 0;
  foreach my $dir (@dirs) {
    my $id = &get_id($dir);
    if ($id && !exists($order{$id})) {
      $order{$id} = $i++;
    }
  }
  unless ($i) { print STDERR "No ROCX data for the specified parameter!"; }

  $_ = sprintf("rocx.pl %s -d %s", $rocx_opts, join(" ", @dirs));
  $_ = `$_`;
  if ($?) { die "Error calling rocx.pl!"; }
  my %h;
  foreach (split(/\n/)) {
    if (my ($rocx, $dir) = /^\s*\d+\s+(\S+)\s+(\S+)\s*$/) {
      $dir = basename($dir);
      if (my ($p) = $dir =~ /_$param(\d[^_]*)/) {
        $dir =~ s/_$param[^_]+/_${param}X/;
        push(@{$h{$dir}}, [$p, $rocx]);
      }
    }
  }
  unless (scalar(keys(%h))) {
    print STDERR "No ROCX data for the specified parameter!";
    exit 1;
  }
  my @unsorted;
  while (my ($key, $value) = each(%h)) {
    my %e;
    $e{NAME} = $key;
    $e{VALUES} = [sort({ $a->[0] <=> $b->[0] } @{$value})];
    if ($labels[$order{$key}]) {
      $e{LABEL} = $labels[$order{$key}];
    } else {
      $e{LABEL} = $key;
      $e{LABEL} =~ s/_/ /g;
    }
    my @range;
    for (1 .. 2) { push(@range, [&INT_MAX, -&INT_MAX]); }
    foreach my $v (@{$e{VALUES}}) {
      for my $i (0 .. 1) {
        $range[$i]->[0] = min($range[$i]->[0], $v->[$i]);
        $range[$i]->[1] = max($range[$i]->[1], $v->[$i]);
      }
    }
    $e{RANGE} = \@range;
    push(@unsorted, \%e);
  }
  @unsorted = sort({ $order{$a->{NAME}} <=> $order{$b->{NAME}} } @unsorted);
  foreach my $o (@order) {
    for my $i (0 .. $#unsorted) {
      my $e = $unsorted[$i];
      if ($e->{NAME} eq $o) { 
        push(@entities, $e); 
        splice(@unsorted, $i, 1);
        last;
      }
    }
  }
  push(@entities, @unsorted);
}

sub get_id {
  my ($dir) = @_;
  $dir = basename($dir);
  if (my ($p) = $dir =~ /_$param(\d[^_]*)/) {
    $dir =~ s/_$param[^_]+/_${param}X/;
    return $dir;
  }
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { die "Error executing '$cmd'!\n"; }
  return $?;
}
