#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use POSIX;
use My::Utils qw(remove_suffix trim min max);


=pod
=head1 NAME

  csoptviz.pl - Visualizes the optimization table of CSOPT

=head1 SYNOPSIS

  csoptviz.pl [OPTIONS] --infile INFILE

  OPTIONS:
    -i, --infile INFILE       CSOPT optimization table
    -o, --outfile OUTFILE     Output file [default: INFILE.pdf]
    -p, --param PARAM+        Parameters to be plotted [default: lltrain]
    -t, --title TITLE         Title of the plots [default: PARAM]
    -x, --x-range RANGE       GNU plot x-range specification
    -y, --y-range RANGE+      GNU plot y-range specification
    -k, --keep                Keep plot files [default: false]
    -h, --help                Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $infile;
my $outfile;
my $basename;
my @plot_params;
my $xrange;
my @yrange;
my $title;
my $keep;

my @entries;
my %range;
my @params;
my $nparams;
my $show_ll;

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


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "i|infile=s"      => \$infile,
  "o|outfile=s"     => \$outfile,
  "p|params=s{1,}"  => \@plot_params,
  "t|title=s"       => \$title,
  "x|x-range=s"     => \$xrange,
  "y|y-range=s{1,}" => \@yrange,
  "k|keep!"         => \$keep,
  "h|help"          => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($infile) { pod2usage("Not input file provided!"); }
unless (@plot_params) { @plot_params = qw/lltrain/; }
unless ($outfile) { $outfile = "$infile.pdf"; }
$basename = &remove_suffix($outfile);


### Parse the optimization table ###


print "Reading optimization table ...\n";
&parse_table;


### Create the plot ###


print "Creating '$outfile' ... \n";
my $cmd;
if ($plot_params[0] eq "lltrain") { $cmd= &plot_scatter; }
else { die "Invalid plot parameter!"; }
open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) {
  foreach my $s (qw/dat gp/) {
    &exec("rm -f $basename.$s");
  }
}
&exec("ps2pdf $basename.ps $basename.pdf");
&exec("rm -f $basename.ps");
&exec("pdfcrop $basename.pdf $outfile");
print "Done!\n";


### Plotting ###


sub plot_scatter {
  unless ($show_ll) { 
    print "Optimization table does not list the likelihood!\n";
    exit 1;
  }

  open FOUT, "> $basename.dat" or die "Can't write to '$basename.dat'!";
  printf(FOUT "#%7s %8s %8s\n", "lltrain", "llval", "rocx");
  foreach my $e (@entries) {
    printf(FOUT "%8.4f %8.4f %8.4f\n", $e->{LLTRAIN}, $e->{LLVAL}, $e->{ROCX});

  }
  
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$basename.ps"
    #set multiplot layout 1, 2
    set key off
    @{[&plot_scatter_cmd(0 , "llval")]}/;
#    @{[&plot_scatter_cmd(1 , "rocx")]}/;
  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub plot_scatter_cmd {
  my ($i, $label) = @_;
  my $j = $i + 2;
  my $t = $title ? $title : "lltrain vs. $label";
  my $col = $opts{COLORS}->[$i + 1];
  my $yr;
  if ($yrange[$i]) { 
    $yr = $yrange[$i]; }
  else { 
    $yr = sprintf("[%g:%g]", $i ? @{$range{ROCX}} : @{$range{LLVAL}});
  }

  my $cmd = qq/
    set title "$t"
    set xlabel "lltrain"
    set ylabel "$label"
    set yrange $yr/;
  if ($xrange) { $cmd .= qq/
    set xrange $xrange/;
  }
  $cmd .= qq/
    f(x) = a*x+b
    fit f(x) "$basename.dat" using 1:$j via a, b
    plot "$basename.dat" using 1:$j notitle with points ps 3 pt 1 lw 2 lc rgb "$col", / .
      qq/f(x) notitle with lines lw 2 lt 1 lc rgb "red", / . 
  $cmd =~ s/^\s+//mg;
  return $cmd;
}


### Utils ###


sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { print STDERR "Error executing '$cmd'!\n"; exit 1; }
  return $?;
}

sub invalid_format {
  die "Invalid format: '$infile' line $.!\n";
}

sub parse_table {
  open FIN, "< $infile" or die "Can't read from '$infile'!";
  # parse header
  $_ = <FIN>;
  /^#\s+r\s+optimizing\s+(\S+(?:\s+\S+)*?)\s+(lltrain\s+llval)?\s+ROC5\s+\+\/\-%\s*$/ or &invalid_format;
  @params = split(/\s+/, $1);
  $nparams = scalar(@params);
  $show_ll = defined($2);
  <FIN>;

  # parse content
  my $k = 0;
  my $r = 0;
  my $opt = "";
  my @params_cur;
  while (<FIN>) {
    last if /^-/;
    my %e;
    /^(.+?)\s+(\S+\s+\S+)?\s+(\S+)\s+(\S+)$/ or &invalid_format;
    $e{ROCX} = $3;
    $e{GAIN} = $4;
    if ($show_ll) {
      $2 or &invalid_format;
      @_ = split(/\s+/, $2);
      $e{LLTRAIN} = $_[0];
      $e{LLVAL} = $_[1];
    }
    my @f = split(/\s+/, trim($1));
    unless (scalar(@f)) {
      unless (@params_cur) { &invalid_format; }
      $e{K} = $k;
      $e{R} = $r;
      $e{OPT} = $opt;
      $e{PARAMS} = [@params_cur];
    } else {
      if (scalar(@f) < $nparams) { &invalid_format; }
      @params_cur = @f[scalar(@f) - $nparams ... $#f];
      $e{PARAMS} = [@params_cur];
      for (1 .. $nparams) { pop(@f); }
      if (scalar(@f) == 3) {
        ($k, $r, $opt) = @f; 
      } elsif (scalar(@f) != 0) { &invalid_format; }
      $e{K} = $k;
      $e{R} = $r;
      $e{OPT} = $opt;
    }
    push(@entries, \%e);
  }
  close FIN;
  scalar(@entries) or &invalid_format;

  # filter entries with invalid rocx score
  my @filtered;
  foreach my $e (@entries) {
    if ($e->{ROCX} =~ /^[+-]?\d+\.\d+$/ && $e->{ROCX} > 0) { push(@filtered, $e); }
  }
  @entries = @filtered;

  # get the range of each field
  $range{LLTRAIN} = [&DBL_MAX, -&DBL_MAX];
  $range{LLVAL} = [&DBL_MAX, -&DBL_MAX];
  $range{ROCX} = [&DBL_MAX, -&DBL_MAX];
  $range{GAIN} = [&DBL_MAX, -&DBL_MAX];
  $range{PARAMS} = [];
  for (1 .. $nparams) {
    push(@{$range{PARAMS}}, [&DBL_MAX, -&DBL_MAX]);
  }
  foreach my $e (@entries) {
    $range{LLTRAIN}->[0] = min($range{LLTRAIN}->[0], $e->{LLTRAIN});
    $range{LLTRAIN}->[1] = max($range{LLTRAIN}->[1], $e->{LLTRAIN});
    $range{LLVAL}->[0] = min($range{LLVAL}->[0], $e->{LLVAL});
    $range{LLVAL}->[1] = max($range{LLVAL}->[1], $e->{LLVAL});
    $range{ROCX}->[0] = min($range{ROCX}->[0], $e->{ROCX});
    $range{ROCX}->[1] = max($range{ROCX}->[1], $e->{ROCX});
    $range{GAIN}->[0] = min($range{GAIN}->[0], $e->{GAIN});
    $range{GAIN}->[1] = max($range{GAIN}->[1], $e->{GAIN});
    for my $i (0 .. $nparams - 1) {
      $range{PARAMS}->[$i]->[0] = min($range{PARAMS}->[$i]->[0], $e->{PARAMS}->[$i]);
      $range{PARAMS}->[$i]->[1] = max($range{PARAMS}->[$i]->[1], $e->{PARAMS}->[$i]);
    }
  }
}
