#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename);


=pod

=head1 NAME

  benchviz.pl - Visualizes CS-BLAST benchmark results.

=head1 SYNOPSIS

  benchviz.pl [OPTIONS] --dir BENCHDIR+

  OPTIONS:
 
    -d, --dir BENCHDIR+     List of directories containing benchmark result files
    -o, --out OUTBASE       Output basename [default: ./]
    -p, --plot PLOT+        List of plots to be created [default: tpfp wtpfp rocx]
    -t, --title TITLE       Title of the plots [default: ]
    -l, --label LABEL+      List of labels to be used instead of directory names
    -s, --sort SORT+        List of directories defining the plot order
    -k, --keep              Keep plot files [default: false]
    -h, --help              Show this help message

=head1 AUTHOR
  
  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###

my @dirs;
my $outbase = "./";
my @plots;
my $title;
my @labels;
my @order;
my $keep;

my @entities;

my %opts = (
  WIDTH           => 12,
  HEIGHT          => 7,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 4,
  LINEWIDTHE      => 3,
  COLORS          => [
    "#696969",  # dimgray
    "#FF0000",  # red
    "#1E90FF",  # dodgerblue
    "#32CD32",  # limegreen
    "#C71585",  # purple
    "#FFA500",  # orange
    "#FF00FF",  # fuchsia
    "#8470ff",  # slate blue
    "#696969",  # 20 % error rate
  ]
);


### Initialization ###


GetOptions(
  "d|dir=s{1,}"   => \@dirs,
  "o|out=s"       => \$outbase,
  "p|plot=s{1,}"  => \@plots,
  "t|title=s"     => \$title,
  "l|label=s{1,}" => \@labels,
  "s|sort=s{1,}"  => \@order,
  "k|keep"        => \$keep,
  "h|help"        => sub { pod2usage(2); }
) or pod2usage(1);
unless (@dirs) { pod2usage("No benchmark directory provided!"); }
unless (@plots) { @plots = qw/tpfp wtpfp rocx/; }
&get_entities;
unless (@entities) { pod2usage("No plot data found in the specified benchmark directories!"); }


### Create plots ###


print "Creating plots...\n";
foreach my $plot (@plots) {
  my $b = 0;
  foreach my $e (@entities) {
    if ($e->{PLOTS}->{$plot}) {
      $b = 1;
      last;
    }
  }
  unless ($b) { next; }
  my $plotbase = sprintf("%s%s%s", 
    $outbase, $outbase =~ /\/$/ ? "" : "_", &get_name($plot));
  print "  - $plotbase.pdf\n";
  open FOUT, "> $plotbase.gp" or die "Can't write to '$plotbase.gp'!";
  print FOUT &plot($plot, "$plotbase.ps");
  close FOUT;
  if (&exec("gnuplot $plotbase.gp")) { next; }
  unless ($keep) { &exec("rm $plotbase.gp"); }
  if (&exec("ps2pdf $plotbase.ps $plotbase.pdf")) { next; }
  &exec("rm $plotbase.ps");
  &exec("pdfcrop $plotbase.pdf $plotbase.pdf");
}
print "Done!\n";


### Plotting ###


sub plot {
  my ($plot, $outfile) = @_;
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$outfile"

    set style line 1 linetype 2 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[0]"
    set style line 2 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[1]"
    set style line 3 linetype 3 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[1]"
    set style line 4 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[2]"
    set style line 5 linetype 3 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[2]"
    set style line 6 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[5]"
    set style line 7 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[6]"
    set style line 8 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[7]"
    set style line 9 linetype 0 linewidth $opts{LINEWIDTHE} linecolor rgb "$opts{COLORS}->[8]"/;

  if ($plot eq "tpfp") { $cmd .= qq/
    set key top left
    set xrange [1:4000]
    set yrange [0:10000]
    set log x
    set grid
    set xlabel "FP"
    set ylabel "TP"
    set label "1%" at 35,5000 textcolor ls 9
    set label "10%" at 350,5000 textcolor ls 9
    set label "20%" at 800,5000 textcolor ls 9/;
  } elsif ($plot eq "wtpfp") { $cmd .= qq/
    set key top right
    set xrange [1:600]
    set yrange [0:600]
    set log x
    set grid
    set xlabel "weighted FP"
    set ylabel "weighted TP" 
    set label "1%" at 3,390 textcolor ls 9
    set label "10%" at 30,390 textcolor ls 9
    set label "20%" at 70,390 textcolor ls 9/;
  } elsif ($plot eq "rocx") { $cmd .= qq/
    set key top right
    set xrange [0:1.0]
    set yrange [0:0.5]
    set grid
    set xlabel "ROC5"
    set ylabel "Fraction of queries"/;
  } 
  $cmd .= qq/
    plot /;
  if ($plot eq "tpfp" || $plot eq "wtpfp") { 
    $cmd .= qq/99*x notitle with lines ls 9, 9*x notitle with lines ls 9, 4*x notitle with lines ls 9,/;
  }
  $cmd .= &cmd_curves($plot);
  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub cmd_curves {
  my ($plot) = @_;
  my @curves;
  for my $i (0 .. $#entities) {
    my $e = $entities[$i];
    if ($e->{PLOTS}->{$plot}) {
      push(@curves, sprintf('"%s" title "%s" with lines ls %d', 
          &get_datafile($plot, $e->{DIR}), $e->{LABEL}, $i + 1));
    }
  }
  return join(", ", @curves);
}


### Miscellanious ###


sub get_entities {
  my @unsorted;
  for my $i (0 .. $#dirs) {
    my %e;
    $e{DIR} = $dirs[$i];
    $e{NAME} = basename($dirs[$i]);
    $e{ROCX} = &get_rocx($dirs[$i]);    
    if ($labels[$i]) {
      $e{LABEL} =  $labels[$i]; 
    } else {
      $e{LABEL} = basename($e{DIR}); 
      $e{LABEL} =~ s/_/ /g;
    }
    if ($e{ROCX}) { $e{LABEL} .= sprintf(": %.3f", $e{ROCX}); }
    my $has_file = 0;
    foreach my $p (@plots) {
      my $h =  -f sprintf("%s/%s", $dirs[$i], &get_datafile($p)) ? 1 : 0;
      if ($h) { $has_file = 1; }
      $e{PLOTS}->{$p} = $h;
    }
    if ($has_file) { push(@unsorted, \%e); }
  }
  @entities = ();
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

sub get_datafile {
  my ($plot, $dir) = @_;
  $plot .= ".dat";
  if ($dir) { $plot = sprintf("%s/%s", $dir, $plot); }
  return $plot;
}

sub get_name {
  my ($plot) = @_;
  return $plot;
}

sub get_rocx {
  my ($dir) = @_;
  my $filename = &get_datafile("rocx", $dir);
  unless (-f $filename) { return undef; }
  my $rocx = 0;
  my $x_last = 0;
  my $y_last = 1;
  open(FIN, "< $filename") or die "$filename: $!";
  while (<FIN>) {
    if (my ($x, $y) = /^([^\#]\S*)\s+(\S+)\s*$/) {
      $rocx += ($x - $x_last) * 0.5 * ($y + $y_last);
      ($x_last, $y_last) = ($x, $y);
    }
  }
  return $rocx;
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd &> /dev/null")) { print STDERR "Error executing '$cmd'!\n"; }
  return $?;
}
