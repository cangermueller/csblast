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

  trainsetviz.pl - Visualizes the diversity of a training set

=head1 SYNOPSIS

  trainsetviz.pl [OPTIONS] -i INFILE

  OPTIONS:
    -i, --infile INFILE       Input training set
    -o, --outfile OUTFILE     Output file [default: INFILE.pdf]
    -t, --title TITLE         Title of the plot [default: ]
    -k, --keep                Keep data and plot files [default: false]
    -u, --xmin XMIN           Minimum value of the x axis
    -U, --xmax XMAX           Maximum value of the x axis
    -v, --ymin YMIN           Minimum value of the y axis
    -V, --ymax YMAX           Maximum value of the y axis
    -b, --bin [0;inf[         Bin size [default: 0.5]
    -N, --num                 Number of training samples [default = 10000]
    -h, --help                Show this help message


=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###

my $infile;
my $outfile;
my $keep;
my $title    = "";
my @xrange;
my @yrange;
my $bin_size = 0.5;
my $neff_max = 20.0;
my $nsamples = 10000;

my $outbase;
my @entities;

my %opts = (
  WIDTH           => 10,
  HEIGHT          => 10,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
  LINEWIDTH       => 2,
  LINETYPE        => 1,
  POINTSIZE       => 2,
  POINTTYPE       => 7,
  CURVETYPE       => "boxes",
  COLORS          => [
    "#32CD32",  # limegreen
    "#1E90FF",  # dodgerblue
    "#FFA500",  # orange
    "#FF0000",  # red
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
  "i|infile=s"  => \$infile,
  "o|outfile=s" => \$outfile,
  "t|title=s"   => \$title,
  "k|keep!"     => \$keep,
  "u|xmin=f"    => \$xrange[0],
  "U|xmax=f"    => \$xrange[1],
  "v|ymin=f"    => \$yrange[0],
  "V|ymax=f"    => \$yrange[1],
  "b|bin=f"     => \$bin_size,
  "N|num=i"     => \$nsamples,
  "h|help"      => sub { pod2usage(2); }
) or die pos2usage(1);
unless ($infile) { pod2usage("No input files provided!"); }
unless ($outfile) { $outfile = basename($infile) . ".pdf"; }
$outbase = remove_suffix($outfile);


### Create the plot ###


print "Calculating the Neff in $infile...\n";
&exec("cstrainset_neff -i $infile -x $outbase.x -y $outbase.y -N $nsamples");
print "Creating $outfile...\n";
my @xr = (&INT_MAX, -&INT_MAX);
my @yr = (&INT_MAX, -&INT_MAX);
foreach my $t (qw(x y)) {
  my $e = &get_entity("$outbase.$t");
  $xr[0] = min($xr[0], $e->{BINRANGE}->[0]);
  $xr[1] = max($xr[1], $e->{BINRANGE}->[1]);
  $yr[1] = max($yr[1], &get_max($e->{DENS}));
  $e->{FILE} = "$outbase.$t.bin";
  &write_entity($e);
  push(@entities, $e);
}
$entities[0]->{LABEL} = sprintf(" N(c): %5.2f", &get_mean($entities[0]->{VALUES}));
$entities[1]->{LABEL} = sprintf("N(c'): %5.2f", &get_mean($entities[1]->{VALUES}));
unless ($xrange[0]) { $xrange[0] = 0.5; }
unless ($xrange[1]) { $xrange[1] = ($xr[1] + 2) * $bin_size; }
unless ($yrange[0]) { $yrange[0] = 0.0; }
unless ($yrange[1]) { $yrange[1] = $yr[1] + 0.05; }

my $cmd = &plot("$outbase.ps");
open FOUT, "> $outbase.gp" or die "Can't write to '$outbase.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $outbase.gp");
unless ($keep) { 
  &exec("rm -f $outbase.gp $outbase.x $outbase.y"); 
  foreach my $e (@entities) {
    system("rm -f $e->{FILE}");
  }
}
&exec("ps2pdf $outbase.ps $outbase.pdf");
&exec("rm -f $outbase.ps");
&exec("pdfcrop $outbase.pdf $outfile");
print "Done!\n";


### Plotting ###


sub plot {
  my ($outfile) = @_;
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$outfile"

    set title "$title"
    set key top right
    set border 15 back
    set grid

    set xlabel "Neff"
    set xrange [$xrange[0]:$xrange[1]]
    set xtics $bin_size out nomirror

    set ylabel "Frequency"
    set yrange[$yrange[0]:$yrange[1]]
    set ytics 0.1 out nomirror

    set style line 1 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[1]"
    set style line 2 lt 1 lw $opts{LINEWIDTH} pt $opts{POINTTYPE} ps $opts{POINTSIZE} lc rgb "$opts{COLORS}->[3]"
    set style line 3 lt 1 lw 1 lc rgb "black"

    set style fill solid 0.5 border

    plot /;
  $cmd .= sprintf("%s, %s", &cmd_curve($entities[0], 1), &cmd_curve($entities[1], 2));
  $cmd =~ s/^\s+//mg;
  # print "$cmd\n"; exit 0;
  return "$cmd\n";
}

sub cmd_curve {
  my ($e, $i) = @_;
  return qq/"$e->{FILE}" with $opts{CURVETYPE} ls $i title "$e->{LABEL}"/; 
}


### Miscellanious ###


sub get_entity {
  my ($file) = @_;

  my @values;
  my @freq;
  my @dens;
  my @bin_range = (&INT_MAX, -&INT_MAX);
  my $count = 0;

  for (0 .. &get_bin($neff_max)) { push(@freq, 0); }
  my $msg = "Invalid format: '$file'!";
  open FIN, "< $file" or die "Can't read from '$file'!";
  while (<FIN>) {
    if (/^#/) { next; }
    unless (/^\s*(\d+(\.\d+)?)/ && $1 <= $neff_max) { die $msg; }
    push(@values, $1);
    my $b = &get_bin($1);
    $bin_range[0] = min($bin_range[0], $b);
    $bin_range[1] = max($bin_range[1], $b);
    $freq[$b]++;
    $count++;
  }
  close FIN;
  unless (@values) { die $msg; }

  @dens = @freq;
  for my $i ($bin_range[0] .. $bin_range[1]) { $dens[$i] /= $count; }
  return {
    VALUES   => \@values,
    FREQ     => \@freq,
    DENS     => \@dens,
    BINRANGE => \@bin_range,
    COUNT    => $count,
  }
}

sub write_entity {
  my ($e) = @_;

  open FOUT, "> $e->{FILE}" or die "Can't write to '$e->{FILE}'!";
  my $x = $e->{BINRANGE}->[0] * $bin_size + 0.5 * $bin_size;
  for my $i ($e->{BINRANGE}->[0] .. $e->{BINRANGE}->[1]) {
    printf(FOUT "%6.3f    %.5f\n", $x, $e->{DENS}->[$i]);
    $x += $bin_size;
  }
  close FOUT;
}

sub get_bin {
  my ($value) = @_;
  return int($value / $bin_size);
}

sub exec {
  my $cmd = $_[0];
  if (system("$cmd &> /dev/null")) { die "Error executing '$cmd'!\n"; }
  return $?;
}

sub get_max {
  my ($a) = @_;
  my $max = $a->[0];
  for my $i (1 .. $#{$a}) { $max = max($max, $a->[$i]); }
  return $max;
}

sub get_mean {
  my ($a) = @_;
  my $sum = 0;
  foreach (@{$a}) { $sum += $_; }
  return $sum / scalar(@{$a});
}
