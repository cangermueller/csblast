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
    -p, --plot X Y            Plot X agains Y [default: lltrain llval]
    -t, --title TITLE         Title of the plots [default: PARAM]
    -x, --x-range RANGE       GNU plot x-range specification
    -y, --y-range RANGE       GNU plot y-range specification
    -r, --rounds [1,inf[      Rounds to be plotted
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
my @plot_xy;
my $xrange;
my $yrange;
my @rounds;
my $title;
my $keep;

my @groups;
my %params;
my $nparams;
my $has_ll;

my %opts = (
  WIDTH           => 12,
  HEIGHT          => 7,
  FONT            => "bold",
  FONTSIZE        => 14,
  LINEWIDTH_SCALE => 2,
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
);


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "i|infile=s"      => \$infile,
  "o|outfile=s"     => \$outfile,
  "p|plot=s{2}"     => \@plot_xy,
  "t|title=s"       => \$title,
  "x|x-range=s"     => \$xrange,
  "y|y-range=s"     => \$yrange,
  "r|rounds=i{1,}"  => \@rounds,
  "k|keep!"         => \$keep,
  "h|help"          => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($infile) { pod2usage("Not input file provided!"); }
unless (@plot_xy) { @plot_xy = qw/lltrain llval/; }
unless ($outfile) { $outfile = "$infile.pdf"; }
$basename = &remove_suffix($outfile);


### Parse the optimization table ###


print "Reading optimization table ...\n";
&parse_table;
my %h;
foreach my $g (@groups) { $h{$g->[0]->{"r"}} = 1; }
delete $h{0};
unless (@rounds) {
  @rounds = keys(%h);
} else {
  my @rf;
  foreach my $r (@rounds) {
    if (defined($h{$r})) { push(@rf, $r); }
  }
  @rounds = @rf;
}


### Create the plot ###


print "Creating '$outfile' ... \n";

my $cmd;
if ($plot_xy[0] eq "lltrain" &&  $plot_xy[1] eq "llval") {
  $cmd = &plot_lltrain_llval;
} elsif (defined($params{$plot_xy[0]}) && defined($groups[0]->[0]->{$plot_xy[1]})) {
  $cmd = &plot_xy;
} else {
  printf(STDERR "'%s' can't be plotted against '%s'!\n", @plot_xy);
  exit 1;
}

open FOUT, "> $basename.gp" or die "Can't write to '$basename.gp'!";
print FOUT $cmd;
close FOUT;
&exec("gnuplot $basename.gp");
unless ($keep) {
  foreach my $s (qw/*dat gp/) {
    &exec("rm -f $basename.$s");
  }
}
&exec("ps2pdf $basename.ps $basename.pdf");
unless ($keep) {
  &exec("rm -f $basename.ps");
}
&exec("pdfcrop $basename.pdf $outfile");
print "Done!\n";


### Plotting ###


sub plot_lltrain_llval {
  unless ($has_ll) { 
    print "Optimization table does not list the likelihood!\n";
    exit 1;
  }

  # write data file using the upper x% of entries
  my @entries;
  foreach my $g (@groups) {
    foreach my $e (@{$g}) { push(@entries, $e); }
  }
  @entries = sort({ $a->{"lltrain"} <=> $b->{"lltrain"} } @entries);
  @entries = @entries[0.4 * scalar(@entries) .. $#entries];
  my %range = &get_range(\@entries);
  open FOUT, "> $basename.dat" or die "Can't write to '$basename.dat'!";
  printf(FOUT "#%7s %8s %8s\n", "lltrain", "llval", "score");
  foreach my $e (@entries) {
    printf(FOUT "%8.4f %8.4f %8.4f\n", $e->{"lltrain"}, $e->{"llval"}, $e->{"score"});
  }
  my $xr = $xrange ? $xrange : sprintf("[%g:%g]", @{$range{"lltrain"}});
  my $yr = $yrange ? $yrange : sprintf("[%g:%g]", @{$range{"llval"}});
  
  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$basename.ps"
    set key off
    set border 15 back
    set xrange $xr
    set yrange $yr
    set grid
    @{[&plot_lltrain_llval_cmd(0 , "llval")]}/;
  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub plot_lltrain_llval_cmd {
  my ($i, $label) = @_;
  my $j = $i + 2;
  my $t = $title ? $title : "lltrain vs. $label";
  my $col = $opts{COLORS}->[$i + 1];

  my $cmd = qq/
    set title "$t"
    set xlabel "lltrain"
    set ylabel "$label"
    f(x) = a*x+b
    fit f(x) "$basename.dat" using 1:$j via a, b
    plot "$basename.dat" using 1:$j notitle with points ps 3 pt 1 lw 2 lc rgb "$col", / .
      qq/f(x) notitle with lines lw 2 lt 1 lc rgb "red"/;
  $cmd =~ s/^\s+//mg;
  return $cmd;
}

sub plot_xy {
  my @plotr;
  my @entries;
  my $opt = $groups[0]->[0];
  foreach my $g (@groups[1 .. $#groups]) {
    my $e = $g->[0];
    if ($e->{"opt"} eq $plot_xy[0]) {
      unless (defined($plotr[$e->{"r"}])) {
        push(@{$plotr[$e->{"r"}]}, $opt);
      }
      push(@{$plotr[$e->{"r"}]}, $e); 
    }
    if ($e->{"score"} > $opt->{"score"}) {
      $opt = $e;
    }
  }

  foreach my $r (@rounds) {
    if (defined($plotr[$r])) {
      push(@entries, @{$plotr[$r]});
    }
  }
  my %range = &get_range(\@entries);
  my $nplots = scalar(@rounds);

  my $cmd = qq/
    set terminal postscript enhanced color linewidth $opts{LINEWIDTH_SCALE} font "$opts{FONT}" $opts{FONTSIZE}
    set output "$basename.ps"
    set multiplot layout $nplots, 1
    set key off
    set border 15 back
    /;
  foreach my $r (@rounds) {
    $cmd .= &plot_xy_cmd($r, $plotr[$r], %range);
  }
  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub plot_xy_cmd {
  my ($round, $entries, %range) = @_;
  my $p = $params{$plot_xy[0]};

  # write data file
  my $file = sprintf("%s.%02d.dat", $basename, $round);
  open FOUT, "> $file" or die "Can't write to '$file'!";
  printf(FOUT "#%7s %8s\n", @plot_xy);
  my @en = sort({ $a->{"params"}->[$p] <=> $b->{"params"}->[$p] } @{$entries});
  foreach my $e (@en) {
    printf(FOUT "%8s %8s\n", $e->{"params"}->[$p], $e->{$plot_xy[1]});
  }
  close FOUT;

  #my %range = &get_range($entries);
  my $xr = $xrange ? $xrange : sprintf("[%g:%g]", @{$range{"params"}->[$p]});
  my $yr = $yrange ? $yrange : sprintf("[%g:%g]", @{$range{$plot_xy[1]}});

  my $cmd = qq/
    set title "Round $round: $plot_xy[0] vs $plot_xy[1]"
    set xlabel "$plot_xy[0]"
    set ylabel "$plot_xy[1]"
    set xrange $xr
    set yrange $yr
    set grid
    plot "$file" using 1:2 notitle with lines lw 3 lc rgb "$opts{COLORS}->[2]", \\
         "$file" using 1:2 title "$plot_xy[1]" with points ps 1 pt 7 lw 3 lc rgb "$opts{COLORS}->[0]"/;
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
  /^#\s+r\s+optimizing\s+(\S+(?:\s+\S+)*?)\s+(lltrain\s+llval)?\s+score\s+\+\/\-%\s*$/ or &invalid_format;
  @_ = split(/\s+/, $1);
  for my $i (0 .. $#_) { $params{$_[$i]} = $i; }
  $nparams = scalar(@_);
  $has_ll = defined($2);
  <FIN>;

  # parse content
  my $k = 0;
  my $r = 0;
  my $opt = "";
  my @params_cur;
  my @group_cur;
  while (<FIN>) {
    last if /^-/;
    my %e;
    /^(.+?)\s+(\S+\s+\S+)?\s+(\S+)\s+(\S+)$/ or last;
    $e{"score"} = $3;
    $e{"gain"} = $4;
    if ($has_ll) {
      $2 or &invalid_format;
      @_ = split(/\s+/, $2);
      $e{"lltrain"} = $_[0];
      $e{"llval"} = $_[1];
    }
    my @f = split(/\s+/, trim($1));
    unless (scalar(@f)) {
      unless (@params_cur) { &invalid_format; }
      $e{"k"} = $k;
      $e{"r"} = $r;
      $e{"opt"} = $opt;
      $e{"params"} = [@params_cur];
    } else {
      if (scalar(@f) < $nparams) { &invalid_format; }
      if (@group_cur) {
        push(@groups, [@group_cur]);
        @group_cur = ();
      }
      @params_cur = @f[scalar(@f) - $nparams ... $#f];
      $e{"params"} = [@params_cur];
      for (1 .. $nparams) { pop(@f); }
      if (scalar(@f) == 3) {
        ($k, $r, $opt) = @f; 
      } elsif (scalar(@f) != 0) { &invalid_format; }
      $e{"k"} = $k;
      $e{"r"} = $r;
      $e{"opt"} = $opt;
    }
    push(@group_cur, \%e);
  }
  if (@group_cur) { push(@groups, [@group_cur]); }
  close FIN;
  scalar(@groups) or &invalid_format;

  # filter entries with invalid score
  my @groups_f;
  foreach my $g (@groups) {
    my @gf;
    foreach my $e (@{$g}) {
      if ($e->{"score"} =~ /^[+-]?\d+\.\d+$/ && $e->{"score"} > 0) { push(@gf, $e); }
    }
    if (scalar(@gf)) { push(@groups_f, [@gf]); }
  }
  @groups = @groups_f;
}

sub get_range {
  my ($entries) = @_;
  my %range;

  my @flds = qw/score gain/;
  if ($has_ll) { push(@flds, qw/lltrain llval/); }
  foreach my $f (@flds) {
    my @r = (&DBL_MAX, -&DBL_MAX);
    foreach my $e (@{$entries}) {
      $r[0] = min($r[0], $e->{$f});
      $r[1] = max($r[1], $e->{$f});
    }
    $range{$f} = [@r];
  }
  for my $i (0 .. $nparams - 1) {
    my @r = (&DBL_MAX, -&DBL_MAX);
    foreach my $e (@{$entries}) {
      $r[0] = min($r[0], $e->{"params"}->[$i]);
      $r[1] = max($r[1], $e->{"params"}->[$i]);
    }
    $range{"params"}->[$i] = [@r];
  }
  
  return %range;
}
