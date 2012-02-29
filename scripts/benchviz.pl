#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename);
use Cwd qw(abs_path);


=pod

=head1 NAME

  benchviz.pl - Visualizes CS-BLAST benchmark results.

=head1 SYNOPSIS

  benchviz.pl [OPTIONS] --dir BENCHDIR+

  OPTIONS:
 
    -d, --dir BENCHDIR+     List of directories containing benchmark result files
    -o, --out OUTBASE       Output basename [def: ./]
    -p, --plot PLOT+        List of plots to be created [default: tpfp wtpfp rocx evalue]
    -t, --title TITLE       Title of the plots [def: ]
    -l, --label LABEL+      List of labels to be used instead of directory names
        --db DB             Database used for scaling axis
        --iter INT          Number of CSI-BLAST iterations [def: 1]
        --no-under          Substitute underscore in labels [def: false]
        --[no-]sort         Sort entities by their ROC score [def: true]
    -k, --keep              Keep plot files [def: false]
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
my $db;
my $iter;
my $no_under;
my $sort = 1;
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
  "db=s"          => \$db,
  "iter=i"        => \$iter,
  "no-under"      => \$no_under,
  "sort!"         => \$sort,
  "k|keep"        => \$keep,
  "h|help"        => sub { pod2usage(2); }
) or pod2usage(1);
unless (@dirs) { pod2usage("No benchmark directory provided!"); }
unless (@plots) { @plots = qw/tpfp wtpfp rocx evalue/; }
&get_entities;
unless (@entities) { pod2usage("No plot data found in the specified benchmark directories!"); }
unless ($db) {
  if (abs_path($dirs[0]) =~ /(scop20_1\.7._(opt|test))/) { $db = $1; }
  else { $db = "scop20_1.73_test"; }
}
unless ($iter) {
  $iter = 1;
  foreach my $d (@dirs) {
    if (basename($d) =~ /(?:\A|_)j(\d)(?:\z|_)/ && $1 > $iter) { $iter = $1; }
  }
}


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
    set style line 3 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[2]"
    set style line 4 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[3]"
    set style line 5 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[4]"
    set style line 6 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[5]"
    set style line 7 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[6]"
    set style line 8 linetype 1 linewidth $opts{LINEWIDTH} linecolor rgb "$opts{COLORS}->[7]"
    set style line 9 linetype 0 linewidth $opts{LINEWIDTHE} linecolor rgb "$opts{COLORS}->[8]"/;

  if ($plot eq "tpfp") { $cmd .= qq/
    set key top left reverse invert Left
    set log x
    set grid
    set xlabel "FP"
    set ylabel "TP"/;
    if ($db eq "scop20_1.73_opt" && $iter == 1) { $cmd .= qq/
      set xrange [1:4000]
      set yrange [0:10000]
      set label "1%" at 40,7500 textcolor ls 9
      set label "10%" at 500,7500 textcolor ls 9
      set label "20%" at 1150,7500 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_opt" && $iter == 2) { $cmd .= qq/
      set xrange [1:4000]
      set yrange [0:12000]
      set label "1%" at 40,7500 textcolor ls 9
      set label "10%" at 500,7500 textcolor ls 9
      set label "20%" at 1150,7500 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_opt") { $cmd .= qq/
      set xrange [1:4000]
      set yrange [0:14000]
      set label "1%" at 40,7500 textcolor ls 9
      set label "10%" at 500,7500 textcolor ls 9
      set label "20%" at 1150,7500 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 1) { $cmd .= qq/
      set xrange [1:6000]
      set yrange [0:10000]
      set label "1%" at 35,5000 textcolor ls 9
      set label "10%" at 350,5000 textcolor ls 9
      set label "20%" at 800,5000 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 2) { $cmd .= qq/
      set xrange [1:6000]
      set yrange [0:20000]
      set label "1%" at 40,7500 textcolor ls 9
      set label "10%" at 400,7500 textcolor ls 9
      set label "20%" at 1100,7500 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test") { $cmd .= qq/
      set xrange [1:6000]
      set yrange [0:25000]
      set label "1%" at 40,7500 textcolor ls 9
      set label "10%" at 400,7500 textcolor ls 9
      set label "20%" at 1100,7500 textcolor ls 9/;
    } elsif ($db eq "scop20_1.75_opt") { $cmd .= qq/
      set xrange [1:4000]
      set yrange [0:8000]
      set label "1%" at 45,6500 textcolor ls 9
      set label "10%" at 450,6500 textcolor ls 9
      set label "20%" at 1000,6500 textcolor ls 9/;
    } else { $cmd .= qq/
      set xrange [1:4000]
      set yrange [0:10000]
      set label "1%" at 50,7500 textcolor ls 9
      set label "10%" at 500,7500 textcolor ls 9
      set label "20%" at 1100,7500 textcolor ls 9/;
    }

  } elsif ($plot eq "wtpfp") { $cmd .= qq/
    set key top left reverse invert Left
    set log x
    set grid
    set xlabel "weighted FP"
    set ylabel "weighted TP"/;
    if ($db eq "scop20_1.73_opt" && $iter == 1) { $cmd .= qq/
      set xrange [1:200]
      set yrange [0:250]
      set label "1%" at 1.2,175 textcolor ls 9
      set label "10%" at 12,175 textcolor ls 9
      set label "20%" at 30,175 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_opt" && $iter == 2) { $cmd .= qq/
      set xrange [1:200]
      set yrange [0:400]
      set label "1%" at 2,275 textcolor ls 9
      set label "10%" at 20,275 textcolor ls 9
      set label "20%" at 50,275 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_opt") { $cmd .= qq/
      set xrange [1:200]
      set yrange [0:500]
      set label "1%" at 2,275 textcolor ls 9
      set label "10%" at 20,275 textcolor ls 9
      set label "20%" at 50,275 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 1) { $cmd .= qq/
      set xrange [1:600]
      set yrange [0:800]
      set label "1%" at 4,550 textcolor ls 9
      set label "10%" at 40,550 textcolor ls 9
      set label "20%" at 95,550 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 2) { $cmd .= qq/
      set xrange [1:600]
      set yrange [0:1200]
      set label "1%" at 4,550 textcolor ls 9
      set label "10%" at 40,550 textcolor ls 9
      set label "20%" at 95,550 textcolor ls 9/;
    } elsif ($db eq "scop20_1.73_test") { $cmd .= qq/
      set xrange [1:600]
      set yrange [0:1200]
      set label "1%" at 4,550 textcolor ls 9
      set label "10%" at 40,550 textcolor ls 9
      set label "20%" at 95,550 textcolor ls 9/;
    } elsif ($db eq "scop20_1.75_opt") { $cmd .= qq/
      set xrange [1:600]
      set yrange [0:800]
      set label "1%" at 2,275 textcolor ls 9
      set label "10%" at 20,275 textcolor ls 9
      set label "20%" at 50,275 textcolor ls 9/;
    } else { $cmd .= qq/
      set xrange [1:600]
      set yrange [0:600]
      set label "1%" at 2.5,350 textcolor ls 9
      set label "10%" at 25,350 textcolor ls 9
      set label "20%" at 60,350 textcolor ls 9/;
    }

  } elsif ($plot eq "fdr") { $cmd .= qq/
    set key top left reverse invert Left
    set grid
    set xlabel "FDR"
    set ylabel "TPR"/;
    if ($db eq "scop20_1.73_opt") { $cmd .= qq/
      set xrange [0:1]
      set yrange [0:1]/;
    } elsif ($db eq "scop20_1.73_test") { $cmd .= qq/
      set xrange [0:1]
      set yrange [0:1]/;
    } elsif ($db eq "scop20_1.75_opt") { $cmd .= qq/
      set xrange [0:1]
      set yrange [0:1]/;
    } else { $cmd .= qq/
      set xrange [0:1]
      set yrange [0:1]/;
    }

  } elsif ($plot eq "rocx") { $cmd .= qq/
    set key top right invert
    set grid
    set xlabel "ROC"
    set ylabel "Fraction of queries"
    set xrange [0:1.0]/;
    if ($db eq "scop20_1.73_opt" && $iter == 1) { $cmd .= qq/
      set yrange [0:0.6]/;
    } elsif ($db eq "scop20_1.73_opt" && $iter == 2) { $cmd .= qq/
      set yrange [0:0.7]/;
    } elsif ($db eq "scop20_1.73_opt") { $cmd .= qq/
      set yrange [0:0.8]/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 1) { $cmd .= qq/
      set yrange [0:0.6]/;
    } elsif ($db eq "scop20_1.73_test" && $iter == 2) { $cmd .= qq/
      set yrange [0:0.75]/;
    } elsif ($db eq "scop20_1.73_test") { $cmd .= qq/
      set yrange [0:0.8]/;
    } elsif ($db eq "scop20_1.75_opt") { $cmd .= qq/
      set yrange [0:0.5]/;
    } else { $cmd .= qq/
      set yrange [0:0.5]/;
    }

  } elsif ($plot eq "evalue") { $cmd .= qq/
    set key top left reverse invert Left
    set grid
    set xlabel "Reported E-value"
    set ylabel "Actual E-value"
    set xrange [1e-4:1e2]
    set yrange [1e-4:1e3]
    set logscale x
    set logscale y
    set tics format "%.0e"/;

  } elsif ($plot eq "pvalue") { $cmd .= qq/
    set key top left reverse invert Left
    set grid
    set xlabel "Reported P-value"
    set ylabel "Actual P-value"
    set xrange [1e-7:1]
    set yrange [1e-7:1]
    set logscale x
    set logscale y
    set tics format "%.0e"/;
  } 

  $cmd .= qq/
    plot /;
  if ($plot eq "tpfp" || $plot eq "wtpfp") { 
    $cmd .= qq/99*x notitle with lines ls 9, 9*x notitle with lines ls 9, 4*x notitle with lines ls 9,/;
  } elsif ($plot eq "evalue" || $plot eq "pvalue") {
    $cmd .= qq/x notitle with lines ls 9,/;
  }
  $cmd .= &cmd_curves($plot);
  $cmd =~ s/^\s+//mg;
  return "$cmd\n";
}

sub cmd_curves {
  my ($plot) = @_;
  my @curves;
  my $ls = 2;
  for my $i (0 .. $#entities) {
    my $e = $entities[$i];
    if ($e->{PLOTS}->{$plot}) {
      push(@curves, sprintf('"%s" title "%s" with lines ls %d', 
          &get_datafile($plot, $e->{DIR}), ($plot eq "evalue" ? "" : sprintf("%.3f: ", $e->{ROCX})) . $e->{LABEL},
          basename($e->{DIR}) =~ /^(blast|psiblast)/ ? 1 : $ls++));
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
      if ($no_under) { $e{LABEL} =~ s/_/ /g; }
    } else {
      $e{LABEL} = basename($e{DIR}); 
      $e{LABEL} =~ s/_/ /g;
    }
    my $has_file = 0;
    foreach my $p (@plots) {
      my $h =  -f sprintf("%s/%s", $dirs[$i], &get_datafile($p)) ? 1 : 0;
      if ($h) { $has_file = 1; }
      $e{PLOTS}->{$p} = $h;
    }
    if ($has_file) { push(@unsorted, \%e); }
  }
  if ($sort) {
    @entities = sort({ $a->{ROCX} <=> $b->{ROCX} } @unsorted);
  } else {
    @entities = @unsorted;
  }
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
