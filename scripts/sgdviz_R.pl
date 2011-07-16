#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use My::Utils qw(remove_suffix filename);

=pod
=head1 NAME

  sgdviz.pl - Visualizes parameters of an SGD optimization.

=head1 SYNOPSIS

  sgdviz.pl [OPTIONS] --infile INFILE

  OPTIONS:
  -i, --infile INFILE     SGD log file.
  -o, --outfile OUTFILE   Output file [default: INFILE.pdf]
  -p, --param PARAM+      Parameters to be plotted [default: ll-train,ll-val,prior,neff]
  -n, --name NAME         Name of the SGD model.
  -k, --keep              Keep data and plot file [default: false]
  -h, --help              Show this help message.

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $infile;
my $outfile;
my $basename;
my @params;
my $name;
my $keep;
my @lines;


### Initialization ###


GetOptions(
  "i|infile=s"     => \$infile,
  "o|outfile=s"    => \$outfile,
  "p|params=s{1,}" => \@params,
  "n|name=s"       => \$name,
  "keep!"          => \$keep,
  "h|help"         => sub { pod2usage(2); }
) or die pos2usage(1);
unless ($infile) { pod2usage("Not input file provided!"); }
if ($outfile) { 
  $basename = remove_suffix($outfile); 
} else {
  $basename = filename($infile);
  $outfile = "$basename.pdf";
}
unless (@params) { @params = qw(ll-train prior ll-val neff); }
unless ($name) { $name = filename($infile); }
@params = map({
    if ($_ eq "ll-train") { 2; }
    elsif ($_ eq "prior") { 3; }
    elsif ($_ eq "ll-val") { 4; }
    elsif ($_ eq "neff") { 5; }
    else { die "Invalid plot parameter!"; }
  } @params);
print "Creating '$outfile'...\n";


### Write data file ###


#Epoch                  LL-Train    Prior            LL-Val     Neff
#1     [==============]   2.3253  +1.0918  -0.0325   2.1410  14.5634
open FIN, "< $infile" or die "Can't read from '$infile'!";
while (<FIN>) {
  if (my @line = /^(\d+)\s+\[=+\]\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) { 
    $line[2] *= -1;
    push(@lines, \@line); 
  } elsif (scalar(@lines)) { last; }
}
close FIN;
unless (scalar(@lines)) { die "Invalid input format!"; }

@lines = sort({$a->[0] <=> $b->[0]} @lines);
open FOUT, "> $basename.sgd" or die "Can't write to '$basename.sgd'";
printf(FOUT "#%9s%10s%10s%10s%10s\n", "Epoch", "LL-Train", "Prior", "LL-Val", "Neff");
foreach my $line (@lines) {
  foreach my $l (@{$line}) { printf(FOUT "%10s", $l); }
  print(FOUT "\n");
}
close FOUT;

open FOUT, "> $basename.R" or die "Can't write to $basename.R";
print FOUT <<END;
library(calibrate)

opt.infile <- NA
opt.outfile <- NA
opt.title <- NA
opt.cols <- NA
opt.cols.total <- 5
opt.cols.neff <- 5
opt.cols.prior <- 3
opt.cols.y2 <- c(opt.cols.prior, opt.cols.neff)

opt.size <- c(8, 8)
opt.legend.cex <- 1.0
opt.linewidth <- 1.8
opt.linewidth.y2 <- 1.2
opt.labels <- c("", "LL-Train", "Prior", "LL-Val", "Neff")
opt.colors <- c("", "green", "blue", "red", "darkmagenta")
opt.linetypes <- c(0, 1, 6, 5, 2)
opt.pch <- rep(1, opt.cols.total)
opt.abline.linewidth <- 2.0
opt.abline.linetype <- 3

opt.args <- commandArgs(TRUE)
if (length(opt.args) < 4) {
  print("plot.R INFILE OUTFILE TITLE COL+")
  quit(status=1)
}
opt.infile <- opt.args[1]
opt.outfile <- opt.args[2]
opt.title <- opt.args[3]
opt.cols.s <- opt.args[4:length(opt.args)]
opt.cols.n <- length(opt.cols.s)
opt.cols.has <- array(F, opt.cols.total)
opt.cols <- c()
for (i in 1:opt.cols.n) {
  opt.cols[i] <- as.integer(opt.cols.s[i])
  opt.cols.has[opt.cols[i]] <- T
}

data.lab <- read.table(opt.infile, comment.char="", nrow=1)
data.lab <- data.lab[1,2:ncol(data.lab)]
data <- read.table(opt.infile, skip=1)
plot.x <- (data[,1])

pdf(opt.outfile, title=opt.title, opt.size[1], opt.size[2])

for (i in 1:length(opt.cols.y2)) {
  c <- opt.cols.y2[i]
  if (opt.cols.has[c]) {
    plot(plot.x, data[,c], type="l", lwd=opt.linewidth.y2, xaxt="n", yaxt="n",
        xlab="SGE epoch", ylab="", lty=opt.linetypes[c], col=opt.colors[c], 
        pch=opt.pch[c], ylim=c(min(data[,c]), max(data[,c])))
    par(new=T)
  }
}
for (i in length(opt.cols.y2):1) {
  c <- opt.cols.y2[i]
  if (opt.cols.has[c]) {
    axis(4, pretty(c(min(data[,c]), max(data[,c]))), col=opt.colors[c])
    break
  }
}   

plot.y <- c(min(data[,2], data[,4]), max(data[,2], data[,4]))
for (i in 1:opt.cols.n) { 
  c <- opt.cols[i]
  if (! c %in% opt.cols.y2) {
    par(new=T)
    plot(plot.x, data[,c], type="l", xlab="", ylab="", axes=F,
        lwd=opt.linewidth, lty=opt.linetypes[c], col=opt.colors[c], 
        pch=opt.pch[c], ylim=plot.y)
    max.x <- 0
    max.y <- 0
    for (i in 1:length(data[,1])) {
      if (data[i, c] > max.y) {
        max.x <- data[i,1]
        max.y <- data[i,c]
      }
    }
    abline(h=max.y, lty=opt.abline.linetype, lwd=opt.abline.linewidth,
        col=opt.colors[c])
    m <- max(plot.x) - 0.08 * (max(plot.x) - min(plot.x))
    textxy(min(m, max.x), max.y, 
        sprintf("%.3f", max.y), cx=1.0, dcol=opt.colors[c])
  }
}
axis(1, at=plot.x, labels=plot.x)
axis(2, pretty(plot.y))
if (length(opt.cols.y2) == 0) box()

legend.labels <- c()
legend.colors <- c()
legend.linetypes <- c()
for (i in 1:opt.cols.n) {
  c <- opt.cols[i]
  legend.labels[i] <- opt.labels[c]
  legend.colors[i] <- opt.colors[c]
  legend.linetypes[i] <- opt.linetypes[c]
}

legend("right", hor=F, cex=opt.legend.cex, bty="n",
	legend.labels, col=legend.colors, lty=legend.linetypes, lwd=opt.linewidth)
END
close FOUT;
`R --vanilla --slave -f $basename.R --args $basename.sgd $outfile $name @params &> /dev/null`;
if ($?) { die "Error calling '$basename.R'!"; }
unless ($keep) { `rm -f $basename.dat $basename.sgd $basename.R`; }
print "Done!\n";
