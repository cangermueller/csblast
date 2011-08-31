#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename dirname);
use List::Util qw(shuffle);


=pod
=head1 NAME

  csbenchset.pl - Creates a benchmark set by partitioning a database into k optimization and test sets.

=head1 SYNOPSIS

  csbenchset.pl -d DATABASE -k [1,inf[

  OPTIONS:
    -d, --db-dir DATABASE         Path to the database containing the sequence files to be partitioned
    -K, --folds [1,inf[           Number of folds to be created
    -e, --ext  EXT                File extension of the sequence files [default: seq]
    -o, --out-dir OUTPUT-DIR      Output directory [default: DATABASE]

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $db_dir;
my $K;
my $ext = "seq";
my $out_dir;

my %sfolds;
my @folds;
my $nseqs;


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "d|db-dir=s"  => \$db_dir,
  "K|folds=i"   => \$K,
  "e|ext=s"     => \$ext,
  "o|out-dir=s" => \$out_dir,
  "h|help"      => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($db_dir) { pod2usage("No database provided!"); }
unless (-d $db_dir) { pod2usage("Database does not exist!"); }
$db_dir =~ s/\/$//;
unless ($K) { pod2usage("No number of folds provided!"); }
if ($K < 2) { pod2usage("The number of folds must be greater or equal two!"); }
unless ($out_dir) { $out_dir = dirname($db_dir); }
$out_dir =~ s/\/$//;


### Do the job ###


print "Globbing sequence files ...\n";
$nseqs = 0;
for my $file (glob("$db_dir/*.$ext")) {
  if (`head -n 1 $file` =~ /^>\S+\s+(\w+\.\w+).\w+\.\w+\s*/) {
    push(@{$sfolds{$1}}, basename($file));
    $nseqs++;
  } else { die "Invalid file format: '$file'!"; }
}
print "$nseqs sequence files globbed!\n";

print "Shuffling sequences ...\n";
my @sfold_ids = shuffle(keys(%sfolds));
my $nseqs_per_fold = $nseqs / $K;
my $i = 0;
for my $k (0 .. $K - 1) {
  $folds[$k] = [];
  while ($i < scalar(@sfold_ids) && ($k == $K - 1 || 
         $nseqs_per_fold - scalar(@{$folds[$k]}) > 0.5 * scalar(@{$sfolds{$sfold_ids[$i]}}))) {
    push(@{$folds[$k]}, @{$sfolds{$sfold_ids[$i++]}});
  }
}

for my $k (0 .. $K - 1) {
  printf "Creating fold %02d ...\n", $k + 1;
  my $out_base = sprintf("$out_dir/%s_%02d", basename($db_dir), $k + 1);
  &create_set("${out_base}_test", [$k]);
  my @opt = (0 .. $K - 1);
  splice(@opt, $k, 1);
  &create_set("${out_base}_opt", \@opt);
}
print "Done!\n";


### Utils ###


sub create_set {
  my ($dir, $fold_ids) = @_;
  &exec("rm -rf $dir; mkdir $dir");
  for my $k (@{$fold_ids}) {
    for my $f (@{$folds[$k]}) {
      my $file = "$db_dir/$f";
      &exec("cat $file >> $dir/db");
      $file =~ s/$ext$//;
      &exec("cp $file* $dir");
    }
  }
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd")) {
    die "Error executing '$cmd'!"; 
  }
}
