#!/usr/bin/perl -w

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename dirname);
use List::Util qw(shuffle);
use My::Config;


=pod
=head1 NAME

  db_benchset.pl - Creates a benchmark set by randomly partitioning a database into K blocks containing
                   sequences of different folds and by using the first block as optimization set and 
                   all remaining blocks as test set. If --crossval is set, K optimization and test sets 
                   are created where block i is used as test set and the blocks [K] \ {i} as optimization set.

=head1 SYNOPSIS

  db_benchset.pl -d DATABASE -K [1,inf[

  OPTIONS:
    -d, --db-dir DATABASE         Path to the database containing the sequence files to be partitioned
    -K, --folds [1,inf[           Number of folds blocks
    -c, --crossval                Use cross-validation scheme
    -e, --ext  EXT                File extension of the sequence files [default: seq]
    -o, --out-base OUT-BASE       Output basename [default: DATABASE]

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut


### Variables ###


my $db_dir;
my $K;
my $ext = "seq";
my $crossval;
my $outbase;

my %sfolds;
my @folds;
my $nseqs;
my %cfg = new My::Config;

### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "d|db-dir=s"  => \$db_dir,
  "K|folds=i"   => \$K,
  "e|ext=s"     => \$ext,
  "c|crossval!" => \$crossval,
  "o|outbase=s" => \$outbase,
  "h|help"      => sub { pod2usage(2); }
) or die pod2usage(1);
unless ($db_dir) { pod2usage("No database provided!"); }
unless (-d $db_dir) { pod2usage("Database does not exist!"); }
$db_dir =~ s/\/$//;
unless ($K) { pod2usage("No number of folds provided!"); }
if ($K < 2) { pod2usage("The number of folds must be greater or equal two!"); }
unless ($outbase) { $outbase = $db_dir; }
$outbase =~ s/\/$//;


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

if ($crossval) {
  for my $k (0 .. $K - 1) {
    printf "Creating fold %02d ...\n", $k + 1;
    my $outdir = sprintf("%s_%02d", $outbase,  $k + 1);
    &create_set("${outdir}_test", [$k]);
    my @opt = (0 .. $K - 1);
    splice(@opt, $k, 1);
    &create_set("${outdir}_opt", \@opt);
  }
} else {
  print "Creating ${outbase}_opt ...\n";
  &create_set("${outbase}_opt", [0]);
  print "Creating ${outbase}_test ...\n";
  &create_set("${outbase}_test", [1 .. $K]);
}

print "Done!\n";


### Utils ###


sub create_set {
  my ($outbase, $fold_ids) = @_;
  my $db = $outbase;
  my $db_file = "${outbase}_db";
  # copy sequence files and create databases file
  &exec("rm -rf $db; mkdir $db");
  for my $k (@{$fold_ids}) {
    for my $f (@{$folds[$k]}) {
      my $file = "$db_dir/$f";
      &exec("cat $file >> $db_file");
      $file =~ s/$ext$//;
      &exec("cp $file* $db");
    }
  }
  # create blast databases files
  my $title = basename($outbase);
  &exec("$cfg{'blast_path'}/formatdb -i $db_file -p T -t '$title' -l /dev/null");
}

sub exec {
  my ($cmd) = @_;
  if (system("$cmd")) {
    die "Error executing '$cmd'!"; 
  }
}
