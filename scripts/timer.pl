#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Time::HiRes qw(gettimeofday tv_interval);
use File::Basename qw(basename);
use File::Temp qw(tempfile);
use My::Config;


=pod

=head1 NAME
 
  timer.pl - Measure time for CSBLAST pseudocounts generation before a BLAST call.

=head1 SYNOPSIS

  timer.pl [OPTIONS] -s SEQ-DIR -D CONTEXT-DATA+

  OPTIONS:
    -s, --seq-dir DIR           Directory containing sequence files
    -e, --ext EXT               Sequence files extension [def: seq]
    -n, --num [1,inf[           Number of sequence to be used for time measurement [def: 50]
    -L, --length [1,inf[+       List of sequence lengths to be used [default: 100 200 400 800 1600]
    -D, --context-data FILE+    Models with context-data
    -x, --pc-admix ]0,1]+       Pseudocount admix coefficients for context-specific pseudocounts [default: 0.90]
    -z, --pc-neff [1,inf[+      Target Neff for pseudocounts admixture [def: ]
    -o, --out-dir DIR           Directory for output files
    -p, --parallel              Activate parallelization [default: true]
    -h, --help                  Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut

 
### Variables ###


my $seq_dir;
my $seq_ext    = "seq";
my $seq_num    = 50;
my @seq_length;
my @models;
my @pc_x;
my @pc_z;
my $parallel   = 1;
my $out_dir;

my %cfg = new My::Config;
my @pc_args;
my @seqs;
my @seq_files;
my %time;
my $sep = "    ";
my $bin;


### Initialization ###


Getopt::Long::Configure("no_ignore_case");
GetOptions(
  "s|seq-dir=s"          => \$seq_dir,
  "e|ext=s"              => \$seq_ext,
  "n|num=i"              => \$seq_num,
  "L|length=i{1,}"       => \@seq_length,
  "D|context-data=s{1,}" => \@models,
  "x|pc-admix=f{1,}"     => \@pc_x,
  "z|pc-neff=f{1,}"      => \@pc_z,
  "o|out-dir=s"          => \$out_dir,
  "p|parallel!"          => \$parallel,
  "h|help"               => sub { pod2usage(2); }
) or pod2usage(1);
unless ($seq_dir) { pod2usage("Sequence directory not provided!"); }
unless (@models) { pod2usage("No context-data provided!"); }
foreach my $m (@models) {
  unless (-f $m) { pod2usage("'$m' does not exist!"); }
}
unless (@seq_length) { @seq_length = (100, 200, 400, 800, 1600); }
@seq_length = sort({ $a <=> $b } @seq_length);
unless (scalar(@pc_x) || scalar(@pc_z)) { @pc_x = (0.9); }
foreach my $x (@pc_x) { push(@pc_args, "-x $x"); }
foreach my $z (@pc_z) { push(@pc_args, "-z $z"); }
$out_dir =~ s/\/$//;
$bin = $parallel ? "csblast" : "csblast_debug";


### Prepare sequences ###


print "Preparing sequences for time measuring...\n";
foreach my $file (glob("$seq_dir/*$seq_ext")) {
  my $seq = "";
  open FIN, "< $file" or die "Can't read from '$file'!";
  my $h;
  while (<FIN>) {
    chomp;
    if (/^(>.+)$/) {
      if ($h) { last; }
      else { $h = $1; }
    } elsif ($h) {
      s/^\s+//;
      s/\s+$//;
      $seq .= $_;
    }
  }
  close FIN;
  my $l = $seq_length[$#seq_length];
  if (length($seq) < $l / 5) {  next; }
  $seq = $seq x (int($l / length($seq)) + 1);
  push(@seqs, [$h, $seq]);
  if (scalar(@seqs) == $seq_num) { last; }
}
if (scalar(@seqs) < $seq_num) { die "Sequence directory does not contain sufficient sequences!"; }
($_, my $base) = tempfile("timerseqXXXXX", DIR => "/tmp");
close $_;
system("rm -f $base");
for my $i (0 .. $seq_num - 1) {
  push(@seq_files, "$base.$i.seq");
}


### Time measurement ###


printf("\n#%8s$sep%8s$sep%8s$sep%s\n", "TIME", "L", "ADMIX", "MODEL");
foreach my $L (@seq_length) {
  for my $i (0 .. $seq_num - 1) {
    open FOUT, "> $seq_files[$i]" or die "Can't write to '$seq_files[$i]'!";
    printf(FOUT "%s\n%s", $seqs[$i]->[0], substr($seqs[$i]->[1], 0, $L));
    close FOUT;
  }
  foreach my $m (@models) {
    foreach my $pc (@pc_args) {
      my $time = 0;
      foreach my $seq_file (@seq_files) {
        my $cmd = qq/$bin -D $m -i $seq_file $pc --blast-path $cfg{"blast_path"} --emulate &> \/dev\/null/;
        my $t0 = [gettimeofday];
        system($cmd);
        my $t1 = [gettimeofday];
        if ($? != 0) { die "Error calling csblast: '$!'"; }
        my $td = tv_interval($t0, $t1);
        $time += $td;
      }
      $time /= scalar(@seqs);
      $_ = $pc;
      s/^-//;
      s/ //;
      push(@{$time{$m}->{$_}}, $time);
      printf("%9.4f$sep%8d$sep%8s$sep%s\n", $time, $L, $_, $m);
    }
  }
}
foreach my $seq_file (@seq_files) { system("rm -f $seq_file"); }


### Write datafiles ###


print "\n";
if ($out_dir) {
  foreach my $m (keys(%time)) {    
    foreach my $pc (keys(%{$time{$m}})) {
      my $file = sprintf("%s/%s_%s.dat", $out_dir, basename($m), $pc);
      open FOUT, "> $file" or die "Can't write to '$file'!";
      printf(FOUT "#%7s$sep%8s\n", "L", "TIME");
      for my $i (0 .. $#{$time{$m}->{$pc}}) {
        printf(FOUT "%8d$sep%8.4f\n", $seq_length[$i], $time{$m}->{$pc}->[$i]);
      }
      close FOUT;
    }
  }
}
print "Done!\n";
