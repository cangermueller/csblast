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
 
  csruntime.pl - Measure time for CSBLAST pseudocounts generation before a BLAST call.
                 For creating the publication results, csruntime was used.

=head1 SYNOPSIS

  csruntime.pl [OPTIONS] -s SEQ-GLOB -D CONTEXT-DATA+ [CSBLAST OPTIONS]

  OPTIONS:
    -s, --seq-glob DIR          Glob of input sequences
    -n, --num [1,inf[           Number of sequence to be used for time measurement [def: 50]
    -L, --length [1,inf[+       List of sequence lengths to be used [default: 100 200 400 800 1600]
    -D, --context-data FILE+    Models with context-data
    -x, --pc-admix ]0,1]+       Pseudocount admix coefficients for context-specific pseudocounts [def: 0.90]
    -z, --pc-neff [1,inf[+      Target Neff for pseudocounts admixture [def: off]
    -o, --out-dir DIR           Directory for output files [def: .]
    -b, --binary FILE           CSBLAST binary file [def: csblast]
    -d, --db-file FILE          Database to be used for calling blastpgp [def: off]
    -h, --help                  Show this help message

=head1 AUTHOR

  Angermueller Christof
  angermue@in.tum.de

=cut

 
### Variables ###


my $seq_glob;
my $seq_num = 50;
my @seq_length;
my @models;
my @pc_x;
my @pc_z;
my $binary = "csblast";
my $db_file;
my @csblast_opts;
my $out_dir;

my %cfg = new My::Config;
my @pc_args;
my @seqs;
my @seq_files;
my %time;
my $sep = "    ";
my $bin;


### Initialization ###


Getopt::Long::Configure("no_ignore_case", "pass_through");
GetOptions(
  "s|seq-glob=s"         => \$seq_glob,
  "n|num=i"              => \$seq_num,
  "L|length=i{1,}"       => \@seq_length,
  "D|context-data=s{1,}" => \@models,
  "x|pc-admix=f{1,}"     => \@pc_x,
  "z|pc-neff=f{1,}"      => \@pc_z,
  "o|out-dir=s"          => \$out_dir,
  "b|binary=s"           => \$binary,
  "d|db-file=s"          => \$db_file,
  "h|help"               => sub { pod2usage(2); }
) or pod2usage(1);
@csblast_opts = @ARGV;
unless ($seq_glob) { pod2usage("Sequence glob not provided!"); }
unless (@models) { pod2usage("No context-data provided!"); }
foreach my $m (@models) {
  unless (-f $m) { pod2usage("'$m' does not exist!"); }
}
unless (@seq_length) { @seq_length = (100, 200, 400, 800, 1600); }
@seq_length = sort({ $a <=> $b } @seq_length);
unless (scalar(@pc_x) || scalar(@pc_z)) { @pc_x = (0.9); }
foreach my $x (@pc_x) { push(@pc_args, "-x $x"); }
foreach my $z (@pc_z) { push(@pc_args, "-z $z"); }
if ($out_dir) { $out_dir =~ s/\/$//; }


### Prepare sequences ###


print "Preparing sequences for time measuring...\n";
foreach my $file (glob("$seq_glob")) {
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
  if (length($seq) < $l / 5) {  next; } # sequence is shorter than one fith of the maximal length
  $seq = $seq x (int($l / length($seq)) + 1); # elongate the sequence
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
  # Write sequences with length $L
  for my $i (0 .. $seq_num - 1) {
    open FOUT, "> $seq_files[$i]" or die "Can't write to '$seq_files[$i]'!";
    printf(FOUT "%s\n%s", $seqs[$i]->[0], substr($seqs[$i]->[1], 0, $L));
    close FOUT;
  }
  foreach my $m (@models) {
    foreach my $pc (@pc_args) {
      my $time = 0;
      foreach my $seq_file (@seq_files) {
        my $cmd = sprintf(qq/$binary -D $m -i $seq_file $pc --blast-path $cfg{"blast_path"} %s %s &> \/dev\/null/,
          ($db_file ? "-d $db_file" : "--emulate"),
          join(" ", @csblast_opts));
        my $t0 = [gettimeofday];
        system($cmd);
        my $t1 = [gettimeofday];
        if ($? != 0) { 
          print STDERR "Error calling '$cmd'!\n";
          exit $?
        }
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
