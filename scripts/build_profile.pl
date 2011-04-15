#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Math::Random::OO::Normal;
use File::Temp qw(tempfile);
use File::Basename qw(basename fileparse);


=pod 
=head1 NAME

 build_profile.pl - Samples an alignment from a given set of alignments and builds a count profile

=head1 SYNOPSIS

 build_profile.pl [OPTIONS] --input-ali INPUT-ALI+

 OPTIONS:
  -i, --input-ali INPUT-ALI+    Input alignments (A2M, A3M, FASTA)
  -o, --output-prf OUTPUT-PRF   Output profile [default: INPUT-ALI.prf]
  -m, --min-neff MIN-NEFF       Minimum Neff
  -x, --max-neff MAX-NEFF       Maximum Neff
  -s, --sigma SIGMA             Sigma for sampling from gaussian distribution 
                                [default: equal distribution]
  -f, --filter                  Filter diverity if the minimum Neff is higher than MAX-NEFF 
                                [default: true]
  -h, --help                    Shows this help message

=head1 AUTHOR

 Angermueller Christof
 angermue@in.tum.de

=cut


### Variables ###


my @ali_files;
my @ali_neff;
my $out_prf;
my $min_neff;
my $max_neff;
my $sigma;
my $filter = 1;
my $tmp_file;


### Initialization ###


GetOptions( 
    "f|filter!"         => \$filter,
    "h|help"            => sub { pod2usage(2); },
    "i|input-ali=s{1,}" => \@ali_files,
    "m|min-neff=f"      => \$min_neff,
    "x|max-neff=f"      => \$max_neff,
    "o|output-prf=s"    => \$out_prf,
    "s|sigma=f"         => \$sigma
) or pod2usage(1);

unless (scalar(@ali_files)) { pod2usage("Input a3m files missing!"); }
if (defined($min_neff) && $min_neff < 1) { pod2usage("Neff must greater or equal 1!"); }
if (defined($max_neff) && $max_neff < 1) { pod2usage("Neff must greater or equal 1!"); }
if (defined($max_neff) && defined($min_neff) && $min_neff > $max_neff) { pod2usage("Neff range invalid!"); }


### Compute Neff in alignments ###


print("Computing Neff...\n");
foreach my $f (@ali_files) {
    push(@ali_neff, [$f, &neff($f)]);
}
@ali_neff = sort({ $a->[1] <=> $b->[1] } @ali_neff);
unless (defined($min_neff)) { 
    $min_neff = $ali_neff[0]->[1]; 
    if (defined($max_neff) && $min_neff > $max_neff) { $min_neff = $max_neff; }
}
unless (defined($max_neff)) { 
    $max_neff = $ali_neff[$#ali_neff]->[1]; 
    if (defined($min_neff) && $min_neff > $max_neff) { $max_neff = $min_neff; }
}


### Sample alignemnt ###


my $i;
if ($max_neff < $ali_neff[0]->[1]) {
    $i = 0;
    if ($filter) {        
        (undef, $tmp_file) = tempfile("filteredXXXX", SUFFIX => ".a3m", TMPDIR => 1);
        `hhfilter -i $ali_neff[0]->[0] -o $tmp_file -neff $max_neff`;
        unshift(@ali_neff, [$tmp_file, &neff($tmp_file)]);
    }

} elsif ($min_neff > $ali_neff[$#ali_neff]->[1]) { 
    $i = $#ali_neff; 

} else {
    my $neff;
    if ($sigma) {
        $neff = Math::Random::OO::Normal->new(0.5 * ($max_neff + $min_neff), $sigma)->next;
    } else {
        $neff = $min_neff + rand() * ($max_neff - $min_neff);
    }
    $i = 0;
    while ($i < $#ali_neff) {
        if (abs($neff - $ali_neff[$i + 1]->[1]) <= abs($neff - $ali_neff[$i]->[1])) { $i++; }
        else { last; }
    }
}


### Creating count profile ###


printf("Creating count profile from %s...\n", basename($ali_neff[$i]->[0]));
unless ($out_prf) { $out_prf = sprintf("%s.prf", fileparse($ali_neff[$i]->[0], qr/\..+$/)); }
system("csbuild -i $ali_neff[$i]->[0] -o $out_prf -x 0");
if ($tmp_file) { system("rm -f $tmp_file"); }
printf("Count profile created with Neff %1.2f!\n", $ali_neff[$i]->[1]);
       




### Functions ###

sub neff {
    my ($file) = @_;
    `hhmake -i $file -o /dev/null` =~ /exp\(entropy\)\s*=\s*(\S+)/ or die "Error computing Neff for $file";
    return $1;
}
