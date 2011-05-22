#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
# use Math::Random::OO::Normal;
use File::Temp qw(tempfile);
use File::Basename qw(basename fileparse);
use My::Utils qw(remove_suffix suffix filename);


=pod 
=head1 NAME

 build_profile.pl - Samples an alignment from a given set of alignments and builds a count profile

=head1 SYNOPSIS

 build_profile.pl [OPTIONS] --input-ali INPUT-ALI+

 OPTIONS:
  -i, --input-ali INPUT-ALI+    Input alignments (A2M, A3M, FASTA)
  -o, --outbase OUTPUT-BASE     Output basename [default: INPUT-ALI]
  -m, --min-neff MIN-NEFF       Minimum Neff
  -x, --max-neff MAX-NEFF       Maximum Neff
  -d, --mode-neff MODE-NEFF     Mode Neff when sampling from the triangular distribution
                                [default: equal distribution]
  -f, --filter                  Filter diverity if the minimum Neff is higher than MAX-NEFF 
                                [default: no filtering]
  -h, --help                    Shows this help message

=head1 AUTHOR

 Angermueller Christof
 angermue@in.tum.de

=cut


### Variables ###


my @ali_files;
my @ali_neff;
my $outbase;
my $min_neff;
my $max_neff;
my $mode_neff;
my $filter;
my $tmp_file;


### Initialization ###


GetOptions( 
    "i|input-ali=s{1,}" => \@ali_files,
    "o|outbase=s"       => \$outbase,
    "m|min-neff=f"      => \$min_neff,
    "x|max-neff=f"      => \$max_neff,
    "d|mode-neff=f"     => \$mode_neff,
    "f|filter!"         => \$filter,
    "h|help"            => sub { pod2usage(2); }
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


my $i;  # Index of the selected alignment
if ($max_neff < $ali_neff[0]->[1]) { 
    if ($filter) {        
        (undef, $tmp_file) = tempfile(sprintf("%sXXXX", &filename($ali_neff[0]->[0])), SUFFIX => ".a3m", TMPDIR => 1);
        my $neff = &rand;
        `hhfilter -i $ali_neff[0]->[0] -o $tmp_file -neff $neff`;
        if (-s $tmp_file) { 
            unshift(@ali_neff, [$tmp_file, &neff($tmp_file)]); 
            $i = 0;
        }
    }
} elsif ($min_neff <= $ali_neff[$#ali_neff]->[1]) {
    my $neff = &rand;
    $i = 0;
    while ($i < $#ali_neff) {
        if (abs($neff - $ali_neff[$i + 1]->[1]) <= abs($neff - $ali_neff[$i]->[1])) { $i++; }
        else { last; }
    }
}


### Creating count profile ###


if (defined($i)) {
    printf("Creating count profile from %s...\n", basename($ali_neff[$i]->[0]));
    unless ($outbase) { $outbase = remove_suffix($ali_neff[$i]->[0], qr/\..+$/); }
    system("csbuild -i $ali_neff[$i]->[0] -o $outbase.prf -x 0");
    system(sprintf("cp $ali_neff[$i]->[0] $outbase%s", suffix($ali_neff[$i]->[0])));
    if ($tmp_file) { system("rm -f $tmp_file"); }
    printf("Count profile created with Neff %1.2f!\n", $ali_neff[$i]->[1]);
} else {
    print("There is no alignment whose Neff lies within the specified range!\n");
}





### Functions ###


sub neff {
    my ($file) = @_;
    `hhmake -i $file -o /dev/null` =~ /exp\(entropy\)\s*=\s*(\S+)/ or die "Error computing Neff for $file";
    return $1;
}

sub rand {
    if ($mode_neff) { return &rand_triangular($min_neff, $mode_neff, $max_neff); }
    else { return $min_neff + CORE::rand() * ($max_neff - $min_neff); }
    # Math::Random::OO::Normal->new(0.5 * ($max_neff + $min_neff), $sigma)->next;
}

# Draws a random number from the triangular distribution defined by $min, $max, $mode
sub rand_triangular {
    my ($min, $mode, $max) = @_;
    my $r = CORE::rand();
    if ($r <= ($mode - $min) / ($max - $min)) { return $min + sqrt($r * ($max - $min) * ($mode - $min)); }
    else { return $max - sqrt((1 - $r) * ($max - $min) * ($max - $mode)); }
}
