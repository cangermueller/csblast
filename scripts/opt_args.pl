#!/usr/bin/perl -w

use YAML::Tiny;


$_ = YAML::Tiny->read($ARGV[0]) or exit 1;
my @args;
foreach my $arg (keys(%{$_->[0]})) {
  push(@args, sprintf("-%s <%%= %s %%>", $arg, $arg));
}
print join(" ", @args), "\n";

