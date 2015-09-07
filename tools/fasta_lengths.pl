#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: fasta_lengths.pl fasta";
	exit;
}

my %tally = ();
my $len = 0;
readline;   # skip the first line
while (<>) {
    if (/^>/) {
        $tally{$len}++;
        $len = 0;
    } else {
        chomp;
        $len += length;
    }
}

my @lens = sort { $b <=> $a } keys %tally;
for my $l (@lens) { say "$l\t$tally{$l}"; }

__END__

get the distribution of fasta entry lengths
