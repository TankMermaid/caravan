#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: fastq_lengths.pl fastq";
	exit;
}

my %tally = ();
while (<>) {
	next unless $. % 4 == 2;

	chomp;
	$tally{(length)}++;
}

my @lens = sort { $b <=> $a } keys %tally;
for my $l (@lens) { say "$l\t$tally{$l}"; }

__END__

get the distribution of read lengths
