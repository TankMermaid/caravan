#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: count_barcodes.pl fastq";
	exit;
}

my %tally = ();
while (<>) {
	next unless $. % 4 == 1;

	chomp;
	s{^.*#}{};
	s{\/[12]$}{};
	$tally{$_} += 1;
}

my @codes = sort { $tally{$b} <=> $tally{$a} || $a cmp $b } keys %tally;
for my $code (@codes) { say "$code\t$tally{$code}" };

__END__
