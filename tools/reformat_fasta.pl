#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "concatenate and make each entry two lines\nusage: reformat_fasta.pl fasta [...]";
	exit;
}

my $seq = "";
while (<>) {
	chomp;
	if (/^>/) {
		say $seq unless $. == 1;
		$seq = "";
		say;
	} else {
		$seq .= $_;
	}
}

__END__