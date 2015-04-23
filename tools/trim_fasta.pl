#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "trim sequences \nusage: trim_fasta.pl N fasta [...]";
	exit;
}

my $length = scalar shift @ARGV;
die "need positive length" unless $length > 0;

my $line = <>;
chomp $line;
die "not a fasta" unless (substr $line, 0, 1) eq ">";
say $line;

my $seq = "";
while (<>) {
	chomp;
	if (/^>/) {
		say (substr $seq, 0, $length);
		$seq = "";
		say;
	} else {
		$seq .= $_;
	}
}
say (substr $seq, 0, $length);

__END__
