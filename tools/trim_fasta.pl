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

my $name = "";
my $seq = "";
while (<>) {
	chomp;
	if (/^>/) {
		output_entry($name, $seq) unless $. == 1;
		$name = $_;
		$seq = "";
	} else {
		$seq .= $_;
	}
}
output_entry($name, $seq);

sub output_entry {
	if ((length $seq) > $length) {
		say $name;
		say (substr $seq, 0, $length);
	}
}

__END__
