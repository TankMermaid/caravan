#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

my $in =  'ABCDGHMNRSTUVWXYabcdghmnrstuvwxy';
my $out = 'TVGHCDKNYSAABWXRtvghcdknysaabwxr';

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: reverse_complement.pl [seqs or fastas or lists of seqs]";
	exit;
}

if (-r $ARGV[0]) {
	while (<>) {
		if (/^>/) {
			chomp;
			say;
		} else {
			chomp;
			eval "tr/${in}/${out}/";
			say scalar reverse;
		}
	}
} else {
	while (my $seq = shift @ARGV) {
		eval "\$seq =~ tr/${in}/${out}/";
		$seq = reverse($seq);
		say $seq;
	}
}

__END__

e.g., you could call

rc.pl ACGT TATA GGTT
rc.pl newline_separated_seqs.txt
rc.pl foo.fasta