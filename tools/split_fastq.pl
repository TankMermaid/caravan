#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: split_fastq.pl size(Gb) fastq";
	exit;
}

my $size = shift @ARGV;
my $in_fn = shift @ARGV;
my $suffix = 0;

$size *= 2**30;

open(my $in, $in_fn);

my $out_fn = "${in_fn}.${suffix}";
open(my $out, '>', $out_fn) or die "can't open $out_fn";

while (<$in>) {
	print $out $_;

	if (($. % 4 == 0) && (tell($out) > $size)) {
		close($out);
		$suffix++;
		my $out_fn = "${in_fn}.${suffix}";
		open($out, '>', $out_fn) or die "can't open $out_fn";
	}
}

__END__

split_fastq.pl
author: Scott W Olesen (swo@mit.edu)

split a fastq into separate files based on the size (in Gb) of the desired
output files.