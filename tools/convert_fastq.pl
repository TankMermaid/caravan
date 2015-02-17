#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: convert_fastq.pl fastq";
	exit;
}

while (<>) {
    $_ = "+\n" if $. % 4 == 3;
    tr/A-h/"-I/ if $. % 4 == 0;
    print;
}

__END__

convert_fastq.pl
author: Scott W Olesen (swo@mit.edu)

convert a fastq with Illumina 1.3-style Phred scores (A-h) to Illumina 1.8-style
Phred score ("-I). remove all information from the + lines.