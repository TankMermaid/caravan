#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: extract_index_reads.pl fastq";
	exit;
}

while (<>) {
    next unless $. % 4 == 1;
    chomp;
    s/^@/>/;
    s/#/\n/;
    s{\/[12]$}{\n};
    print;
}

__END__

extract_index_reads.pl
author: Scott W Olesen (swo@mit.edu)

for a fastq with index reads in the @ line, put them in a new fasta.

e.g., a fastq might have a line @MISEQ:1:1111:18276:25622#TTACTATA/1
that encodes an index read TTACTATA, which this script spits out in
a new file with entries like

>MISEQ:1:1111:18276:25622
TTACTATA
