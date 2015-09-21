#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: fastq2fasta.pl fastq";
	exit;
}

while (<>) {
  if ($. % 4 == 1) {
    s/^@/>/;
    print;
  } elsif ($. % 4 == 2) {
    print;
  }
}

__END__

fastq2fasta.pl
author: Scott W Olesen (swo@mit.edu)

change a fastq into fasta format
