#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
  say "usage: tsv2fasta.pl [files]";
  exit;
}

while (<>) {
  my @fields = split;
  say ">" . $fields[0];
  say $fields[1];
}


__END__

Convert a tab-separated file (two fields per line) into a fasta
file: first line is header, second line is sequence
