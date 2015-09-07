#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: check_intersect.pl forward_fastq reverse_fastq";
	exit;
}

my $for_fn = shift @ARGV;
my $rev_fn = shift @ARGV;

open(my $for_fh, $for_fn);
open(my $rev_fh, $rev_fn);

while (!eof($for_fh) and !eof($rev_fh)) {
  my $for_l = <$for_fh>;
  my $rev_l = <$rev_fh>;

  next unless $. % 4 == 1;

  # just grab the first part of line
  $for_l =~ /(^@[^\s\/]+)/;
  my $for_id = $1;
  $rev_l =~ /(^@[^\s\/]+)/;
  my $rev_id = $1;

  die "ids don't match: '$for_id' vs. '$rev_id' at fastq line $." unless $for_id eq $rev_id;
}

say "headers all match";

__END__
