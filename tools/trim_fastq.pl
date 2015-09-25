#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;
use Switch;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: trim_fastq.pl N fasta [...]";
	exit;
}

my $length = scalar shift @ARGV;
die "need positive length" unless $length > 0;

my $id;
my $good = 1;
while (<>) {
  chomp;
  switch ($. % 4) {
    case 1 { 
      $id = $_;
      $good = 1;
    } case 2 {
      if (length $_ < $length) {
        $good = 0;
      } else {
        say $id;
        say substr $_, 0, $length;
      }
    } case 3 {
      say "+" if $good;
    } case 0 {
      say substr $_, 0, $length if $good;
    }
  }
}

__END__
