#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

my $i = 1;

while (<>) {
  if ($. % 2 == 1) {
    die "need two-line format" unless /^>/;
    say ">read$i";
    $i++;
  } else {
    die "need two-line format" unless not /^>/;
    print;
  }
}

__END__
