#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

my $i = 1;

while (<>) {
  if ($. % 4 == 1) {
    say "\@read$i";
    $i++;
  } else {
    print;
  }
}

__END__
