#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

my $i = 1;

while (<>) {
  my $pos = $. % 4;
  if ($pos == 1) {
    say "\@read$i";
    $i++;
  } elsif ($pos == 3) {
    say "+";
  } else {
    print;
  }
}

__END__
