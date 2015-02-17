#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: count_by_otu.pl otu_table";
	exit;
}

while (<>) {
    # ignore the header line
    next if $. == 1;

    # split into fields and drop the otu name
    my @fields = split;
    shift @fields;

    my $sum = 0;
    grep { $sum += $_ } @fields;
    print "$sum\n";
}

__END__
