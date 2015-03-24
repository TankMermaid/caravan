#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

while (<>) {
    chomp;
    my @fields = split /;/;
    my $id = shift @fields;
    my $size = shift @fields;
    die unless $size =~ /size=/;

    die unless shift @fields eq "";
    next unless shift @fields eq "+";
    die unless shift @fields eq "Root";
    die "not 100% for $_" unless shift @fields eq "100%";

    my @out = ();
    while (@fields) {
        my $tax = shift @fields;
        my $conf = shift @fields;
        $conf =~ s/\%//;
        last if $conf < 50;
        push @out, "$tax($conf\%)";
    }
    say "$id;$size;\t" . (join ",", @out);
}


__END__

convert an RDP allrank file to a utax format