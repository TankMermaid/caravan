#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;

if (grep($_ eq "-h" || $_ eq "--help", @ARGV)) {
	say "usage: check_primers.pl primer fastq";
	exit;
}

my $primer = shift @ARGV;

# prepare a regex library
my %code = qw{W [AT] S [CG] M [AC] K [GT] R [AG] Y [CT]
	B [CGT] D [AGT] H [ACT] V [ACG] N [ACGT] };
my $check = join '|', keys %code;

# grab the primer and convert to a sensible regex
my $regex = $primer
$regex =~ s/($check)/$code{$1}/g;

my $entries = 0;
my $matches = 0;
while (<>) {
	next unless $. % 4 == 2;
	$matches += 1 if /${primer}/;
	$entries += 1;
}

say "of $entries entries, $matches match the primer $primer (aka regex $regex)"

__END__
