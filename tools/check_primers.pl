#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;
use Getopt::Long;
use Pod::Usage;

my $max_entries = 0;
my $verbose = 0;
my $help = 0;

GetOptions("n=i" => \$max_entries
    , 'verbose' => \$verbose
    , 'help' => \$help)
    or pod2usage(2);

pod2usage(2) if $help;

my $primer = shift @ARGV;

# prepare a regex library
my %code = qw{W [AT] S [CG] M [AC] K [GT] R [AG] Y [CT]
	B [CGT] D [AGT] H [ACT] V [ACG] N [ACGT] };
my $check = join '|', keys %code;

# grab the primer and convert to a sensible regex
my $regex = $primer;
$regex =~ s/($check)/$code{$1}/g;

say "searching for primer $primer => regex $regex..." if $verbose;

my $entries = 0;
my $matches = 0;
my $average_position = 0.0;
while (<>) {
	next unless $. % 4 == 2;
	if (/${regex}/) {
      $average_position = ($average_position * $matches + @-) / ($matches + 1);
      $matches += 1;
    }
	$entries += 1;
	last if $max_entries > 0 and $entries >= $max_entries;
}

say "of $entries entries, $matches match the primer (at average position $average_position)"

__END__

=head1 NAME

check_primers

=head1 SYNOPSIS

check_primers.pl PRIMER [options] [file]

    Options
      -help
      -n    number of entries

=head1 DESCRIPTION

This program will read the input fastqs and see how many have the
given primer in the given window. Will do regex matching for non-ACGT
characters.

=cut
