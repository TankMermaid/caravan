#!/usr/bin/env perl

use warnings;
use strict;
use 5.10.1;
use Getopt::Long;
use Pod::Usage;

my $max_entries = 0;
my $verbose = 0;
my $help = 0;
my $rc = 0;

GetOptions("n=i" => \$max_entries
    , 'verbose' => \$verbose
    , 'rc' => \$rc
    , 'help' => \$help)
    or pod2usage(2);

pod2usage(2) if $help;

my $primer = shift @ARGV;
$primer = uc $primer;

if ($rc) {
    $primer =~ tr/ATGCUNYRSWKMBDHV/TACGANRYSWMKVHDB/;
    $primer = reverse($primer);
}

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
my $max_position = 0;
my $min_position = -1;
while (<>) {
	next unless $. % 4 == 2;
	if (/${regex}/) {
      $average_position = ($average_position * $matches + $-[0]) / ($matches + 1);

      $max_position = $-[0] if $-[0] > $max_position;
      $min_position = $-[0] if ($min_position == -1 or $-[0] < $min_position);

      $matches += 1;
    }
	$entries += 1;
	last if $max_entries > 0 and $entries >= $max_entries;
}

my $match_percent = sprintf("%i", $matches / $entries * 100) . "%";
say "of $entries entries, $matches match the primer (i.e., $match_percent)";
if ($matches > 0) {
    say "min/avg/max position: " . sprintf("%i / %.2f / %i", $min_position, $average_position, $max_position);
}

__END__

=head1 NAME

check_primers

=head1 SYNOPSIS

check_primers.pl PRIMER [options] [file]

    Options
      -help
      -n N   check only the first N entries
      -rc    take reverse complement of provided primer

=head1 DESCRIPTION

This program will read the input fastqs and see how many have the
given primer in the given window. It does regex matching for non-ACGT
characters.

=cut
