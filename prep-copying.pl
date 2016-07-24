#!/usr/bin/env perl

## prepand copying information

$#ARGV >= 1 || die "usage: prep-copying.pl <notice> <file 1> [<file 2> ...]\n";

$copying = shift @ARGV;

print "$copying\n";

foreach (@ARGV) {
    print "$_\n";
}
