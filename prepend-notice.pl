#!/usr/bin/env perl

## prepand copying information

$#ARGV >= 1 || die "usage: prepend-notice.pl <notice> <file 1> [<file 2> ...]\n";

$copying = shift @ARGV;

print "#notice: $copying\n";

foreach (@ARGV) {
    $orig = `cat $_`;
    open(FOUT, ">$_");
    print "$copying\n";
    print "$orig";
    close(FOUT);
}
