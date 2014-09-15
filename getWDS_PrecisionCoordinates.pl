#! /usr/bin/perl
# Grab only the precise coordinates from the WDS entires and print them out.
# Ignore entires without these coordinates.

use strict;

unless ((defined($ARGV[0])) and (-f $ARGV[0])) {
  die "Usage: rearrangeWDS.pl <WDS File to fix>\n";
}

open FILE, $ARGV[0] or die "Can't open $ARGV[0].\n";
open NEW, ">WDS.re" or die "Can't open WDS.re.\n";

while(<FILE>) {
  my $s = substr($_, 112);
  if ($s =~ /\d+/) { print NEW $s; }
}
close FILE;
close NEW;
