#!/usr/bin/perl -w
use strict;

open FILE, "<@ARGV" or die $!;
open OUT, ">@ARGV.proportion" or die $!;

my @file=<FILE>;
close FILE;

my $file;

foreach $file(@file) {
    if ($file =~ m/^>/) {
    }
    else{
        print OUT $file;
        }
}
close OUT;

open OUT, "<@ARGV.proportion" or die $!;

my ($chars, $X_chars) = (0,0);



while (<OUT>) {
      $chars += length;
      $X_chars += ($_ =~ tr/X//);
}


my $proportion=(( $X_chars / $chars)*100);

open SUMMARY, ">@ARGV.sequencing.summary" or die $!;

print SUMMARY "ORF statistics for:";
print SUMMARY @ARGV;
print SUMMARY "\n\nTotal characters:";
print SUMMARY $chars;
print SUMMARY "\n\nTotal X's:";
print SUMMARY $X_chars;
print SUMMARY "\n\nProprotion X's:";
print SUMMARY $proportion;
close SUMMARY;
close OUT;

