#!/usr/bin/perl

use warnings;
use strict;

open(INPUT, "<ptrlist.in") or die $!;
open(TEMPLATE, "<ppm_t_ptr.in") or die $!;
open(OUTPUT, ">ppm_t_ptr.f") or die $!;

my @template = <TEMPLATE>;
close TEMPLATE;

for my $ptr (<INPUT>) {
   chomp $ptr;
   for my $line (@template){
      $line =~ s/DTYPE/$ptr/g;
      print OUTPUT $line;
   }
}

close INPUT;
close OUTPUT;
