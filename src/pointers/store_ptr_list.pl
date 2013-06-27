#!/usr/bin/perl

use warnings;
use strict;

open(INPUT, "<ptrlist.in") or die $!;
open(TEMPLATE, "<ppm_t_ptr.in") or die $!;
open(OUTPUT, ">ppm_t_ptr.f") or die $!;

my @template = <TEMPLATE>;
close TEMPLATE;

for (<INPUT>) {
   chomp;
   my ($ptr, $type) = split " ";
   for (@template){
      my $line = $_;
      $line =~ s/DTYPE/$ptr/g;
      $line =~ s/CTYPE/$type/g;
      print OUTPUT $line ;
   }
}

close INPUT;
close OUTPUT;
