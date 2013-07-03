#!/usr/bin/perl

use warnings;
use strict;

open(INPUT, "<ptrlist.in") or die $!;
open(TEMPLATE, "<ppm_t_ptr.f.in") or die $!;
open(OUTPUT, ">ppm_t_ptr.f") or die $!;

my @template = <TEMPLATE>;
close TEMPLATE;

my $debug = 3;
sub check_log_level {
   my $line = $_;
   if (( ($line =~ /!CRITICAL/) and $debug < 1) or
       ( ($line =~ /!ERROR/) and $debug < 2) or
       ( ($line =~ /!WARNING/) and $debug < 3) or
       ( ($line =~ /!INFO/) and $debug < 4) or
       ( ($line =~ /!DEBUG/) and $debug < 5)  ) {
       return 1;
   }
   return 0;
}

for (<INPUT>) {
   chomp;
   my ($ptr, $type) = split " ";
   for (@template){
      my $line = $_;
      $line =~ s/WRITE/!WRITE/ if check_log_level() != 0;
      $line =~ s/DTYPE/$ptr/g;
      $line =~ s/CTYPE/$type/g;
      print OUTPUT $line ;
   }
}

close INPUT;
close OUTPUT;
