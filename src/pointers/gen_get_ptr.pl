#!/usr/bin/perl
use warnings;
use strict;
my @types;
open(TYPEFILE, "<typelist.in") or die;
for (<TYPEFILE>) {
   chomp;
   push @types ,$_;
}
print @types;

open(TEMPLATE, '<pointer.in');
my @template = <TEMPLATE>;
close(TEMPLATE);
open(OUTPUT, '>ppm_pointers.f');

for my $class (@types) {
   for (@template){
      my $line = $_;
      $line =~ s/DTYPE/$class/g;
      print OUTPUT $line;
   }
}
open (PTRTYPES, "<ptrlist.in");
for my $ptr (<PTRTYPES>){
   chomp $ptr;
   push @types, $ptr;
   for (@template){
      s/CLASS\(DTYPE\)/TYPE\(DTYPE\), DIMENSION\(:\), POINTER/;
      s/DTYPE/$ptr/;
      print OUTPUT;
   }
}

close OUTPUT;
open(OUTPUT, '>ppm_pointers_interface.f');
print OUTPUT ' 'x9 . "INTERFACE get_pointer\n";
print OUTPUT ' 'x12 . "MODULE PROCEDURE ";

my @functions = map { "get_ptr_" . $_; } @types;
print OUTPUT join(", &\n" . ' 'x15, @functions) . "\n";

print OUTPUT ' 'x9 . "END INTERFACE get_pointer\n";
