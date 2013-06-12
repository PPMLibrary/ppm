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
      my $line = $_;
      $line =~ s/CLASS\(DTYPE\)/TYPE\(DTYPE\), DIMENSION\(:\), POINTER/;
      $line =~ s/DTYPE/$ptr/;
      print OUTPUT $line;
   }
}

close OUTPUT;
open(OUTPUT, '>ppm_pointers_interface.f');
print OUTPUT ' 'x9 . "INTERFACE get_pointer\n";
print OUTPUT ' 'x12 . "MODULE PROCEDURE ";

my @functions = map { "get_ptr_" . $_; } @types;
print OUTPUT join(", &\n" . ' 'x15, @functions);

print OUTPUT ", &\n";
print OUTPUT ' 'x15 . "get_real1d_pointer, &\n";
print OUTPUT ' 'x15 . "get_real2d_pointer, &\n";
print OUTPUT ' 'x15 . "get_integer2d_pointer, &\n";
print OUTPUT ' 'x15 . "get_integer1d_pointer\n";

print OUTPUT ' 'x9 . "END INTERFACE get_pointer\n";
