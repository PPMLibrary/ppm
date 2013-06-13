#!/usr/bin/perl

use strict;
use warnings;
my $spaces = ' 'x6;
my @tree_types;

# Print the tree type definitions
open TYPEDEF_I, "<tree_abstract_typedef.f.in" or die $!;
open TYPEDEF_O, ">tree_abstract_typedef.f" or die $!;
for (<TYPEDEF_I>){
   print TYPEDEF_O;
}
close TYPEDEF_I;

open PTRLIST, "<ptr_list.in" or die $!;
my @ptrlist = <PTRLIST>;
close PTRLIST;
for (@ptrlist) {
   chomp;
   my $tree_t = $_ . "_tree";
   print TYPEDEF_O $spaces . "TYPE, extends(intrinsic_tree) :: $tree_t\n";
   print TYPEDEF_O $spaces . "   TYPE($_), DIMENSION(:), POINTER :: val\n";
   print TYPEDEF_O $spaces . "END TYPE\n";
   push @tree_types, $tree_t;
}

open TYPELIST, "<typelist.in" or die $!;
my @typelist = <TYPELIST>;
close TYPELIST;
for (@typelist) {
   chomp;
   my $tree_t = $_ . "_tree";
   push @tree_types, $tree_t;
   print TYPEDEF_O $spaces . "TYPE, extends(derived_tree) :: $tree_t\n";
   print TYPEDEF_O $spaces . "   TYPE($_), DIMENSION(:), POINTER :: val\n";
   print TYPEDEF_O $spaces . "END TYPE\n";
}

close TYPEDEF_O;

# Add the types to the check associated function
open ABSTRACT_I, "<tree_abstract.f.in" or die $!;
open ABSTRACT_O, ">tree_abstract.f" or die $!;

$spaces .= ' 'x6;
my $select_line = "";
for (@tree_types) {
   $select_line .= $spaces . "   TYPE is ($_)\n";
   $select_line .= $spaces . "      select TYPE (newnode)\n";
   $select_line .= $spaces . "      TYPE is ($_)\n";
   $select_line .= $spaces . "         check_associated = associated(treenode%val, &\n";
   $select_line .= $spaces . "            TARGET=newnode%val)\n";
   $select_line .= $spaces . "      END SELECT\n";
}

for (<ABSTRACT_I>) {
   s/^ *! TYPE_SELECT/$select_line/;
   print ABSTRACT_O;
}
