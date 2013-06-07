#!/usr/bin/perl

use warnings;
use strict;

my $root;

open(TYPEDEF, "<", $ARGV[0]) or die;
my @typedef = reverse <TYPEDEF>;

my @typestack;

sub eval_type {
   $_ = $_[0];
   my $eval_line = "";
   if (m/\btype *(\w+)/){
      my $type = $1;
      print "pushing $1\n";
      push @typestack, $type;
      $eval_line = &parse_type();
      return ("type", $eval_line);
   }
   elsif (m/\bend\b/) {
      return ('end', '');
   }
   else {
      print "Eval $_\n";
      $eval_line = $_;
   }
   return ('', $eval_line);
}

sub parse_type {
   my $line = pop @typedef;
   my ($ret,$eval) = &eval_type($line);
   if ($ret eq 'type') {
      my $type = pop @typestack;
      print "End TYPE $type\n";
      $eval .= "SELECT TYPE(type_ptr)\n";
      $eval .= "   CLASS IS($type)\n" . "END SELECT\n";
      $eval =~ s/(\z.)/$1   /g;
   }
   elsif ($ret eq 'end'){
      return "";
   }
   if (@typedef) {
      $eval .= &parse_type();
   }
   return $eval;
}
my $code = &parse_type();
print $code;
