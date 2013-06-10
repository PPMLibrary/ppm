#!/usr/bin/perl

use warnings;
use strict;

my $root;
my ($csize, $isize, $dsize);

open(TYPEDEF, "<", $ARGV[0]) or die;
my @typedef = reverse <TYPEDEF>;

my @typestack;
sub spaces {
   my $base_spaces=" "x15 . " "x(3*($#typestack));
   return $base_spaces;
}


sub parse_type {
   my %ints;
   my %dubs;
   my %chars;
   my %bool;
   my %pointers;
   my @eval;
   # this is the main evaluation loop
   while (@typedef) {
      $_ = pop @typedef;
      #print "Eval $_\n";
      if (m/ *#/) {
         # Ignore comment lines
         next;
      }
      if (m/integer (\w*) *(\w*)?/) {
         if ($2) { $ints{$1} = $2; } else { $ints{$1} = '1'; }
      }
      elsif (m/real (\w*) *(\w*)?/){
         if ($2) { $dubs{$1} = $2; } else { $dubs{$1} = '1'; }
      }
      elsif (m/character *(\w*) *(\w*)?/){
         if ($2) { $chars{$1} = $2; } else { $chars{$1} = '1'; }
      }
      elsif (m/logical *(\w*) *(\w*)?/){
         if ($2) { $bool{$1} = $2; } else { $bool{$1} = '1'; }
      }
      elsif (m/pointer *(\w*) *(\w*)/){
         print "Parse pointer\n";
         $pointers{$1} = $2;
      }
      elsif (m/\btype *(\w+)/){
         my $type = $1;
         #print "pushing $1\n";
         push @typestack, $type;
         my @subtype = &parse_type();
         $type = pop @typestack;
         #print "pop $type\n";
         if (@typestack){
            for (@subtype) {
               if ($_){
               $_  = join("\n   ", split('\n', $_)) . "\n";
               $_ = "\n" . &spaces() ."SELECT TYPE(type_ptr)\n" .&spaces().
                      "   CLASS IS($type)\n" .
                      $_ . &spaces(). "END SELECT\n";
                   }
            }
         }
         else {
            $root = $type;
         }
         for (0..2){
            $eval[$_] .= $subtype[$_];
         }
      }
      elsif (m/\bend\b/) {
         print "End TYPE \n";
         last;
      }
      else {
         die  "Failed to evaluate expr: ".$_;
      }
   }
   if (@typestack){
      my @ret = &eval_type(\%ints,\%dubs,\%chars,\%bool,\%pointers);
      for (0..2) {
         if ($eval[$_]){ $eval[$_] = $ret[$_] . $eval[$_]; }
         else { $eval[$_] = $ret[$_]; }
      }
   }
   #print $eval[0];
   return ($eval[0],$eval[1], $eval[2]);
}

# Evaluates the code for the type.
sub eval_type {
   #my (%ints, %dubs, %chars, %bool, %pointers) = @_;
   my $calc = &eval_calc(@_);
   my $create = &eval_create(@_);
   my $write = &eval_write(@_);
   #return "Type evaluated\n";
   #return $create;
   return ($calc, $create, $write);
}
sub eval_calc {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $calc_section .= "\n";
   my $tsizeline = &spaces() . "tsize = tsize";

   # Well flag this for later
   $calc_section .= &spaces() . "! Calculate datatype size\n";
   if (%$ints and not $isize) {
      $isize = "CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)\n";
   }
   if (%$dubs and not $dsize) {
      $dsize = "CALL h5tget_size_f(H5T_NATIVE_DOUBLE, dsize, error)\n";
   }
   if ((%$bool or %$chars or %$pointers) and not $csize) {
      $csize = "CALL h5tget_size_f(H5T_NATIVE_CHARACTER, csize, error)\n";
   }
   if ($#typestack == 0){
      if ($isize) { $calc_section .= &spaces() . $isize; }
      if ($dsize) { $calc_section .= &spaces() . $dsize; }
      if ($csize) { $calc_section .= &spaces() . $csize; }
   }

   $calc_section .= "\n";
   $tsizeline .= &sum_type($ints, 'isize');
   $tsizeline .= &sum_type($dubs, 'dsize');
   $tsizeline .= &sum_type($bool, 'csize');
   $tsizeline .= &sum_type($chars, 'csize');
   if (%$pointers) {
      $tsizeline .= " + csize*32*". keys(%$pointers);
   }

   $calc_section .= "\n";
   $calc_section .= $tsizeline . "\n\n";
   return $calc_section;
}
sub sum_type {
   my ($hash, $var) = @_;
   my $count = 0;
   my $cline = $var . "*(";
   my $tsizeline = "";
   for my $key (keys %$hash) {
      if ($$hash{$key} =~ /\d/){
         $count += $$hash{$key};
      }
      else {
         $cline .= $$hash{$key} . "+";
      }
      #$calc_section .= $key . ' ' . $$hash{$key};
   }
   if ($count > 0 or not $cline eq $var . "*("){
      $tsizeline .= " + ". $cline . $count . ")";
   }
   return $tsizeline;
}
sub print_attr {
   my ($type, $names, $vn) = @_;
   my $create_section;
   for my $var (keys %$names){
      if ($type eq "pointer"){
         $create_section .= print_ptr($var);
      }
      elsif ($$names{$var} eq '1'){
         $create_section .= print_el ($var, $type, $vn);
      }
      else {
         $create_section .= print_arr($type, $var, $vn."*".$$names{$var});
      }
   }
   $create_section .= "\n";
   return $create_section;
}
sub print_el {
   my $type = $_[1];
   my $create_section = "";
   if ($type eq "logical"){
      $type = "H5T_NATIVE_CHARACTER"
   }
   $create_section .= &spaces();
   $create_section .= "CALL h5tinsert_f(dtype_id, \"".$_[0]."\", offset, &\n";
   $create_section .= &spaces() . "      " . $type .", error)\n";
   $create_section .= &spaces() . "offset = offset + ". $_[2] . "\n";
   return $create_section;
}
sub print_arr {
   my ($type, $var, $num) = @_;
   my $create_section = "";
   $create_section .= &spaces() . "dims = (/$num/)\n";
   if ($type eq "H5T_NATIVE_CHARACTER") {
      $create_section .= &spaces() . "CALL h5tcreate_f(H5T_STRING_F, $num, &\n";
      $create_section .= &spaces() . "    array_id, error)\n";
   }
   else {
      $type =~ s/logical/H5T_NATIVE_CHARACTER/;
      $create_section .= &spaces() . "CALL h5tarray_create_f(".$type.", rank, &\n";
      $create_section .= &spaces() . "    dims, array_id, error)\n";
   }
   $create_section .= &spaces() . "CALL h5tinsert_f(dtype_id, \"".$var."\", offset, &\n";
   $create_section .= &spaces() . "    array_id, error)\n";
   $create_section .= &spaces() . "offset = offset + (". $num. ")\n";
   #$numarr += 1;
   return $create_section;
}
sub print_ptr {
   my $name = $_[0];
   my $create_section = "";
   $create_section .= &spaces() . "CALL h5tcreate_f(H5T_STRING_F, 32*csize, &\n";
   $create_section .= &spaces() . "    array_id, error)\n";
   $create_section .= &spaces() . "CALL h5tinsert_f(dtype_id, \"".$name."\", offset, &\n";
   $create_section .= &spaces() . "    array_id, error)\n";
   $create_section .= &spaces() . "offset = offset + (32*csize)\n";
   return $create_section;
}
sub eval_create {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $create_section;
   $create_section .= &spaces() . "! Insert the members\n";

   $create_section .= &spaces() . "! Integer members\n";
   $create_section .= print_attr("H5T_NATIVE_INTEGER", $ints, "isize");

   $create_section .= &spaces() . "! Real members\n";
   $create_section .= print_attr("H5T_NATIVE_DOUBLE", $dubs, "dsize");

   $create_section .= &spaces() . "! Character members\n";
   $create_section .= print_attr("H5T_NATIVE_CHARACTER", $chars, "csize");

   $create_section .= &spaces() . "! Logical members\n";
   $create_section .= print_attr("logical", $bool, "csize");

   $create_section .= &spaces() . "! Pointer members\n";
   $create_section .= print_attr("pointer", $pointers, "csize*32");
   return $create_section;
}

sub eval_write {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $write_section = "";
   # Now we generate the code for the write function
   for my $map ($ints, $dubs, $chars, $bool) {
      for my $key (keys %$map) {
         $write_section .= &spaces() . "CALL write_attribute(dset_id, \'$key\', &\n";
         if ($$map{$key} eq '1'){
            $write_section .= &spaces() . "    type_ptr%$key)\n";
         }
         else {
            $write_section .= &spaces() . "    type_ptr%$key, $$map{$key})\n";
         }
      }
   }
   for my $ptr (keys %$pointers) {
      $write_section .= &spaces() . "IF (associated(type_ptr%$ptr)) THEN\n";
      $write_section .= &spaces() . "   pointer_addr = get_pointer(type_ptr%$ptr)\n";
      #$write_section .= &spaces() . "   CALL store_$$pointers{$ptr}(cpfile_id, &\n";
      $write_section .= &spaces() . "   CALL store_type(cpfile_id, &\n";
      $write_section .= &spaces() . "       pointer_addr, type_ptr%$ptr)\n";
      $write_section .= &spaces() . "ELSE\n";
      $write_section .= &spaces() . "   pointer_addr = \"00000000000000000000000000000000\"\n";
      $write_section .= &spaces() . "ENDIF\n";
      $write_section .= &spaces() . "CALL write_attribute(dset_id, \"$ptr\", pointer_addr, 32)\n";
   }
   return $write_section;
}
my ($calc, $create, $write) = &parse_type();
print $calc;
#print $create;
#print $write;
#print $root;

my $filename = $ARGV[1];
my $parenttype;
open TEMPLATE, "+<", "type.f.in" or die $!;
open OUTFILE, ">", $filename or die $!;
print OUTFILE &spaces() . "! File autogenerated by my super script\n";

for (<TEMPLATE>) {

   # if its not an omitted line, print
   if (not ((m/isize/ and (not $isize))
            or (m/dsize/ and (not $dsize))
            or (m/csize/ and (not $csize))
         #or (m/dims/ and ($numarr == 0))
         #or ((m/array_id/ or m/rank/) and ($numarr == 0 and not %pointers))
         # or (m/pointer_addr/ and (not %pointers))
         #or (m/subsize/ and ($numarr == 0))
          or ((/[sg]et_size/ or /PARENT/i) and (not $parenttype)))) {

      # Make the template substitutions
      $_ =~ s/DTYPE/$root/;
      $_ =~ s/ *!CREATE_STUB/$create/;
      $_ =~ s/ *!CALCULATE_STUB/$calc/;
      s/ *!WRITE_STUB/$write/;
      #s/ *!READ_STUB/$read_section/;
      print OUTFILE;
   }
}
close TEMPLATE;
close OUTFILE;
