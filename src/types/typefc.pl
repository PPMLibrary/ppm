#!/usr/bin/perl

package TypeCompiler;
use warnings;
use strict;

my $root;
my ($csize, $isize, $dsize);
my @allocatables;
my @derived_types;
my $arrays = 0;
my $numpointers = 0;

open(TYPEDEF, "<", $ARGV[0]) or die;
my @typedef = reverse <TYPEDEF>;

my @typestack;
sub spaces {
   my $base_spaces=" "x15 . " "x(3*($#typestack));
   return $base_spaces;
}


sub parse_type {
   my ($type) = @_;
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
      if (m/integer +(\w+) *(\w*)?/) {
         if ($2) { $ints{$1} = $2; $arrays++; } else { $ints{$1} = '1'; }
      }
      elsif (m/real +(\w+) *(\w*)?/){
         if ($2) { $dubs{$1} = $2; $arrays++; } else { $dubs{$1} = '1'; }
      }
      elsif (m/character +(\w+) *(\w*)?/){
         if ($2) { $chars{$1} = $2; $arrays++; } else { $chars{$1} = '1'; }
      }
      elsif (m/logical +(\w+) *(\w*)?/){
         if ($2) { $bool{$1} = $2; $arrays++; } else { $bool{$1} = '1'; }
      }
      elsif (m/pointer +(\w+) *(\w*)?/){
         $pointers{$1} = $2;
         $numpointers++;
      }
      elsif (/allocatable/) {
         print "allocatable list\n";
         while (@typedef) {
            $_ = pop @typedef;
            chomp;
            last if m/\bend\b/;
            m/(\w+)/;
            push @allocatables, $1;
         }
      }
      elsif (m/\btype +(\w+)/){
         my $type = $1;
         #print "pushing $1\n";
         push @typestack, $type;
         push @derived_types, $type;
         my @subtype = &parse_type($type);
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
         for (0..3){
            $eval[$_] .= $subtype[$_];
         }
      }
      elsif (m/\bend\b/) {
         last;
      }
      else {
         die  "Failed to evaluate expr: ".$_;
      }
   }
   if (@typestack){
      my @ret = &eval_type(\%ints,\%dubs,\%chars,\%bool,\%pointers);
      for (0..3) {
         if ($eval[$_]){ $eval[$_] = $ret[$_] . $eval[$_]; }
         else { $eval[$_] = $ret[$_]; }
      }
   }
   #print $eval[0];
   return ($eval[0],$eval[1], $eval[2],$eval[3]);
}

# Evaluates the code for the type.
sub eval_type {
   #my (%ints, %dubs, %chars, %bool, %pointers) = @_;
   my $calc = &eval_calc(@_);
   my $create = &eval_create(@_);
   my $write = &eval_write(@_);
   my $recover = &eval_recover(@_);
   #return "Type evaluated\n";
   #return $create;
   return ($calc, $create, $write, $recover);
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
   $create_section .= spaces() . "CALL h5tcreate_f(H5T_STRING_F, 32*csize, &\n";
   $create_section .= spaces() . "    array_id, error)\n";
   $create_section .= spaces() . "CALL h5tinsert_f(dtype_id, \"".$name."\", offset, &\n";
   $create_section .= spaces() . "    array_id, error)\n";
   $create_section .= spaces() . "offset = offset + (32*csize)\n";
   return $create_section;
}
sub eval_create {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $create_section = "";

   $create_section .= spaces() . "! Integer members\n";
   $create_section .= print_attr("H5T_NATIVE_INTEGER", $ints, "isize");

   $create_section .= spaces() . "! Real members\n";
   $create_section .= print_attr("H5T_NATIVE_DOUBLE", $dubs, "dsize");

   $create_section .= spaces() . "! Character members\n";
   $create_section .= print_attr("H5T_NATIVE_CHARACTER", $chars, "csize");

   $create_section .= spaces() . "! Logical members\n";
   $create_section .= print_attr("logical", $bool, "csize");

   $create_section .= spaces() . "! Pointer members\n";
   $create_section .= print_attr("pointer", $pointers, "csize*32");
   return $create_section;
}

sub eval_write {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $write_section = "";
   # Now we generate the code for the write function
   for my $map ($ints, $dubs, $chars, $bool) {
      for my $key (keys %$map) {
         # attribute level debug statements
         # $write_section .= &spaces() . "WRITE (*,*) \"$key\"\n"; # dbg
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
      #$write_section .= &spaces() . "WRITE (*,*) \"$ptr\"\n"; # dbg
      $write_section .= &spaces() . "IF (associated(type_ptr%$ptr)) THEN\n";
      $write_section .= &spaces() . "   pointer_addr = get_pointer(type_ptr%$ptr)\n";
      #$write_section .= &spaces() . "   call write_tree(pointer_data%itree)\n"; # dbg
      #$write_section .= &spaces() . "   CALL store_$$pointers{$ptr}(cpfile_id, &\n";
      $write_section .= &spaces() . "   CALL store_type(cpfile_id, &\n";
      $write_section .= &spaces() . "       pointer_addr, type_ptr%$ptr)\n";
      $write_section .= &spaces() . "ELSE\n";
      $write_section .= &spaces() . "   pointer_addr = \"00000000000000000000000000000000\"\n";
      # $write_section .= &spaces() . "   WRITE (*,*) \"   is null\"\n"; # dbg
      $write_section .= &spaces() . "ENDIF\n";
      $write_section .= &spaces() . "CALL write_attribute(dset_id, \"$ptr\", pointer_addr, 32)\n";
   }
   return $write_section;
}
sub eval_recover {
   my ($ints, $dubs, $chars, $bool, $pointers) = @_;
   my $recover_section = "";
   # Now we generate the code for the write function
   for my $map ($ints, $dubs, $chars, $bool) {
      for my $key (keys %$map) {
         # attribute level debug statements
         # $write_section .= &spaces() . "WRITE (*,*) \"$key\"\n";
         $recover_section .= &spaces() . "CALL read_attribute(dset_id, \'$key\', &\n";
         if ($$map{$key} eq '1'){
            $recover_section .= &spaces() . "    type_ptr%$key)\n";
         }
         else {
            $recover_section .= &spaces() . "    type_ptr%$key, $$map{$key})\n";
         }
      }
   }
   for my $ptr (keys %$pointers) {
      #$write_section .= &spaces() . "WRITE (*,*) \"$ptr\"\n"; # dbg
      $recover_section .= &spaces() . "CALL read_attribute(dset_id, \"$ptr\", pointer_addr, 32)\n";
      $recover_section .= &spaces() . "WRITE (*,*) \"pointer addr is \", pointer_addr\n";
      $recover_section .= &spaces() . "IF (pointer_addr == \"00000000000000000000000000000000\") THEN\n";
      $recover_section .= &spaces() . "   type_ptr%$ptr => null()\n";
      $recover_section .= &spaces() . "ELSE\n";
      #$recover_section .= &spaces() . "   type_ptr%$ptr => recover_type(cpfile_id, pointer_addr, type_ptr%$ptr)\n";
      $recover_section .= &spaces() . "   type_ptr%$ptr => recover_$pointers->{$ptr}(cpfile_id, pointer_addr, type_ptr%$ptr)\n";
      $recover_section .= &spaces() . "   IF (.NOT. associated(type_ptr%$ptr)) THEN\n";
      $recover_section .= &spaces() . "      write(*,*) \"Failed recovering $ptr\"\n";
      $recover_section .= &spaces() . "   ENDIF\n";
      $recover_section .= &spaces() . "ENDIF\n";
   }
   return $recover_section;
}
my ($calc, $create, $write,$recover) = &parse_type();

# Template for the function to recover a specific derived type
# This is necessary since you must return bottom up when returning a pointer
# in order to prevent losing information
# This way we allocate the proper type in the right pointer type
my $recover_function_template =<< "END_TEMPLATE";
            CLASS(CTYPE) FUNCTION recover_CTYPE(cpfile_id, type_ptr_id, type_ptr)
               POINTER :: recover_CTYPE
               CLASS(CTYPE), POINTER :: type_ptr
               CLASS(DTYPE), POINTER :: abstr_type_ptr
               CLASS(derived_tree), POINTER :: tree_node_ptr
               INTEGER(HID_T) :: cpfile_id, dset_id, group_id, &
                  recover_type_id
               CHARACTER(LEN=*) :: type_ptr_id
               INTEGER(HSIZE_T) :: offset = 0, isize
               INTEGER(HSIZE_T), DIMENSION(1) :: dims = (/0/)
               INTEGER :: type_num, error
               recover_CTYPE => type_ptr
               CALL lookup_pointer(type_ptr_id, pointer_data%dtree, tree_node_ptr)
               IF (associated(tree_node_ptr)) THEN
                  SELECT TYPE(tree_node_ptr)
                  CLASS is (DTYPE_tree)
                     abstr_type_ptr => tree_node_ptr%val
                  END SELECT
               ELSE
                  CALL h5gopen_f(cpfile_id, 'DTYPE', &
                     group_id, error)
                  CALL h5dopen_f(group_id, type_ptr_id, dset_id, &
                      error)

                  CALL h5tget_size_f(H5T_NATIVE_INTEGER, isize, error)
                  CALL h5tcreate_f(H5T_COMPOUND_F, isize, recover_type_id, error)
                  CALL h5tinsert_f(recover_type_id, "type_num", offset, &
                     H5T_NATIVE_INTEGER, error)
                  CALL h5dread_f(dset_id, recover_type_id, type_num, dims,&
                     error)
                  !ALLOCATE_BLOCK
                  IF (associated(abstr_type_ptr)) THEN
                     CALL read_DTYPE(cpfile_id, dset_id, abstr_type_ptr)
                  ELSE
                     WRITE (*,*) "CTYPE Allocation Failed"
                  ENDIF
                  CALL h5dclose_f(dset_id, error)
                  CALL h5gclose_f(group_id, error)
               ENDIF
               SELECT TYPE(abstr_type_ptr)
               CLASS is (CTYPE)
                  recover_CTYPE => abstr_type_ptr
               CLASS DEFAULT
                  WRITE (*,*) "Failed recovery of CTYPE"
               END SELECT
            END FUNCTION recover_CTYPE
END_TEMPLATE

# We use an if then block to decide which of our allocatable type we
# should make our recover object
my $allocate_stmt = "";
my $select_stmt = spaces() . "SELECT TYPE(type_ptr)\n";
while (my ($ind, $var) = each @allocatables) {
   $allocate_stmt .= spaces() . "IF (type_num == $ind) THEN\n";
   $allocate_stmt .= spaces() . "   ALLOCATE(${var}::abstr_type_ptr)\n";
   $allocate_stmt .= spaces() . "ENDIF\n";
   $select_stmt .= spaces() . "CLASS is ($var)\n";
   $select_stmt .= spaces() . "   type_num = $ind\n";
}
$select_stmt .= spaces() . "CLASS DEFAULT\n";
$select_stmt .= spaces() . "   WRITE (*,*) \"CTYPE: What are you?????\"\n";
$select_stmt .= spaces() . "END SELECT\n";

# Generate the functions for each derived type
sub generate_recover_function {
   my ($ctype) = @_;
   my @template = split '\n', $recover_function_template;
   for (@template){
      s/CTYPE/$ctype/g;
      s/ *!ALLOCATE_BLOCK/$allocate_stmt/;
   }
   return join ("\n", @template) . "\n";
}
my $recover_functions;
for (@derived_types) {
   $recover_functions .= generate_recover_function($_);
}
#print $recover_functions;

my $filename = $ARGV[1];
my $template = "type.f.in";
my $parenttype;
if ($ARGV[2]) { $template = $ARGV[2]; }
open TEMPLATE, "+<", $template or die $!;
open OUTFILE, ">", $filename or die $!;
print OUTFILE &spaces() . "! File autogenerated by my super script\n";

sub treetype {
   local $/ = "_";
   my $treetype = $root;
   chomp $treetype;
   return $treetype . "_tree";
}
my $treetype = treetype();
#print $treetype;

for (<TEMPLATE>) {

   # if its not an omitted line, print
   if (not ((m/dsize/ and (not $dsize))
            or (m/csize/ and (not $csize))
            or (m/dims/ and ($arrays == 0))
            or ((m/array_id/ or m/rank/) and ($arrays == 0) and ($numpointers==0))
            or (m/pointer_addr/ and ($numpointers == 0))
         #   or (m/!MTYPE/ and ($numtypes == 1))
         #   or (m/!STYPE/ and ($numtypes > 1))
         #or (m/subsize/ and ($numarr == 0))
          or ((/[sg]et_size/ or /PARENT/i) and (not $parenttype)))) {

      # Make the template substitutions
      $_ =~ s/ *!RECOVER_FUNCTIONS/$recover_functions/;
      $_ =~ s/DTYPE/$root/g;
      $_ =~ s/TREE_TYPE/$treetype/;
      $_ =~ s/ *!CREATE_STUB/$create/;
      $_ =~ s/ *!CALCULATE_STUB/$calc/;
      s/ *!WRITE_STUB/$write/;
      s/ *!SELECT_TYPE/$select_stmt/;
      s/ *!READ_STUB/$recover/;
      #s/ *!ALLOCATE_STUB/$allocate_stmt/ if (@allocatables);
      print OUTFILE;
   }
}
close TEMPLATE;
close OUTFILE;
