require 'erb'

module Funit

  TEST_RUNNER = ERB.new( %q{
    ! TestRunner.f - runs fUnit test suites
    !
    ! <%= File.basename $0 %> generated this file on <%= Time.now %>.

    PROGRAM TestRunner

      <% test_suites.each do |test_suite| -%>
      USE <%= File.basename(test_suite) %>_fun
      <% end -%>

      IMPLICIT NONE

      <% if use_mpi -%>
      INCLUDE 'mpif.h'
      <% end -%>

      INTEGER, DIMENSION(<%=test_suites.size%>) :: numTests, numAsserts, numAssertsTested, numFailures
      CHARACTER(LEN=100)                        :: log_file_name
      INTEGER                                   :: log = 20
      INTEGER                                   :: comm
      INTEGER                                   :: rank
      INTEGER                                   :: size
      INTEGER                                   :: i

      <% if use_mpi -%>
      INTEGER :: mpiinfo
      CALL MPI_Init(mpiinfo)
      <% end -%>

      rank = 0
      size = 1
      comm = 0

      <% if use_mpi -%>
      CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, mpiinfo)
      CALL MPI_Comm_size(MPI_COMM_WORLD, size, mpiinfo)
      COMM = MPI_COMM_WORLD
      <% end -%>

      WRITE(log_file_name,'(A,I0,A)') 'test_runner.', rank, '.log'
      OPEN(log, FILE=log_file_name, ACTION='WRITE')
      WRITE(log,*) "Starting new test run..."

      <% test_suites.each_with_index do |test_suite,i| -%>
      IF (rank .EQ. 0) THEN
         WRITE(*,*)
         WRITE(*,*) "<%= File.basename(test_suite) %> test suite:"
      END IF
      WRITE(log,*)
      WRITE(log,*) "<%= File.basename(test_suite) %> test suite:"

      CALL test_<%= File.basename(test_suite) %> &
        ( numTests(<%= i+1 %>), numAsserts(<%= i+1 %>), numAssertsTested(<%= i+1 %>), &
          numFailures(<%= i+1 %>), log, rank, comm)

      WRITE(*,1) rank, numAssertsTested(<%= i+1 %>), numAsserts(<%= i+1 %>), &
         numTests(<%= i+1 %>)-numFailures(<%= i+1 %>), numTests(<%= i+1 %>)
      WRITE(log,1) rank, numAssertsTested(<%= i+1 %>), numAsserts(<%= i+1 %>), &
        numTests(<%= i+1 %>)-numFailures(<%= i+1 %>), numTests(<%= i+1 %>)

      <%= i+1 %> format('[',i0,'] Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')

      <% if use_mpi -%>
      CALL MPI_Barrier(MPI_COMM_WORLD, mpiinfo)
      <% end -%>

      <% end -%>

      <% if use_mpi -%>
      CALL MPI_Finalize(mpiinfo)
      <% end -%>

      IF (rank .EQ. 0) THEN
         WRITE(*,*)
         WRITE(*,'(a/)') "==================================[ SUMMARY ]==================================="
      END IF
      WRITE(log,*)
      WRITE(log,'(a/)') "==================================[ SUMMARY ]==================================="
      <% max_length = test_suites.empty? ? 0 : test_suites.max.length -%>
      <% test_suites.each_with_index do |test_suite,i| -%>

      IF (rank .EQ. 0) THEN
        DO i=1,<%= OUTPUT_INDENT %>
          WRITE(*,'(A)',advance='no') " "
        END DO
        WRITE(*,'(A)',advance='no') "<%= File.basename(test_suite) %>"
        DO i=1,<%= OUTPUT_WIDTH - 2 * OUTPUT_INDENT - File.basename(test_suite).length - 7 %>
          WRITE(*,'(A)',advance='no') " "
        END DO
      END IF

      DO i=1,<%= OUTPUT_INDENT %>
        WRITE(log,'(A)',advance='no') " "
      END DO
      WRITE(log,'(A)',advance='no') "<%= File.basename(test_suite) %>"
      DO i=1,<%= OUTPUT_WIDTH - 2 * OUTPUT_INDENT - File.basename(test_suite).length - 7 %>
        WRITE(log,'(A)',advance='no') " "
      END DO

      IF ( numFailures(<%= i+1 %>) == 0 ) THEN
         IF (rank .EQ. 0) THEN
            WRITE(*,*) " passed"
         END IF
         WRITE(log,*) " passed"
      ELSE
         IF (rank .EQ. 0) THEN
            WRITE(*,*) " failed <<<<<<"
         END IF
         WRITE(log,*) " failed <<<<<<"
      END IF
      <% end -%>

      if (rank .eq. 0) THEN
         WRITE(*,*)
      END IF
      WRITE(log,*)

      IF ( sum(numFailures) /= 0 ) STOP 1

    END PROGRAM TestRunner
    }.gsub(/^    /,''), nil, '-' ) # turn off newlines for <% -%>

  MAKEFILE = ERB.new( %q{
    # makefile to compile TestRunner.f
    #
    # <%= File.basename $0 %> generated this file on <%= Time.now %>.

    OBJ=<%= required_objects.join(' ') %>

    all:testrunner

    testrunner: $(OBJ)
    <%= "\t#{ENV['FC']} $(OBJ) #{ENV['FCFLAGS']} #{ENV['LDFLAGS']} #{sourceflag}" %> -o TestRunner

    <% file_dependencies.each do |source,dep| -%>
    <%= "#{source.sub(/\.f/i,'.o')}: #{source} #{dep.map{ |d| d.sub(/\.f/i,'.o') }.join(' ')}" %>
    <%= "\t(cd #{File.dirname(source)}; #{ENV['CPP']} #{ENV['DEFINE']} -DFUNIT_TEST -traditional-cpp -P #{File.basename(source)} | #{ENV['CG']} > __#{File.basename(source)} && #{ENV['FC']} #{ENV['FCFLAGS']} -I#{Dir.pwd+'/'+ENV['MODULES']} #{sourceflag} -c __#{File.basename(source)} -o #{File.basename(source, '.*')}.o)" %>
    <% end -%>
  }.gsub(/^    /,''), nil, '-' ) # turn off newlines for <% -%>

  def requested_modules(module_names)
    if module_names.empty?
      module_names = Dir["**/test/*.fun"].each{ |mod| mod.chomp! ".fun" }
    end
    module_names
  end

  def funit_exists?(module_name)
    File.exists? "#{module_name}.fun"
  end

  def parse_command_line

    module_names = requested_modules(ARGV)

    if module_names.empty?
      raise "   *Error: no test suites found in this directory"
    end

    module_names.each do |mod|
      unless funit_exists?(mod)
        error_message = <<-FUNIT_DOES_NOT_EXIST
 Error: could not find test suite #{mod}.fun
 Test suites available in this directory:
 #{requested_modules([]).join(' ')}

 Usage: #{File.basename $0} [test names (w/o .fun suffix)]
        FUNIT_DOES_NOT_EXIST
        raise error_message
      end
    end

  end

  def write_test_runner( test_suites, use_mpi )
    File.open("TestRunner.f", "w") do |file|
      file.puts TEST_RUNNER.result(binding)
    end
  end

  def syntax_error( message, test_suite )
    raise "\n   *Error: #{message} [#{test_suite}.fun:#$.]\n\n"
  end

  #def warning( message, test_suite )
  #  $stderr.puts "\n *Warning: #{message} [#{test_suite}.fun:#$.]"
  #end

  def compile_tests(test_suites,prog_source_dirs=['.'])
    print_sub("compile")

    print_started("computing dependencies")
    print_started("")
    sourceflag = ''
    if ENV['FSFLAG'] then
      sourceflag = prog_source_dirs.map{|pd| ENV['FSFLAG']+pd }.join(' ')
    end
    dependencies = Fortran::Dependencies.new(:search_paths=>prog_source_dirs)

    print_done

    print_started("locating sources")

    dependencies.source_file_dependencies('TestRunner.f')
    file_dependencies = dependencies.file_dependencies
    required_objects = file_dependencies.values.flatten.uniq.map{|s|s.sub(/\.f/i,'.o')}
    required_objects << 'TestRunner.o'

    print_done

    print_started("writing makefile")
    File.open("makeTestRunner", "w") {|file| file.puts MAKEFILE.result(binding)}
    print_done

    print_started("compiling")
    compile = "make -f makeTestRunner > test_compile.log"
    raise "Compile failed." unless system compile
    system "rm test_compile.log"
    print_done
  end

# pretty printing

  OUTPUT_WIDTH = 80
  OUTPUT_INDENT = 8

  def print_heading(text, line)
    left = OUTPUT_WIDTH - (text.length + 4)
    if (left % 2 == 0) then
      left = left / 2
      right = left
    else
      left = (left - 1) / 2
      right = left + 1
    end
    left.times { print line }
    print "[ " + text + " ]"
    right.times { print line }
    print "\n"
  end

  def print_title(text)
    print "\n\n"
    print_heading(text, '=')
    print "\n"
  end

  def print_sub(text)
    print "\n"
    print_heading(text, '-')
    print "\n"
  end

  @current_position = 0

  def print_started(text)
    print "\n"
    OUTPUT_INDENT.times { print " " }
    print text
    @current_position = OUTPUT_INDENT + text.length
  end

  def print_done(text="done!")
    space = (OUTPUT_WIDTH - @current_position) - text.length - OUTPUT_INDENT
    if (space < 1) then
      space = 1
    end
    space.times { print " " }

    puts text

    @current_position = 0
  end

end

#--
# Copyright 2006-2007 United States Government as represented by
# NASA Langley Research Center. No copyright is claimed in
# the United States under Title 17, U.S. Code. All Other Rights
# Reserved.
#
# This file is governed by the NASA Open Source Agreement.
# See License.txt for details.
#++
