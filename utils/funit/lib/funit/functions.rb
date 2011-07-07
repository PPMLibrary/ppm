require 'erb'

module Funit

  TEST_RUNNER = ERB.new( %q{
    ! TestRunner.f - runs fUnit test suites
    !
    ! <%= File.basename $0 %> generated this file on <%= Time.now %>.

    program TestRunner

      <% test_suites.each do |test_suite| -%>
      use <%= File.basename(test_suite) %>_fun
      <% end -%>
      
      implicit none

      <% if use_mpi -%>
      INCLUDE 'mpif.h'
      <% end -%>

      integer, dimension(<%=test_suites.size%>) :: numTests, numAsserts, numAssertsTested, numFailures
      character(len=100)                        :: log_file_name
      integer                                   :: log = 20
      integer                                   :: comm
      integer                                   :: rank
      integer                                   :: size

      <% if use_mpi -%>
      integer :: mpiinfo
      call mpi_init(mpiinfo)
      <% end -%>

      rank = 0
      size = 1
      comm = 0

      <% if use_mpi -%>
      call mpi_comm_rank(MPI_COMM_WORLD, rank, mpiinfo)
      call mpi_comm_size(MPI_COMM_WORLD, size, mpiinfo)
      comm = MPI_COMM_WORLD
      <% end -%>

      write(log_file_name,'(A,I0,A)') 'test_runner.', rank, '.log'
      OPEN(log, FILE=log_file_name, ACTION='WRITE')
      write(log,*) "Starting new test run..."

      <% test_suites.each_with_index do |test_suite,i| -%>
      if (rank .eq. 0) then
         write(*,*)
         write(*,*) "<%= File.basename(test_suite) %> test suite:"
      end if
      write(log,*)
      write(log,*) "<%= File.basename(test_suite) %> test suite:"

      call test_<%= File.basename(test_suite) %> &
        ( numTests(<%= i+1 %>), numAsserts(<%= i+1 %>), numAssertsTested(<%= i+1 %>), &
          numFailures(<%= i+1 %>), log, rank, comm)

      write(*,1) rank, numAssertsTested(<%= i+1 %>), numAsserts(<%= i+1 %>), &
         numTests(<%= i+1 %>)-numFailures(<%= i+1 %>), numTests(<%= i+1 %>)
      write(log,1) rank, numAssertsTested(<%= i+1 %>), numAsserts(<%= i+1 %>), &
        numTests(<%= i+1 %>)-numFailures(<%= i+1 %>), numTests(<%= i+1 %>)

      <%= i+1 %> format('[',i0,'] Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')

      <% if use_mpi -%>
      call mpi_barrier(MPI_COMM_WORLD, mpiinfo)
      <% end -%>

      <% end -%>

      <% if use_mpi -%>
      call mpi_finalize(mpiinfo)
      <% end -%>

      if (rank .eq. 0) then
         write(*,*)
         write(*,'(a)') "======================[ SUMMARY ]======================"
      end if
      write(log,*)
      write(log,'(a)') "======================[ SUMMARY ]======================"
      <% max_length = test_suites.empty? ? 0 : test_suites.max.length -%>
      <% test_suites.each_with_index do |test_suite,i| -%>

      if (rank .eq. 0) then
         write(*,'(a<%=max_length+2%>)',advance="no") " <%= File.basename(test_suite) %>:"
      end if

      write(log,'(a<%=max_length+2%>)',advance="no") " <%= File.basename(test_suite) %>:"

      if ( numFailures(<%= i+1 %>) == 0 ) then
         if (rank .eq. 0) then
            write(*,*) " passed"
         end if
         write(log,*) " passed"
      else
         if (rank .eq. 0) then
            write(*,*) " failed   <<<<<"
         end if
         write(log,*) " failed   <<<<<"
      end if
      <% end -%>
 
      if (rank .eq. 0) then
         write(*,*)
      end if
      write(log,*)

      if ( sum(numFailures) /= 0 ) stop 1

    end program TestRunner
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
    <%= "\t(cd #{File.dirname(source)}; #{ENV['CPP']} #{ENV['DEFINE']} -DFUNIT_TEST -traditional-cpp -P #{File.basename(source)} __#{File.basename(source)} && #{ENV['FC']} #{ENV['FCFLAGS']} -I#{Dir.pwd+'/'+ENV['MODULES']} #{sourceflag} -c __#{File.basename(source)} -o #{File.basename(source, '.*')}.o && rm __#{File.basename(source)})" %>
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

  def warning( message, test_suite )
    $stderr.puts "\n *Warning: #{message} [#{test_suite}.fun:#$.]"
  end

  def compile_tests(test_suites,prog_source_dirs=['.'])
    puts "computing dependencies"

    sourceflag = ''
    if ENV['FSFLAG'] then
      sourceflag = prog_source_dirs.map{|pd| ENV['FSFLAG']+pd }.join(' ')
    end
    dependencies = Fortran::Dependencies.new(:search_paths=>prog_source_dirs)

    puts "locating associated source files and sorting for compilation"
    dependencies.source_file_dependencies('TestRunner.f')
    file_dependencies = dependencies.file_dependencies
    required_objects = file_dependencies.values.flatten.uniq.map{|s|s.sub(/\.f/i,'.o')}
    required_objects << 'TestRunner.o'

    File.open("makeTestRunner", "w") {|file| file.puts MAKEFILE.result(binding)}

    print "compiling..."
    compile = "make -f makeTestRunner > test_compile.log"
    raise "Compile failed." unless system compile
    system "rm test_compile.log"
    puts " done!"
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
