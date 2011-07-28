require 'funit'

module Funit
  
  include Assertions # FIXME

  ##
  # Create testsuite wrapper code

  class TestSuite < File

    KEYWORDS = Regexp.union(/^\s*(end\s+)?(setup|teardown|test|init|finalize)/i,
                            Assertions::ASSERTION_PATTERN)
    COMMENT_LINE = /^\s*!/
    
    include Funit #FIXME

    def initialize( suite_name, suite_content, trailing_code, wrap_with_module )
      @line_number = 'blank'
      @suite_name = suite_name
      @suite_content = suite_content
      @trailing_code = trailing_code
      return nil unless funit_exists?(suite_name)
      File.delete(suite_name+"_fun.f") if File.exists?(suite_name+"_fun.f")
      super(suite_name+"_fun.f","w")
      @init, @finalize, @tests, @setup, @teardown = [], [], [], [], []
      header
      @wrap_with_module = wrap_with_module
      module_wrapper if @wrap_with_module
      top_wrapper
      expand
      close
    end

    def header
      puts <<-HEADER
! #{@suite_name}_fun.f - a unit test suite for #{@suite_name}.f
!
! #{File.basename $0} generated this file from #{@suite_name}.fun

      HEADER
    end

    def module_wrapper
      puts <<-MODULE_WRAPPER
module #{@suite_name}_mod
contains
  include '#@suite_name.f'
end module #{@suite_name}_mod

      MODULE_WRAPPER
    end

    def top_wrapper
      puts <<-TOP
module #{@suite_name}_fun

 use #{ @wrap_with_module ? @suite_name+'_mod' : @suite_name }

 implicit none

 logical :: noAssertFailed
 integer :: log

 public :: test_#@suite_name

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0

      TOP
    end

    def expand
      $stderr.print "expanding test suite: #{@suite_name}..."
      funit_contents = @suite_content.split("\n")
      @funit_total_lines = funit_contents.length

      while (line = funit_contents.shift) && line !~ KEYWORDS
        puts line
      end

      funit_contents.unshift line

      puts " contains\n\n"

      while (line = funit_contents.shift)
        case line
        when COMMENT_LINE
          puts line
        when /^\s*init/i
          add_to_init funit_contents
        when /^\s*finalize/i
          add_to_finalize funit_contents
        when /^\s*setup/i
          add_to_setup funit_contents
        when /^\s*teardown/i
          add_to_teardown funit_contents
        when /^\s*Xtest\s+(\w+)/i
          ignore_test($1,funit_contents)
        when /^\s*test\s+(\w+)/i
          a_test($1,funit_contents)
        when /^\s*test/i
          syntax_error "no name given for test", @suite_name
        when /^\s*end\s+(setup|teardown|test)/i
          syntax_error "no matching #$1 for an #$&", @suite_name
        when Assertions::ASSERTION_PATTERN
          syntax_error "#$1 assertion not in a test block", @suite_name
        else
          puts line
        end
      end
      $stderr.puts "done."
    end

    def add_to_init funit_contents
      while (line = funit_contents.shift) && line !~ /end\s+init/i
        @init.push line
      end
    end
    
    def add_to_finalize funit_contents
      while (line = funit_contents.shift) && line !~ /end\s+finalize/i
        @finalize.push line
      end
    end
    
    def add_to_setup funit_contents
      while (line = funit_contents.shift) && line !~ /end\s+setup/i
        @setup.push line
      end
    end

    def add_to_teardown funit_contents
      while (line = funit_contents.shift) && line !~ /end\s+teardown/i
        @teardown.push line
      end
    end

    def ignore_test test_name, funit_contents
      warning("Ignoring test: #{test_name}", @suite_name)
      line = funit_contents.shift while line !~ /end\s+test/i
    end

    def a_test test_name, funit_contents
      @test_name = test_name
      @tests.push test_name
      syntax_error("test name #@test_name not unique",@suite_name) if (@tests.uniq!)

      puts " subroutine #{test_name}\n\n"

      num_of_asserts = 0
  
      while (line = funit_contents.shift) && line !~ /(end\s+test|contains)/i
        case line
        when COMMENT_LINE
          puts line
        when Assertions::ASSERTION_PATTERN
          @line_number = @funit_total_lines - funit_contents.length
          num_of_asserts += 1
          puts send( $1.downcase, line )
        else
          puts line
        end
      end

      warning("no asserts in test", @suite_name) if num_of_asserts == 0

      puts "\n  numTests = numTests + 1\n\n"

      if line =~ /contains/i
        puts line
        while (line = funit_contents.shift) && line !~ /end\s+test/
          puts line
        end
      end

      puts " end subroutine #{test_name}\n\n"
    end

    def close
      puts "\n subroutine funit_init"
      puts @init
      puts " end subroutine funit_init\n\n"
      
      puts "\n subroutine funit_setup"
      puts @setup
      puts "  noAssertFailed = .true."
      puts " end subroutine funit_setup\n\n"

      puts "\n subroutine funit_teardown"
      puts @teardown
      puts " end subroutine funit_teardown\n\n"
      
      puts "\n subroutine funit_finalize"
      puts @finalize
      puts " end subroutine funit_finalize\n\n"

      puts <<-NEXTONE

 subroutine test_#{@suite_name}( nTests, nAsserts, nAssertsTested, nFailures, lfh )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures
  integer :: lfh

  continue

  log = lfh  

  call funit_init

      NEXTONE

      @tests.each do |test_name|
        puts "\n  write(log,*) 'setting up...'"
        puts "  call funit_setup"
        puts "  write(log,*) 'Entering #{test_name}...'\n"
        puts "  call #{test_name}"
        puts "  write(log,*) 'Leaving #{test_name}...'\n"
        puts "  call funit_teardown"
        puts "  write(log,*) 'cleaned up...'"
      end

      puts <<-LASTONE
  
  call funit_finalize

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_#{@suite_name}

end module #{@suite_name}_fun
      LASTONE
      puts @trailing_code
      super
    end

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
