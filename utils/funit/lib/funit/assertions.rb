
module Funit

  ##
  # Fortran assertion macro definitions

  module Assertions

    ASSERTION_PATTERN =
     /^\s*?(assert_(array_equal|larray_equal|real_equal|false|true|equal_within|equal))\(.*\)/i

    def assert_true(line)
      line.match(/\((.+)\)/)
      @line = line
      @type = 'Assert_True'
      @condition = ".not.(#$1)"
      @message = "\"#$1 is not true\""
      syntax_error("invalid body for #@type",@suite_name) unless $1=~/\S+/
      write_assert
    end

    def assert_false(line)
      line.match(/\((.+)\)/)
      @line = line
      @type = 'Assert_False'
      @condition = "#$1"
      @message = "\"#$1 is not false\""
      syntax_error("invalid body for #@type",@suite_name) unless $1=~/\S+/
      write_assert
    end

    def assert_real_equal(line)
      line.match(/\((.*)\)/)
      @line = line
      expected, actual = *(get_args($1))
      @type = 'Assert_Real_Equal'
      @condition = ".not.( (#{expected} &\n        +2*spacing(real(#{expected})) ) &\n        .ge. &\n        (#{actual}) &\n            .and. &\n     (#{expected} &\n      -2*spacing(real(#{expected})) ) &\n      .le. &\n       (#{actual}) )"
      @message = "\"#{actual} (\", &\n #{actual}, &\n  \") is not\", &\n #{expected},\&\n \"within\", &\n  2*spacing(real(#{expected}))"
      syntax_error("invalid body for #@type",@suite_name) unless $&
      write_assert
    end

    def assert_equal_within(line)
      line.match(/\((.*)\)/)
      @line = line
      expected, actual, tolerance = *(get_args($1))
      @type = 'Assert_Equal_Within'
      @condition = ".not.((#{actual} &\n     +#{tolerance}) &\n     .ge. &\n     (#{expected}) &\n             .and. &\n     (#{actual} &\n     -#{tolerance}) &\n     .le. &\n     (#{expected}) )"
      @message = "\"#{expected} (\",#{expected},\") is not\", &\n #{actual},\"within\",#{tolerance}"
      syntax_error("invalid body for #@type",@suite_name) unless $&
      write_assert
    end

    def assert_equal(line)
      line.match(/\((\w+\(.*\)|[^,]+),(.+)\)/)
      @line = line
      @type = 'Assert_Equal'
      @condition = ".not.(#$1==#$2)"
      @message = "\"#$1 (\",#$1,\") is not\", #$2"
      syntax_error("invalid body for #@type",@suite_name) unless $&
      write_assert
    end

    def assert_array_equal(line)
      line.match(/\(\s*(\w+)\s*,\s*(\w+|\(\/.*\/\))\s*\)/)
      @line = line
      @type = 'Assert_Array_Equal'
      @condition = ".not. all(#$1==#$2)"
      @message = "\"array #$1 is not #$2\""
      syntax_error("invalid body for #@type",@suite_name) unless $&
      write_assert
    end

    def assert_larray_equal(line)
      line.match(/\(\s*(\w+)\s*,\s*(\w+|\(\/.*\/\))\s*\)/)
      @line = line
      @type = 'Assert_Array_Equal'
      @condition = ".not. all(#$1.eqv.#$2)"
      @message = "\"array #$1 is not #$2\""
      syntax_error("invalid body for #@type",@suite_name) unless $&
      write_assert
    end


    ##
    # An argument scanner thanks to James Edward Gray II
    # by way of ruby-talk mailing list.

    def get_args(string)
      scanner = ::StringScanner.new(string)
      result  = scanner.eos? ? [] : ['']
      paren_depth = 0
      until scanner.eos?
        if scanner.scan(/[^(),]+/)
          # do nothing--we found the part of the argument we need to add
        elsif scanner.scan(/\(/)
          paren_depth += 1
        elsif scanner.scan(/\)/)
          paren_depth -= 1
        elsif scanner.scan(/,\s*/) and paren_depth.zero?
          result << ''
          next
        end
        result.last << scanner.matched
      end
      result
    end

    ##
    # Translate the assertion to Fortran.

    def write_assert
      <<-OUTPUT
      ! #@type assertion
      numAsserts = numAsserts + 1
      IF (noAssertFailed) THEN
         WRITE(log,'(A)', ADVANCE='NO') "trying assert #{@line.gsub(/"/, "'")}"
         IF (#@condition) THEN
            WRITE(log,*) " failed!"
            WRITE(*,'(A,I0,A)', ADVANCE='NO') "[", funit_rank, "] "
            WRITE(*,'(A)', ADVANCE='NO') "*#@type failed* in test #@test_name &
            [#{@suite_name}.fun:#{@line_number.to_s}] "
            WRITE(*,*) #@message
            WRITE(log,*) " *#@type failed* in test #@test_name &
            &[#{@suite_name}.fun:#{@line_number.to_s}]"
            WRITE(log,*) "  ", #@message
            WRITE(log,*) ""
            noAssertFailed = .FALSE.
            numFailures    = numFailures + 1
#ifdef __MPI
            CALL MPI_Abort(funit_comm, funit_info)
#endif
         ELSE
            WRITE(log,*) " success!"
            numAssertsTested = numAssertsTested + 1
         ENDIF
      ENDIF
      OUTPUT
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
