
! TestRunner.f - runs fUnit test suites
!
! funit generated this file on 2013-09-04 15:23:18 +0200.

program TestRunner

    use yaser_fun
    
  implicit none

  
  integer, dimension(1) :: numTests, numAsserts, numAssertsTested, numFailures
  character(len=100)                        :: log_file_name
  integer                                   :: log = 20
  integer                                   :: comm
  integer                                   :: rank
  integer                                   :: size
  integer                                   :: i

  
  rank = 0
  size = 1
  comm = 0

  
  write(log_file_name,'(A,I0,A)') 'test_runner.', rank, '.log'
  OPEN(log, FILE=log_file_name, ACTION='WRITE')
  write(log,*) "Starting new test run..."

    if (rank .eq. 0) then
     write(*,*)
     write(*,*) "yaser test suite:"
  end if
  write(log,*)
  write(log,*) "yaser test suite:"

  call test_yaser &
    ( numTests(1), numAsserts(1), numAssertsTested(1), &
      numFailures(1), log, rank, comm)

  write(*,1) rank, numAssertsTested(1), numAsserts(1), &
     numTests(1)-numFailures(1), numTests(1)
  write(log,1) rank, numAssertsTested(1), numAsserts(1), &
    numTests(1)-numFailures(1), numTests(1)

  1 format('[',i0,'] Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')

  
  
  
  if (rank .eq. 0) then
     write(*,*)
     write(*,'(a/)') "==================================[ SUMMARY ]==================================="
  end if
  write(log,*)
  write(log,'(a/)') "==================================[ SUMMARY ]==================================="
    
  if (rank .eq. 0) then
    do i=1,8
      write(*,'(A)',advance='no') " "
    end do
    write(*,'(A)',advance='no') "yaser"
    do i=1,52
      write(*,'(A)',advance='no') " "
    end do
  end if

  do i=1,8
    write(log,'(A)',advance='no') " "
  end do
  write(log,'(A)',advance='no') "yaser"
  do i=1,52
    write(log,'(A)',advance='no') " "
  end do

  if ( numFailures(1) == 0 ) then
     if (rank .eq. 0) then
        write(*,*) " passed"
     end if
     write(log,*) " passed"
  else
     if (rank .eq. 0) then
        write(*,*) " failed <<<<<<"
     end if
     write(log,*) " failed <<<<<<"
  end if
   
  if (rank .eq. 0) then
     write(*,*)
  end if
  write(log,*)

  if ( sum(numFailures) /= 0 ) stop 1

end program TestRunner
