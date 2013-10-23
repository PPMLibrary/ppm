
! TestRunner.f - runs fUnit test suites
!
! funit generated this file on 2013-10-23 14:15:48 +0200.

program TestRunner

    use ppm_module_map_part_fun
    
  implicit none

    INCLUDE 'mpif.h'
  
  integer, dimension(1) :: numTests, numAsserts, numAssertsTested, numFailures
  character(len=100)                        :: log_file_name
  integer                                   :: log = 20
  integer                                   :: comm
  integer                                   :: rank
  integer                                   :: size
  integer                                   :: i

    integer :: mpiinfo
  call mpi_init(mpiinfo)
  
  rank = 0
  size = 1
  comm = 0

    call mpi_comm_rank(MPI_COMM_WORLD, rank, mpiinfo)
  call mpi_comm_size(MPI_COMM_WORLD, size, mpiinfo)
  comm = MPI_COMM_WORLD
  
  write(log_file_name,'(A,I0,A)') 'test_runner.', rank, '.log'
  OPEN(log, FILE=log_file_name, ACTION='WRITE')
  write(log,*) "Starting new test run..."

    if (rank .eq. 0) then
     write(*,*)
     write(*,*) "ppm_module_map_part test suite:"
  end if
  write(log,*)
  write(log,*) "ppm_module_map_part test suite:"

  call test_ppm_module_map_part &
    ( numTests(1), numAsserts(1), numAssertsTested(1), &
      numFailures(1), log, rank, comm)

  write(*,1) rank, numAssertsTested(1), numAsserts(1), &
     numTests(1)-numFailures(1), numTests(1)
  write(log,1) rank, numAssertsTested(1), numAsserts(1), &
    numTests(1)-numFailures(1), numTests(1)

  1 format('[',i0,'] Passed ',i0,' of ',i0,' possible asserts comprising ',i0,' of ',i0,' tests.')

    call mpi_barrier(MPI_COMM_WORLD, mpiinfo)
  
  
    call mpi_finalize(mpiinfo)
  
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
    write(*,'(A)',advance='no') "ppm_module_map_part"
    do i=1,38
      write(*,'(A)',advance='no') " "
    end do
  end if

  do i=1,8
    write(log,'(A)',advance='no') " "
  end do
  write(log,'(A)',advance='no') "ppm_module_map_part"
  do i=1,38
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
