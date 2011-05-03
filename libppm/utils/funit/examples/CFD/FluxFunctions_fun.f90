! FluxFunctions_fun.f90 - a unit test suite for FluxFunctions.f90
!
! funit generated this file from FluxFunctions.fun

module FluxFunctions_fun

 use FluxFunctions

 implicit none

 logical :: noAssertFailed

 public :: test_FluxFunctions

 private

 integer :: numTests          = 0
 integer :: numAsserts        = 0
 integer :: numAssertsTested  = 0
 integer :: numFailures       = 0



real :: leftState, rightState, interfaceFlux

 contains


 subroutine FluxZero

  real :: state
  state = 0
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((Flux(state) &
     +0.00001 ) &
     .ge. &
     ( 0) &
             .and. &
     (Flux(state) &
     -0.00001 ) &
     .le. &
     ( 0) )) then
      print *, " *Assert_Equal_Within failed* in test FluxZero &
              &[FluxFunctions.fun:13]"
      print *, "  ", " 0 (", 0,") is not", &
 Flux(state),"within",0.00001 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine FluxZero

 
 subroutine FluxOne

  real :: state = 1
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((Flux(state) &
     +0.00001 ) &
     .ge. &
     ( 0.5) &
             .and. &
     (Flux(state) &
     -0.00001 ) &
     .le. &
     ( 0.5) )) then
      print *, " *Assert_Equal_Within failed* in test FluxOne &
              &[FluxFunctions.fun:18]"
      print *, "  ", " 0.5 (", 0.5,") is not", &
 Flux(state),"within",0.00001 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine FluxOne


 subroutine RoeAvgZero

  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( ( 0 &
        +2*spacing(real( 0)) ) &
        .ge. &
        (RoeAvg(0.0,0.0) ) &
            .and. &
     ( 0 &
      -2*spacing(real( 0)) ) &
      .le. &
       (RoeAvg(0.0,0.0) ) )) then
      print *, " *Assert_Real_Equal failed* in test RoeAvgZero &
              &[FluxFunctions.fun:22]"
      print *, "  ", "RoeAvg(0.0,0.0)  (", &
 RoeAvg(0.0,0.0) , &
  ") is not", &
  0,&
 "within", &
  2*spacing(real( 0))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_False assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if ( RoeAvg(0.0,0.0)==1 ) then
      print *, " *Assert_False failed* in test RoeAvgZero &
              &[FluxFunctions.fun:23]"
      print *, "  ", " RoeAvg(0.0,0.0)==1  is not false"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine RoeAvgZero


 subroutine RoeAvgKnown

  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( ( 0.5 &
        +2*spacing(real( 0.5)) ) &
        .ge. &
        (RoeAvg(leftState,rightState) ) &
            .and. &
     ( 0.5 &
      -2*spacing(real( 0.5)) ) &
      .le. &
       (RoeAvg(leftState,rightState) ) )) then
      print *, " *Assert_Real_Equal failed* in test RoeAvgKnown &
              &[FluxFunctions.fun:27]"
      print *, "  ", "RoeAvg(leftState,rightState)  (", &
 RoeAvg(leftState,rightState) , &
  ") is not", &
  0.5,&
 "within", &
  2*spacing(real( 0.5))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_True assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( RoeAvg(leftState,rightState) > 0 )) then
      print *, " *Assert_True failed* in test RoeAvgKnown &
              &[FluxFunctions.fun:28]"
      print *, "  ", " RoeAvg(leftState,rightState) > 0  is not true"
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine RoeAvgKnown


 subroutine CentralFluxKnown

  call CentralFlux( leftState, rightState, interfaceFlux )
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((interfaceFlux &
     +0.001 ) &
     .ge. &
     ( 0.25) &
             .and. &
     (interfaceFlux &
     -0.001 ) &
     .le. &
     ( 0.25) )) then
      print *, " *Assert_Equal_Within failed* in test CentralFluxKnown &
              &[FluxFunctions.fun:33]"
      print *, "  ", " 0.25 (", 0.25,") is not", &
 interfaceFlux,"within",0.001 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal_Within assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.((interfaceFlux &
     +0.00000001 ) &
     .ge. &
     ( 0.25) &
             .and. &
     (interfaceFlux &
     -0.00000001 ) &
     .le. &
     ( 0.25) )) then
      print *, " *Assert_Equal_Within failed* in test CentralFluxKnown &
              &[FluxFunctions.fun:34]"
      print *, "  ", " 0.25 (", 0.25,") is not", &
 interfaceFlux,"within",0.00000001 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( 0.25== interfaceFlux )) then
      print *, " *Assert_Equal failed* in test CentralFluxKnown &
              &[FluxFunctions.fun:35]"
      print *, "  ", " 0.25 (", 0.25,") is not",  interfaceFlux 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine CentralFluxKnown


 subroutine RoeFluxExpansionShock

  leftState = -1
  call RoeFlux( leftState, rightState, interfaceFlux )
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( 0.5== interfaceFlux )) then
      print *, " *Assert_Equal failed* in test RoeFluxExpansionShock &
              &[FluxFunctions.fun:41]"
      print *, "  ", " 0.5 (", 0.5,") is not",  interfaceFlux 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine RoeFluxExpansionShock


 subroutine RoeFluxZero

  rightState = 0
  call RoeFlux( leftState, rightState, interfaceFlux )
  ! Assert_Real_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( ( 0 &
        +2*spacing(real( 0)) ) &
        .ge. &
        (interfaceFlux ) &
            .and. &
     ( 0 &
      -2*spacing(real( 0)) ) &
      .le. &
       (interfaceFlux ) )) then
      print *, " *Assert_Real_Equal failed* in test RoeFluxZero &
              &[FluxFunctions.fun:47]"
      print *, "  ", "interfaceFlux  (", &
 interfaceFlux , &
  ") is not", &
  0,&
 "within", &
  2*spacing(real( 0))
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif
  ! Assert_Equal assertion
  numAsserts = numAsserts + 1
  if (noAssertFailed) then
    if (.not.( 0== interfaceFlux )) then
      print *, " *Assert_Equal failed* in test RoeFluxZero &
              &[FluxFunctions.fun:48]"
      print *, "  ", " 0 (", 0,") is not",  interfaceFlux 
      print *, ""
      noAssertFailed = .false.
      numFailures    = numFailures + 1
    else
      numAssertsTested = numAssertsTested + 1
    endif
  endif

  numTests = numTests + 1

 end subroutine RoeFluxZero


 subroutine funit_setup
  leftState  = 0
  rightState = 1
  noAssertFailed = .true.
 end subroutine funit_setup


 subroutine funit_teardown
 end subroutine funit_teardown


 subroutine test_FluxFunctions( nTests, nAsserts, nAssertsTested, nFailures )

  integer :: nTests
  integer :: nAsserts
  integer :: nAssertsTested
  integer :: nFailures

  continue

  call funit_setup
  call FluxZero
  call funit_teardown

  call funit_setup
  call FluxOne
  call funit_teardown

  call funit_setup
  call RoeAvgZero
  call funit_teardown

  call funit_setup
  call RoeAvgKnown
  call funit_teardown

  call funit_setup
  call CentralFluxKnown
  call funit_teardown

  call funit_setup
  call RoeFluxExpansionShock
  call funit_teardown

  call funit_setup
  call RoeFluxZero
  call funit_teardown

  nTests          = numTests
  nAsserts        = numAsserts
  nAssertsTested  = numAssertsTested
  nFailures       = numFailures

 end subroutine test_FluxFunctions

end module FluxFunctions_fun
