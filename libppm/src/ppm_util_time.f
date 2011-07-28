      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_util_time
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Retreives the current time. Uses either MPI_Wtime, 
      !                 etime or the f90 CPU_TIME intrinsic, based on the 
      !                 DEFINES that are set (see ppm_define.h):
      !                     __MPI       uses MPI_Wtime
      !                     __ETIME     uses etime
      !                     none        uses CPU_TIME
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : timing    (F) current CPU/wall clock time
      !
      !  Routines     : etime (C intrinsic)
      !
      !  Remarks      : etime is an C intrinsic - thus NOT standard fortran.
      !                 We therefore also write it in small letters.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_time.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2006/09/04 18:35:01  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.4  2004/06/10 16:20:06  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.3  2004/02/06 13:51:05  ivos
      !  Bugfix: fixed cpp directives for the serial version.
      !
      !  Revision 1.2  2004/01/23 17:22:13  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_time_s(timing)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_time_d(timing)
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      IMPLICIT NONE
#include "ppm_define.h"
      INCLUDE 'ppm_param.h'
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Type kind
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), INTENT(  OUT) :: timing
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
#ifdef __MPI
      INTEGER   :: info
#endif
#ifdef __ETIME
      REAL(MK) :: array(2)
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      REAL(MK), EXTERNAL :: etime
#endif
      
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Call the MPI function MPI_Wtime
      !-------------------------------------------------------------------------
      timing = MPI_Wtime()
#else
      !-------------------------------------------------------------------------
      !  Call the C routine: etime to get the cpu time
      !-------------------------------------------------------------------------
#ifdef __ETIME
      timing = etime(array)
#else
      CALL CPU_TIME(timing)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_time_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_time_d
#endif
