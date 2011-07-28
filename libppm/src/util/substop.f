      !-------------------------------------------------------------------------
      !  Subroutine   :                     substop
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE substop_s(caller,t0,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE substop_d(caller,t0,info)
#endif
      !!! This routine is called at the end of each subroutine.Depending on
      !!! the debug level, it prints the cpu time, and calling routine.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN   ) :: caller
      !!! Character string with the name of the calling subroutine
      REAL(MK)        , INTENT(IN   ) :: t0
      !!! System/cpu time at start of subroutine. (as returned by substart)
      INTEGER         , INTENT(IN   ) :: info
      !!! The final info of the calling subroutine. Is printed if debug.GT.1.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t1
      CHARACTER(LEN=ppm_char)        :: cbuf
      INTEGER                        :: info2 = 0

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  In parallel collect the MIN info of the processors
      !  THIS IS COMMENTED OUT AS IT IS VERY EXPENSIVE.
      !-------------------------------------------------------------------------
!     CALL MPI_AllReduce(info,i,1,MPI_INTEGER,MPI_MIN,comm,info2)
!     info = i
#endif
      !-------------------------------------------------------------------------
      !  Using debugging print the 'leaving <routine>'
      !-------------------------------------------------------------------------
      CALL ppm_util_time(t1)
      IF     (ppm_debug.GT.1) THEN
         WRITE(cbuf,'(A,I2,A,E12.4)') 'leaving with info=',info,'. took: ',t1-t0
         CALL ppm_write(ppm_rank,caller,cbuf,info2)
      ELSEIF (ppm_debug.GT.0) THEN
         WRITE(cbuf,'(A,E12.4)') 'took: ',t1-t0
         CALL ppm_write(ppm_rank,caller,cbuf,info2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE substop_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE substop_d
#endif
