      !-------------------------------------------------------------------------
      !  Subroutine   :                     substop
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is called at the end of each subroutine.
      !                 Depending on the debug level, it prints the cpu time,
      !                 and calling routine.
      !
      !  Input        : caller   (C) Character string identifying the
      !                              calling subroutine
      !                 t0       (F) time at start of calling subroutine
      !                              (as returned by substart)
      !                 info     (I) The final info of the calling
      !                              subroutine. Is printed if debug.GT.1.
      !
      !  Input/output : 
      !
      !  Output       : 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: substop.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2004/07/26 11:49:57  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.6  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.5  2004/06/10 16:20:06  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.4  2004/01/26 12:34:38  ivos
      !  Bugfix: info is now INTENT(IN) and is NOT changed, but included in
      !  the printed message if debug.GT.1.
      !
      !  Revision 1.3  2004/01/23 17:22:14  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.2  2004/01/13 13:02:16  ivos
      !  Added comment header.
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
      SUBROUTINE substop_s(caller,t0,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE substop_d(caller,t0,info)
#endif
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
      REAL(MK)        , INTENT(IN   ) :: t0
      INTEGER         , INTENT(IN   ) :: info
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
