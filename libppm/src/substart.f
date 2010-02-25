      !-------------------------------------------------------------------------
      !  Subroutine   :                     substart
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is called whenever a subroutine is
      !                 started. It initializes the info of that subroutine
      !                 to 0 and prints a debug message if ppm_debug is
      !                 > 1.
      !
      !  Input        : caller    (C) character string with the name of the
      !                               calling subroutine
      !
      !  Input/output : 
      !
      !  Output       : t0        (F) system/cpu time at start of
      !                               subroutine. Only returned if
      !                               ppm_debug .GT. 0
      !                 info      (I) initialized info for the calling
      !                               subroutine.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: substart.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/07/26 11:49:56  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/07/20 16:46:28  ivos
      !  Changed to use ppm_param_success.
      !
      !  Revision 1.3  2004/02/10 14:33:05  walther
      !  Now calling the timing in all debugging cases.
      !
      !  Revision 1.2  2004/01/23 17:22:14  ivos
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
      SUBROUTINE substart_s(caller,t0,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE substart_d(caller,t0,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN   ) :: caller
      REAL(MK)        , INTENT(  OUT) :: t0
      INTEGER         , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                         :: info2
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      info = ppm_param_success
      IF     (ppm_debug.GT.1) THEN
         CALL ppm_util_time(t0)
         CALL ppm_write(ppm_rank,caller,'entering',info2)
      ELSEIF (ppm_debug.GT.0) THEN
         CALL ppm_util_time(t0)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE substart_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE substart_d
#endif

