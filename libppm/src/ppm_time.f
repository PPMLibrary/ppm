      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_time
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Returns the current cpu time. Uses ppm_util_time,
      !                 which uses either MPI_Wtime, f90 CPU_TIME or etime,
      !                 based on the DEFINES that are set (see
      !                 ppm_define.h):
      !                     __MPI      uses MPI_Wtime
      !                     __ETIME    uses etime
      !                     none       uses f90 intrinsic CPU_TIME
      !
      !                 The difference in the returned timing value between
      !                 two subsequent calls gives the elapsed time in
      !                 seconds.
      !
      !  Input        : 
      !
      !  Input/output : 
      !
      !  Output       : timing     (F) current CPU clock time
      !                 info       (I) 0 on success.
      !
      !  Remarks      : This is just a wrapper to ppm_util_time. The user
      !                 is not allowed to call ppm_util_time directly as
      !                 the module_util is private.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_time.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/07/26 07:45:28  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/07/20 16:48:25  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_time_s(timing,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_time_d(timing,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      REAL(MK)               , INTENT(  OUT) :: timing
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_time',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_time',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------
      CALL ppm_util_time(timing)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_time',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_time_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_time_d
#endif
