      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_time
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_time_s(timing,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_time_d(timing,info)
#endif
      !!! Returns the current cpu time. Uses `ppm_util_time`,
      !!! which uses either `MPI_Wtime`, f90 `CPU_TIME` or `etime`,
      !!! based on how PPM was configured (see
      !!! ./configure --help):
      !!!
      !!! The difference in the returned timing value between
      !!! two subsequent calls gives the elapsed time in
      !!! seconds.
      !!!
      !!! [NOTE]
      !!! This is just a wrapper to `ppm_util_time`. The user is not allowed to
      !!! call `ppm_util_time` directly.
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
      !!! Current CPU clock time
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
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
