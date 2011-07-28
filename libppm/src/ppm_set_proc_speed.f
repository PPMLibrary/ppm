      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_set_proc_speed
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine can be used by the user to set the
      !                 relative speeds of the processors (will be used in
      !                 subs2proc for load balancing).
      !
      !  Input        : proc_speed(:) (F) Relative speeds of all processors
      !                                   from 0 to ppm_nproc-1. The numbers
      !                                   must sum up to 1.
      !
      !  Input/output : 
      !
      !  Output       : info       (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_set_proc_speed.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.7  2006/07/18 09:01:37  ivos
      !  Unrolled SUM into DO loops and removed double computations.
      !
      !  Revision 1.6  2004/07/26 07:45:25  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.5  2004/07/20 16:46:01  ivos
      !  Removed unnecessary cpp ifs.
      !
      !  Revision 1.4  2004/07/16 14:47:22  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/05/06 07:45:16  ivos
      !  Updated comment header.
      !
      !  Revision 1.2  2004/02/26 14:46:26  ivos
      !  bugfix: lower bound of assumed-size argument was missing. 
      !  Routine has been tested.
      !
      !  Revision 1.1  2004/02/25 14:43:03  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_set_proc_speed_s(proc_speed,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_set_proc_speed_d(proc_speed,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(0:), INTENT(IN   ) :: proc_speed
      INTEGER                , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)               :: t0
      REAL(ppm_kind_double)  :: rsum
      INTEGER                :: i
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_set_proc_speed',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_set_proc_speed',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (SIZE(proc_speed,1) .LT. ppm_nproc) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &            'proc_speed must be at least of length nproc',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that numbers add up to 1 and are >=0
      !-------------------------------------------------------------------------
      rsum = 0.0_ppm_kind_double
      DO i=0,ppm_nproc-1
          IF (proc_speed(i) .LT. 0.0_MK) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &            'proc_speed must be >= 0 for all processors',__LINE__,info)
             GOTO 9999
          ENDIF
          rsum = rsum + REAL(proc_speed(i),ppm_kind_double)
      ENDDO
      IF ((rsum - 1.0_ppm_kind_double) .GT. ppm_myepsd) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_set_proc_speed', &
     &        'proc_speed must sum up to 1.0',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Set the internal ppm_proc_speed values
      !-------------------------------------------------------------------------
      rsum = rsum - REAL(proc_speed(ppm_nproc-1),ppm_kind_double)
      DO i=0,ppm_nproc-2
#if   __KIND == __SINGLE_PRECISION
          ppm_proc_speed(i) = REAL(proc_speed(i),ppm_kind_double)
#elif __KIND == __DOUBLE_PRECISION
          ppm_proc_speed(i) = proc_speed(i)
#endif
      ENDDO
      ppm_proc_speed(ppm_nproc-1) = REAL(1.0_ppm_kind_double - rsum,MK)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_set_proc_speed',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_set_proc_speed_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_set_proc_speed_d
#endif
