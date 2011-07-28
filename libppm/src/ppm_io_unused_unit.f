      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_unused_unit
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine finds the next available unit number 
      !                 starting from the given one and counting upwards.
      !
      !  Input        : 
      !
      !  Input/output : iUnit          (I) on input: starting number. On
      !                                    output: next free unit number.
      !                                    -1 is returned if no free unit
      !                                    number .LE. 255 was found.
      !
      !  Output       : info           (I) return code. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_unused_unit.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/10/01 16:33:34  ivos
      !  cosmetics.
      !
      !  Revision 1.3  2004/10/01 16:09:02  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.2  2004/07/26 07:45:27  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.1  2004/05/13 11:38:30  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_unused_unit(iUnit,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER             , INTENT(INOUT) :: iUnit
      INTEGER             , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lopen
      INTEGER                          :: localunit
      INTEGER, PARAMETER               :: maxUnit = 255
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_unused_unit',t0,info)

      !-------------------------------------------------------------------------
      !  Check start unit number
      !-------------------------------------------------------------------------
      IF (iUnit .LE. 0) THEN
          iUnit = 10
          info = ppm_error_notice
          CALL ppm_error(ppm_err_argument,'ppm_io_unused_unit',    &
     &        'Unit number out of range was reset to 10.',__LINE__,info)
      ENDIF
      IF (iUnit .GT. maxUnit) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_unused_unit',    &
     &        'Unit number out of range.',__LINE__,info)
          iUnit = -1
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find next free unit number
      !-------------------------------------------------------------------------
      localunit = iUnit
      lopen = .TRUE.
      DO WHILE (lopen .AND. localunit .LE. maxUnit)
          localunit = localunit + 1
          INQUIRE(localunit,OPENED=lopen)
          IF (localunit .EQ. ppm_stdout) lopen = .TRUE.
          IF (localunit .EQ. ppm_stderr) lopen = .TRUE.
          IF (localunit .EQ. ppm_logfile) lopen = .TRUE.
      ENDDO 

      !-------------------------------------------------------------------------
      !  Check if one was found
      !-------------------------------------------------------------------------
      IF (lopen) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_outof_units,'ppm_io_unused_unit',    &
     &        'No unused I/O unit was found.',__LINE__,info)
          iUnit = -1
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return the found unit
      !-------------------------------------------------------------------------
      iUnit = localunit

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_unused_unit',t0,info)
      RETURN
      END SUBROUTINE ppm_io_unused_unit
