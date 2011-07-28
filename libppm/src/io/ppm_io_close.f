      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_close
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_close(iUnit,info)
      !!! This routine closes an IO channel.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_io
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(IN   ) :: iUnit
      !!! IO unit to be closed.
      INTEGER                , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lopen
      INTEGER                          :: iou,iopt
      INTEGER, DIMENSION(1)            :: ldc
      CHARACTER(LEN=ppm_char)          :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_close',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that the unit number is actually open
      !-------------------------------------------------------------------------
      lopen = .FALSE.
      IF (ASSOCIATED(ppm_io_unit)) THEN
          IF (SIZE(ppm_io_unit,1) .GE. iUnit) THEN
              IF (ppm_io_unit(iUnit) .GT. 0) lopen = .TRUE.
          ENDIF
      ENDIF
      IF (.NOT. lopen) THEN
         info = ppm_error_notice
         WRITE(mesg,'(A,I4,A)') 'Requested ppm I/O unit ',iUnit,  &
     &        ' is not open. Exiting.'
         CALL ppm_error(ppm_err_no_unit,'ppm_io_close',mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get internal unit number
      !-------------------------------------------------------------------------
      iou = ppm_io_unit(iUnit)

      !-------------------------------------------------------------------------
      !  Savely close the unit
      !-------------------------------------------------------------------------
      IF (ppm_io_mode(iou) .EQ. ppm_param_io_centralized) THEN
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(iUnit,OPENED=lopen)
              IF (.NOT. lopen) THEN
                 info = ppm_error_notice
                 WRITE(mesg,'(A,I4,A)') 'Requested Fortran I/O unit ',iUnit,  &
     &               ' is not open. Exiting.'
                 CALL ppm_error(ppm_err_no_unit,'ppm_io_close',mesg,   &
     &               __LINE__,info)
                 GOTO 9999
              ENDIF
              CLOSE(iUnit,IOSTAT=info)
              IF (info .NE. 0) THEN
                 info = ppm_error_error
                 WRITE(mesg,'(A,I4)') 'Failed to close open I/O unit ',iUnit
                 CALL ppm_error(ppm_err_close,'ppm_io_close',mesg,  &
     &               __LINE__,info)
              ENDIF
          ENDIF
      ELSE
          INQUIRE(iUnit,OPENED=lopen)
          IF (.NOT. lopen) THEN
             info = ppm_error_notice
             WRITE(mesg,'(A,I4,A)') 'Requested Fortran I/O unit ',iUnit,  &
     &           ' is not open. Exiting.'
             CALL ppm_error(ppm_err_no_unit,'ppm_io_close',mesg,__LINE__,info)
             GOTO 9999
          ENDIF
          CLOSE(iUnit,IOSTAT=info)
          IF (info .NE. 0) THEN
             info = ppm_error_error
             WRITE(mesg,'(A,I4)') 'Failed to close open I/O unit ',iUnit
             CALL ppm_error(ppm_err_close,'ppm_io_close',mesg,  &
     &           __LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Debug output
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
          WRITE(mesg,'(A,I4,A)') 'I/O unit ',iUnit,' closed.'
          CALL ppm_write(ppm_rank,'ppm_io_close',mesg,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Mark the unit as unused in the internal list
      !-------------------------------------------------------------------------
      ppm_io_mode(iou)   = ppm_param_undefined
      ppm_io_format(iou) = ppm_param_undefined
      ppm_io_unit(iUnit) = 0

      !-------------------------------------------------------------------------
      !  If no more unit is open, deallocate the lists
      !-------------------------------------------------------------------------
      IF (MAXVAL(ppm_io_unit) .EQ. 0) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_io_close',    &
     &            'No more open units. Deallocating I/O control lists.',info)
          ENDIF
          iopt = ppm_param_dealloc
          CALL ppm_alloc(ppm_io_format,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_close',     &
     &            'I/O format list PPM_IO_FORMAT',__LINE__,info)
          ENDIF
          CALL ppm_alloc(ppm_io_mode,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_close',     &
     &            'I/O mode list PPM_IO_MODE',__LINE__,info)
          ENDIF
          CALL ppm_alloc(ppm_io_unit,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_io_close',     &
     &            'I/O unit list PPM_IO_UNIT',__LINE__,info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_close',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_close',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
        ENDIF
        IF (iUnit .LE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_close',    &
     &        'Unit number needs to be > 0',__LINE__,info)
          GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_io_close
