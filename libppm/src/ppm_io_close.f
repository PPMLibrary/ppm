      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_close
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine closes an IO channel.
      !
      !  Input        : iUnit          (I) IO unit to be closed.
      !
      !  Input/output :
      !
      !  Output       : info           (I) return code. 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_close.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.6  2004/10/01 16:33:33  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:09:00  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 07:45:26  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.3  2004/07/16 14:46:26  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.2  2004/05/14 14:44:11  ivos
      !  Tested and added argument check for iUnit .LT. 0.
      !
      !  Revision 1.1  2004/05/06 10:44:58  ivos
      !  Preliminary check-in. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_close(iUnit,info)

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
      INTEGER                , INTENT(  OUT) :: info
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
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_close',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF (iUnit .LE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_close',    &
     &        'Unit number needs to be > 0',__LINE__,info)
          GOTO 9999
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
      END SUBROUTINE ppm_io_close
