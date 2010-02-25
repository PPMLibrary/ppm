      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_io_delete
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine deletes a file.
      !
      !  Input        : filename       (C) Name of file to be deleted
      !                 mode           (I) I/O mode. One of:
      !                                       ppm_param_io_distributed
      !                                       ppm_param_io_centralized
      !                                    In the distributed mode, every
      !                                    processor deletes the file
      !                                    locally. In the centralized
      !                                    mode, only processor 0 deletes
      !                                    the file.
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
      !  $Log: ppm_io_delete.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2006/09/04 18:34:49  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.7  2004/10/01 16:33:33  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:09:00  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/26 07:45:26  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/07/16 14:46:26  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/05/14 14:48:19  ivos
      !  Updated header comment.
      !
      !  Revision 1.2  2004/05/13 11:41:33  ivos
      !  Changed to use ppm_io_unused_unit.
      !
      !  Revision 1.1  2004/05/11 14:54:27  ivos
      !  Initial implementation. Not tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_delete(filename,mode,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_io_unused_unit
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)    , INTENT(IN   ) :: filename
      INTEGER             , INTENT(IN   ) :: mode
      INTEGER             , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      LOGICAL                          :: lexists,lopen
      INTEGER                          :: iUnit
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_io_delete',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_io_delete',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF ((mode .NE. ppm_param_io_distributed) .AND.        &
     &    (mode .NE. ppm_param_io_centralized)) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_io_delete',     &
     &        'Invalid mode specified.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that file exists
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          INQUIRE(FILE=filename,EXIST=lexists)
          IF (.NOT. lexists) THEN
             info = ppm_error_notice
             CALL ppm_error(ppm_err_file,'ppm_io_delete',    &
     &           TRIM(filename),__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(FILE=filename,EXIST=lexists)
              IF (.NOT. lexists) THEN
                 info = ppm_error_notice
                 CALL ppm_error(ppm_err_file,'ppm_io_delete',    &
     &               TRIM(filename),__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that file is not connected to a unit
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          INQUIRE(FILE=filename,OPENED=lopen)
          IF (lopen) THEN
             info = ppm_error_warning
             CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &           'Cannot delete an open file',__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              INQUIRE(FILE=filename,OPENED=lopen)
              IF (lopen) THEN
                 info = ppm_error_warning
                 CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &               'Cannot delete an open file',__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Attempt delete
      !-------------------------------------------------------------------------
      IF (mode .EQ. ppm_param_io_distributed) THEN
          !---------------------------------------------------------------------
          !  Find unused unit number
          !---------------------------------------------------------------------
          iUnit = 9
          CALL ppm_io_unused_unit(iUnit,info)
          IF (iUnit .LT. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_open,'ppm_io_delete',    &
     &              'No I/O unit available.',__LINE__,info)
             GOTO 9999
          ENDIF
          !---------------------------------------------------------------------
          !  Delete
          !---------------------------------------------------------------------
          OPEN(iUnit,FILE=filename,STATUS='OLD',ACTION='WRITE',IOSTAT=info)
          CLOSE(iUnit,STATUS='DELETE',IOSTAT=info)
          IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &           TRIM(filename),__LINE__,info)
             GOTO 9999
          ENDIF
      ELSE
          IF (ppm_rank .EQ. 0) THEN
              !-----------------------------------------------------------------
              !  Find unused unit number
              !-----------------------------------------------------------------
              iUnit = 9
              CALL ppm_io_unused_unit(iUnit,info)
              IF (iUnit .LT. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_open,'ppm_io_delete',    &
     &              'No I/O unit available.',__LINE__,info)
                 GOTO 9999
              ENDIF
              !-----------------------------------------------------------------
              !  Delete
              !-----------------------------------------------------------------
              OPEN(iUnit,FILE=filename,STATUS='OLD',ACTION='WRITE',IOSTAT=info)
              CLOSE(iUnit,STATUS='DELETE',IOSTAT=info)
              IF (info .NE. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_delete,'ppm_io_delete',    &
     &               TRIM(filename),__LINE__,info)
                 GOTO 9999
              ENDIF
          ENDIF
      ENDIF
       
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_io_delete',t0,info)
      RETURN
      END SUBROUTINE ppm_io_delete
