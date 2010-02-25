      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_set_unit
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is called by the user to set the IO
      !                 units for stdout, stderr and log file as used by
      !                 ppm_error, ppm_log and ppm_write.
      !
      !  Input        : stdout     (I) OPTIONAL. Unit number for stdout
      !                 stderr     (I) OPTIONAL. Unit number for stderr
      !                 logfile    (I) OPTIONAL. Unit number for logfile
      !
      !  Input/output : 
      !
      !  Output       : info       (I) OPTIONAL. 0 on success.
      !
      !  Remarks      : Any unit can be given a negative number in which
      !                 case the corresponding output is supressed.
      !
      !                 If units are files, the files are automatically
      !                 opened and REPLACED.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_io_set_unit.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.3  2004/07/26 07:45:27  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.2  2004/06/01 09:28:00  ivos
      !  Arguments chaged back to OPTIONAL again to have consistent interface
      !  within ppm_module_io.
      !
      !  Revision 1.1  2004/05/06 07:43:15  ivos
      !  Was ppm_set_unit before.
      !
      !  Revision 1.5  2004/02/18 17:45:45  walther
      !  Removed optional arguments. The routine is now in ppm_module_io.
      !
      !  Revision 1.4  2004/02/12 17:48:20  ivos
      !  Changed log file from REPLACE to APPEND.
      !
      !  Revision 1.3  2004/01/26 12:35:13  ivos
      !  Updated header.
      !
      !  Revision 1.2  2004/01/13 12:36:45  ivos
      !  Bugfix: forgot an ENDIF.
      !
      !  Revision 1.1  2004/01/13 12:29:26  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_io_set_unit(stdout,stderr,logfile,info)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, OPTIONAL, INTENT(IN   ) :: stdout
      INTEGER, OPTIONAL, INTENT(IN   ) :: stderr
      INTEGER, OPTIONAL, INTENT(IN   ) :: logfile
      INTEGER, OPTIONAL, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL                :: isopen
      INTEGER                :: info2
      CHARACTER(LEN=12)      :: cformat
      CHARACTER(LEN=80)      :: filename
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info2 = 0

      !-------------------------------------------------------------------------
      !  Define the file numbering format
      !-------------------------------------------------------------------------
      IF (ppm_nproc.LE.10) THEN
         cformat = '(A,I1.1,A)' 
      ELSEIF (ppm_nproc.LE.100) THEN
         cformat = '(A,I2.2,A)' 
      ELSEIF (ppm_nproc.LE.1000) THEN
         cformat = '(A,I3.3,A)' 
      ELSE
         cformat = '(A,I4.4,A)' 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Test if all units are open and open them if needed
      !-------------------------------------------------------------------------
      IF (PRESENT(stdout)) THEN
          INQUIRE(stdout,OPENED=isopen)
          IF (.NOT. isopen) THEN
              ! just in case all processors share the disk...
              WRITE(filename,cformat) 'ppm_stdout_',ppm_rank,'.out'
              OPEN(stdout,FILE=filename,   &
     &        STATUS='REPLACE',ACTION='WRITE',IOSTAT=info2)
          ENDIF
          ppm_stdout = stdout
      ENDIF

      IF (PRESENT(stderr)) THEN
          INQUIRE(stderr,OPENED=isopen)
          IF (.NOT. isopen) THEN
              ! just in case all processors share the disk...
              WRITE(filename,cformat) 'ppm_stderr_',ppm_rank,'.out'
              OPEN(stderr,FILE=filename,   &
     &        STATUS='REPLACE',ACTION='WRITE',IOSTAT=info2)
          ENDIF
          ppm_stderr = stderr
      ENDIF

      IF (PRESENT(logfile)) THEN
          INQUIRE(logfile,OPENED=isopen)
          IF (.NOT. isopen) THEN
              ! just in case all processors share the disk...
              WRITE(filename,cformat) 'ppm_log_',ppm_rank,'.out'
              OPEN(logfile,FILE=filename,  &
     &        POSITION='APPEND',ACTION='WRITE',IOSTAT=info2)
          ENDIF
          ppm_logfile = logfile
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      IF (PRESENT(info)) info = info2
      RETURN
      END SUBROUTINE ppm_io_set_unit
