      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_log
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine writes a log message to the ppm log
      !                 file (as defined using ppm_set_unit)
      !
      !  Input        : caller     (C) name of calling subroutine
      !                 mesg       (C) log message
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
      !  $Log: ppm_log.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:55  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.5  2004/07/26 07:45:28  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/02/18 14:26:26  ivos
      !  bugfix: now using local msg to write stuff to since mesg is INTENT(IN).
      !
      !  Revision 1.3  2004/02/12 17:47:35  ivos
      !  Now only calls ppm_write if any log file is set.
      !
      !  Revision 1.2  2004/01/22 13:27:56  ivos
      !  Did (1) update of the header, (2) replaced pwrite with ppm_write or
      !  ppm_error calls, (3) inserted validity check of INTENT(IN) arguments
      !  where needed.
      !
      !  Revision 1.1  2004/01/13 12:31:15  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_log(caller,mesg,info)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data, ONLY : ppm_logfile, ppm_rank, ppm_char
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*)        , INTENT(IN   ) :: caller,mesg
      INTEGER                 , INTENT(INOUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(8)                   :: tval
      CHARACTER(LEN=ppm_char)                 :: msg
      
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info = 0
      !-------------------------------------------------------------------------
      !  Get time and date and prepend to message
      !-------------------------------------------------------------------------
      CALL DATE_AND_TIME(values=tval)
      WRITE(msg,'(I4.4,A,5(I2.2,A),A)') tval(1),'-',tval(2),'-',tval(3),   &
     &               ' ',tval(5),':',tval(6),':',tval(7),' --- ',trim(mesg)

      !-------------------------------------------------------------------------
      !  Write the log
      !-------------------------------------------------------------------------
      IF (ppm_logfile .GT. 0) THEN
          CALL ppm_write(ppm_rank,caller,msg,info,ppm_logfile)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_log
