      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_error
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is called whenever an error occurs. It
      !                 prints the error message, writes a log entry and
      !                 terminates the program if the error was fatal.
      !
      !  Input        : errno      (I) error number
      !                 caller     (C) name of calling subroutine
      !                 mesg       (C) error message
      !                 line       (I) line number of error
      !
      !  Input/output : info       (I) On entry: error severity level.
      !                                   ppm_error_fatal
      !                                   ppm_error_error
      !                                   ppm_error_warning
      !                                   ppm_error_notice
      !                                On exit: new error level (maybe
      !                                different from entry level if this
      !                                routine had additional errors.) 
      !
      !  Output       : 
      !
      !  Remarks      : Line number is only included in error message if
      !                 debug level .GT. 0.
      !
      !                 Warnings and Notices are not printed to stderr (but
      !                 only to the log file) unless debug level is .GT. 0
      !                 or .GT. 1, respectively.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_error.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.13  2004/07/26 07:45:25  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.12  2004/06/10 16:20:00  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.11  2004/05/13 17:12:31  ivos
      !  bigfix: added TRIM() around mesg to avoid internal file overflows.
      !
      !  Revision 1.10  2004/01/23 18:12:07  ivos
      !  Well, some bug fix... Guess I should go home.
      !
      !  Revision 1.9  2004/01/23 18:09:39  ivos
      !  Cosmetics to the error messages.
      !
      !  Revision 1.8  2004/01/23 17:24:15  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.7  2004/01/22 13:25:52  ivos
      !  Included line number in error messages of debug level .GT. 0
      !
      !  Revision 1.6  2004/01/21 16:58:15  ivos
      !  Now only writes to stderr and logfile (and never to stdout). Stderr
      !  output further depends on the debug level that is set.
      !
      !  Revision 1.5  2004/01/19 11:49:35  ivos
      !  Changed to use ppm_error.h.
      !
      !  Revision 1.4  2004/01/13 12:27:59  ivos
      !  Added use of log file and stderr.
      !
      !  Revision 1.3  2004/01/12 12:40:36  ivos
      !  Bugfix: added ifdef __MPI around call to MPI_Abort.
      !
      !  Revision 1.2  2004/01/09 15:56:29  ivos
      !  Added error level constants from ppm_param.h instead of hard numbers.
      !
      !  Revision 1.1  2004/01/09 15:29:46  ivos
      !  Initial implementation.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_error(errno,caller,mesg,line,info)
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_log
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: errno,line
      CHARACTER(LEN=*)        , INTENT(IN   ) :: caller,mesg
      INTEGER                 , INTENT(INOUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                                 :: info2
      CHARACTER(LEN=ppm_char)                 :: msg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info2 = 0

      !-------------------------------------------------------------------------
      !  Assemble error message
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (info .EQ. ppm_error_fatal) THEN
              ! fatal error
              WRITE(msg,'(I4.4,A,I5.5,4A)') errno,' <FATAL> at line ', &
         &            line,'. ',TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_error) THEN
              ! error
              WRITE(msg,'(I4.4,A,I5.5,4A)') errno,' <ERROR> at line ', &
         &            line,'. ',TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_warning) THEN
              ! warning
              WRITE(msg,'(I4.4,A,I5.5,4A)') errno,' <WARNING> at line ',&
         &            line,'. ',TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_notice) THEN
              ! notice
              WRITE(msg,'(I4.4,A,I5.5,4A)') errno,' <NOTICE> at line ',&
         &            line,'. ',TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSE
              ! unknown
              WRITE(msg,'(I4.4,A,I5.5,4A)') errno,' <UNKNOWN> at line ',&
         &            line,'. ',TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ENDIF
      ELSE
          IF (info .EQ. ppm_error_fatal) THEN
              ! fatal error
              WRITE(msg,'(I4.4,4A)') errno,' <FATAL> ',                &
         &            TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_error) THEN
              ! error
              WRITE(msg,'(I4.4,4A)') errno,' <ERROR> ',                &
         &            TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_warning) THEN
              ! warning
              WRITE(msg,'(I4.4,4A)') errno,' <WARNING> ',              &
         &            TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSEIF (info .EQ. ppm_error_notice) THEN
              ! notice
              WRITE(msg,'(I4.4,4A)') errno,' <NOTICE> ',               &
         &            TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ELSE
              ! unknown
              WRITE(msg,'(I4.4,4A)') errno,' <UNKNOWN> ',              &
         &            TRIM(ppm_err_mesg(errno)), ': ',TRIM(mesg)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Write the error to stderr if debug level is high enough
      !-------------------------------------------------------------------------
      IF (info .EQ. ppm_error_notice) THEN
          IF (ppm_debug .GT. 1) THEN
              CALL ppm_write(ppm_rank,caller,msg,info2,ppm_stderr)
          ENDIF
      ELSEIF (info .EQ. ppm_error_warning) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,caller,msg,info2,ppm_stderr)
          ENDIF
      ELSE
          CALL ppm_write(ppm_rank,caller,msg,info2,ppm_stderr)
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Always write the error to log file
      !-------------------------------------------------------------------------
      CALL ppm_log(caller,msg,info2)
      
      !-------------------------------------------------------------------------
      !  If error was fatal: terminate program. All processors get killed
      !  by MPI.
      !-------------------------------------------------------------------------
      IF (info .EQ. ppm_error_fatal) THEN
#ifdef __MPI
          CALL MPI_Abort(ppm_comm,1,info2)
#endif
          STOP
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine new error level
      !-------------------------------------------------------------------------
      IF (info2 .GT. info .AND. info2 .LT. 0) info = info2

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_error
