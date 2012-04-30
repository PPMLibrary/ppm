      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_error
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License 
      ! as published by the Free Software Foundation, either 
      ! version 3 of the License, or (at your option) any later 
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_error(errno,caller,mesg,line,info)
      !!! This routine is called whenever an error occurs. It prints the error
      !!! message, writes a log entry and terminates the program if the error
      !!! was fatal.
      !!!
      !!! [NOTE]
      !!! Line number is only included in error message if debug level > 0.
      !!! Warnings and Notices are not printed to stderr (but only to the log
      !!! file) unless debug level is > 0 or > 1, respectively.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_log
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: errno
      !!! Error number
      INTEGER                 , INTENT(IN   ) :: line
      !!! Line number of error
      CHARACTER(LEN=*)        , INTENT(IN   ) :: caller
      !!! Name of calling subroutine
      CHARACTER(LEN=*)        , INTENT(IN   ) :: mesg
      !!! Error message
      INTEGER                 , INTENT(INOUT) :: info
      !!! On entry: error severity level.
      !!!
      !!! * ppm_error_fatal
      !!! * ppm_error_error
      !!! * ppm_error_warning
      !!! * ppm_error_notice
      !!!
      !!! On exit: new error level (maybe different from entry level if this
      !!! routine had additional errors.)
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
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_error
