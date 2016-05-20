      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_io_set_unit
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

      SUBROUTINE ppm_io_set_unit(stdout,stderr,logfile,info)
      !!! This routine is called by the user to set the IO
      !!! units for stdout, stderr and log file as used by
      !!! `ppm_error`, `ppm_log` and `ppm_write`.
      !!!
      !!! [NOTE]
      !!! Any unit can be given a negative number in which case the corresponding
      !!! output is supressed.
      !!!
      !!! [CAUTION]
      !!! If units are files, the files are automatically opened and *replaced*.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, OPTIONAL, INTENT(IN   ) :: stdout
      !!! Unit number for stdout
      INTEGER, OPTIONAL, INTENT(IN   ) :: stderr
      !!! Unit number for stderr
      INTEGER, OPTIONAL, INTENT(IN   ) :: logfile
      !!! Unit number for logfile
      INTEGER, OPTIONAL, INTENT(  OUT) :: info
      !!! Return status, 0 on success
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
              OPEN(stdout,FILE=filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=info2)
          ENDIF
          ppm_stdout = stdout
      ENDIF

      IF (PRESENT(stderr)) THEN
          INQUIRE(stderr,OPENED=isopen)
          IF (.NOT. isopen) THEN
              ! just in case all processors share the disk...
              WRITE(filename,cformat) 'ppm_stderr_',ppm_rank,'.out'
              OPEN(stderr,FILE=filename,STATUS='REPLACE',ACTION='WRITE',IOSTAT=info2)
          ENDIF
          ppm_stderr = stderr
      ENDIF

      IF (PRESENT(logfile)) THEN
         ! logfile = -1, on GNU compiler
         ! Fortran runtime error: Inquire statement identifies an internal file
         IF (logfile.NE.-1) THEN
            INQUIRE(logfile,OPENED=isopen)
            IF (.NOT.isopen) THEN
               ! just in case all processors share the disk...
               WRITE(filename,cformat) 'ppm_log_',ppm_rank,'.out'
               OPEN(logfile,FILE=filename,POSITION='APPEND',ACTION='WRITE',IOSTAT=info2)
            ENDIF
            ppm_logfile=logfile
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      IF (PRESENT(info)) info = info2
      RETURN
      END SUBROUTINE ppm_io_set_unit
