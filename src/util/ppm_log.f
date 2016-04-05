      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_log
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

      SUBROUTINE ppm_log(caller,mesg,info)
      !!! This routine writes a log message to the ppm log file
      !!! (as defined using ppm_set_unit)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data, ONLY : ppm_logfile, ppm_rank, ppm_char
      USE ppm_module_write
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN   ) :: caller
      !!! Name of calling subroutine
      CHARACTER(LEN=*), INTENT(IN   ) :: mesg
      !!! Log message
      INTEGER,          INTENT(INOUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(8) :: tval

      CHARACTER(LEN=ppm_char) :: msg

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      info = 0
      !-------------------------------------------------------------------------
      !  Get time and date and append to message
      !-------------------------------------------------------------------------
      CALL DATE_AND_TIME(values=tval)
      WRITE(msg,'(I4.4,A,5(I2.2,A),A)') tval(1),'-',tval(2),'-',tval(3),   &
      & ' ',tval(5),':',tval(6),':',tval(7),' --- ',TRIM(mesg)

      !-------------------------------------------------------------------------
      !  Write the log
      !-------------------------------------------------------------------------
      IF (ppm_logfile.GT.0) THEN
         CALL ppm_write(ppm_rank,caller,msg,info,ppm_logfile)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_log
