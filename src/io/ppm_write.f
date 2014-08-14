      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_write
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

      SUBROUTINE ppm_write(rank,caller,cbuf,info,iUnit)
      !!! This subroutine enables the user to write from any processor a
      !!! character strings to stdout or a given I/O Unit.
      !!!
      !!! This routine uses the `WRITE` intrinsic routine.

      USE ppm_module_data, ONLY : ppm_stdout, ppm_char
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,           INTENT(IN   ) :: rank
      !!! MPI rank of the calling processor
      CHARACTER(LEN=*),  INTENT(IN   ) :: caller
      !!! Character string describing the name of the calling subroutine
      CHARACTER(LEN=*),  INTENT(IN   ) :: cbuf
      !!! Character string containing the message to be printed
      INTEGER,           INTENT(  OUT) :: info
      !!! Returns 0 on success
      INTEGER, OPTIONAL, INTENT(IN   ) :: iUnit
      !!! UNIT to print to (OPTIONAL). Defaults to stdout if not specified.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER :: iu
      INTEGER :: ios
      INTEGER :: icaller,ibuf

      CHARACTER(LEN=ppm_char) :: cformat
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      ! print to stdout by default
      iu = ppm_stdout
      IF (PRESENT(iUnit)) THEN
         IF (iUnit .GT. 0) iu = iUnit
      ENDIF

      !-------------------------------------------------------------------------
      !  Get length of messages
      !-------------------------------------------------------------------------
      icaller = LEN_TRIM(caller)
      ibuf    = LEN_TRIM(cbuf)

      !-------------------------------------------------------------------------
      !  Define the print format
      !-------------------------------------------------------------------------
      IF     (rank.LT.0) THEN
         cformat = '(4A)'
         !----------------------------------------------------------------------
         !  Do the print
         !----------------------------------------------------------------------
         IF (iu .GE. 0) THEN
            WRITE(iu,cformat,IOSTAT=ios)'(', caller(1:icaller), ') : ', cbuf(1:ibuf)
         ENDIF
      ELSE
         cformat = '(A,I0,4A)'
         !----------------------------------------------------------------------
         !  Do the print
         !----------------------------------------------------------------------
         IF (iu .GE. 0) THEN
            WRITE(iu,cformat,IOSTAT=ios)'[',rank,'] (', caller(1:icaller), ') : ', cbuf(1:ibuf)
         ENDIF
      ENDIF

      info = ios

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_write
