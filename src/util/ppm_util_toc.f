      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_util_toc
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

      SUBROUTINE ppm_tstats_toc(id,step,diff_t,info,verbose)
      !!! Calls ppm_util_time, pop the last tic out of the buffer
      !!! and returns the difference.
      !!! Optionally, print the results on stdout.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_util_time
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,            INTENT(IN   ) :: id
      INTEGER,            INTENT(IN   ) :: step
      REAL(MK),           INTENT(  OUT) :: diff_t
      !!! Difference in time between tic and toc
      INTEGER,            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      LOGICAL, OPTIONAL,  INTENT(IN   ) :: verbose
      !!! Difference in time between tic and toc
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK) :: t0,t1

      INTEGER, DIMENSION(3) :: ldu

      !!! Current CPU clock time
      CHARACTER(LEN=ppm_char) :: caller = 'ppm_tstats_toc'

      info = 0
      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------
      IF (ppm_tstats(id)%times(step).EQ.0.0_MK) THEN
         fail("never has been ticked before",ppm_err_sub_failed,ppm_error=ppm_error_fatal)
      ENDIF
      CALL ppm_util_time(t1)

      ! difference between current time and last entry in the tic buffer
      diff_t = t1-ppm_tstats(id)%times(step)
      ppm_tstats(id)%times(step) = diff_t

      IF (PRESENT(verbose)) THEN
         IF (verbose) THEN
            stdout_f('(A,A,E17.7,A)','ppm_tstats(id)%label'," took ",diff_t," seconds")
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_tstats_toc
