      !-------------------------------------------------------------------------
      !  Subroutine   :                     substart
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE substart_s(caller,t0,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE substart_d(caller,t0,info)
#endif
      !!! This routine is called whenever a subroutine is started. It
      !!! initializes the info of that subroutine to 0 and prints a debug
      !!! message if ppm_debug is > 1.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(IN   ) :: caller
      !!! Character string with the name of the calling subroutine
      REAL(MK),         INTENT(  OUT) :: t0
      !!! System/cpu time at start of subroutine. Only returned if
      !!! ppm_debug .GT. 0
      INTEGER,          INTENT(  OUT) :: info
      !!! Initialized info for the calling subroutine.

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER :: info2

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      info = ppm_param_success
      IF      (ppm_debug.GT.1) THEN
         CALL ppm_util_time(t0)
         CALL ppm_write(ppm_rank,caller,'entering',info2)
      ELSE IF (ppm_debug.GT.0) THEN
         CALL ppm_util_time(t0)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE substart_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE substart_d
#endif

