      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_util_time
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
      SUBROUTINE ppm_util_time_s(timing)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_time_d(timing)
#endif
      !!! Returns the current cpu time. Uses `ppm_util_time`,
      !!! which uses either `MPI_Wtime`, f90 `CPU_TIME` or `etime`,
      !!! based on how PPM was configured (see
      !!! ./configure --help):
      !!!
      !!! [NOTE]
      !!! etime is an C intrinsic - thus NOT standard fortran.
      !!! We therefore also write it in small letters.
      !!!
      !!! [WARNING]
      !!! This routine should not be called by the PPM client developer. Use
      !!! `ppm_time` instead.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'ppm_param.h'
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Type kind
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), INTENT(  OUT) :: timing
      !!! Current CPU/wall clock time
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
#ifdef __ETIME
      REAL(MK) :: array(2)
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      REAL(MK), EXTERNAL :: etime
#endif

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Call the MPI function MPI_Wtime
      !-------------------------------------------------------------------------
      timing = MPI_Wtime()
#else
      !-------------------------------------------------------------------------
      !  Call the C routine: etime to get the cpu time
      !-------------------------------------------------------------------------
#ifdef __ETIME
      timing = etime(array)
#else
      CALL CPU_TIME(timing)
#endif
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_time_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_time_d
#endif
