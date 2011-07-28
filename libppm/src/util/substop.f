      !-------------------------------------------------------------------------
      !  Subroutine   :                     substop
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Lab (ETH Zurich), 
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
      SUBROUTINE substop_s(caller,t0,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE substop_d(caller,t0,info)
#endif
      !!! This routine is called at the end of each subroutine.Depending on
      !!! the debug level, it prints the cpu time, and calling routine.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_util_time
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
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
      REAL(MK)        , INTENT(IN   ) :: t0
      !!! System/cpu time at start of subroutine. (as returned by substart)
      INTEGER         , INTENT(IN   ) :: info
      !!! The final info of the calling subroutine. Is printed if debug.GT.1.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                       :: t1
      CHARACTER(LEN=ppm_char)        :: cbuf
      INTEGER                        :: info2 = 0

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  In parallel collect the MIN info of the processors
      !  THIS IS COMMENTED OUT AS IT IS VERY EXPENSIVE.
      !-------------------------------------------------------------------------
!     CALL MPI_AllReduce(info,i,1,MPI_INTEGER,MPI_MIN,comm,info2)
!     info = i
#endif
      !-------------------------------------------------------------------------
      !  Using debugging print the 'leaving <routine>'
      !-------------------------------------------------------------------------
      CALL ppm_util_time(t1)
      IF     (ppm_debug.GT.1) THEN
         WRITE(cbuf,'(A,I2,A,E12.4)') 'leaving with info=',info,'. took: ',t1-t0
         CALL ppm_write(ppm_rank,caller,cbuf,info2)
      ELSEIF (ppm_debug.GT.0) THEN
         WRITE(cbuf,'(A,E12.4)') 'took: ',t1-t0
         CALL ppm_write(ppm_rank,caller,cbuf,info2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE substop_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE substop_d
#endif
