      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_util_tic
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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

      SUBROUTINE ppm_util_tic(info)
      !!! Puts the current cpu time in a buffer
      !!! Calls ppm_util_time,
      !!! which uses either `MPI_Wtime`, f90 `CPU_TIME` or `etime`,
      !!! based on how PPM was configured (see
      !!! ./configure --help):
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)               :: ldu
      
      info = 0
      !-------------------------------------------------------------------------
      !  Call ppm_util_time
      !-------------------------------------------------------------------------

      ! increment buffer counter
      ppm_btic%idx = ppm_btic%idx + 1

      ! check that the array in ppm_btic%tic is large enough - reallocate if not
      IF (ppm_btic%nbuff.LT.ppm_btic%idx) THEN
          ldu(1) = ppm_btic%idx + 10
          CALL ppm_alloc(ppm_btic%tic,ldu,&
              ppm_param_alloc_grow_preserve,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_tic','tic',__LINE__,info)
              GOTO 9999
          ENDIF
          ppm_btic%nbuff = ppm_btic%idx + 10
      ENDIF

      ! add a timing entry in the buffer
      CALL ppm_util_time(ppm_btic%tic(ppm_btic%idx))

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_util_tic
