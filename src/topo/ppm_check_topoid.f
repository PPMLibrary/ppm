      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_check_topoid
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

      SUBROUTINE ppm_check_topoid(topoid,valid,info)
      !!! This subroutine is used to check a given topology ID for its validity.
      !!!
      !!! [NOTE]
      !!! It is recommended to not use this routine too often, it should be
      !!! typically placed into DEBUG code and only at critical places in the
      !!! code

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_topo_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN   )        :: topoid
      !!! Topology ID to be checked
      LOGICAL , INTENT(  OUT)        :: valid
      !!! Returns `TRUE` if the given meshid is valid and defined,
      !!! `FALSE` otherwise.
      INTEGER , INTENT(  OUT)        :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_check_topoid',t0,info)
      valid = .FALSE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Validity check
      !-------------------------------------------------------------------------

      IF ((topoid .GE. 1) .AND. (topoid .LE. SIZE(ppm_topo)) .AND. &
     &            (ppm_topo(topoid)%t%isdefined)) THEN
          valid = .TRUE.
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_check_topoid',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_check_topoid',  &
     &            'Please call ppm_init first!',__LINE__,info)
              valid = .FALSE.
              GOTO 8888
          ENDIF
          IF ((topoid.GT.SIZE(ppm_topo)) .OR. &
     &            (topoid.LT.1)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_topoid', &
     &             'topoid indexing outside ppm_topo',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
          IF (.NOT. ASSOCIATED(ppm_topo(topoid)%t)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_topoid', &
     &             'ppm_topo(topoid) pointer not associated',&
     &              __LINE__, info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_check_topoid
