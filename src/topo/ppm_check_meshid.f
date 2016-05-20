      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_check_meshid
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

      SUBROUTINE ppm_check_meshid(topoid,meshid,valid,info)
      !!! This subroutine is used to check a given mesh ID for its validity.
      !!!
      !!! [NOTE]
      !!! It is recommended to not use this routine too often, it should be
      !!! typically placed into DEBUG code and only at critical places in the
      !!! code
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER , INTENT(IN   )        :: meshid
      !!! Mesh ID to be checked.
      INTEGER , INTENT(IN   )        :: topoid
      !!! Topology ID on which the mesh is defined.
      LOGICAL , INTENT(  OUT)        :: valid
      !!! Returns `TRUE` if the given meshid is valid and defined,
      !!! `FALSE` otherwise.
      INTEGER , INTENT(  OUT)        :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0
      LOGICAL                        :: topo_ok
      TYPE(ppm_t_topo), POINTER      :: topo

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_check_meshid',t0,info)
      valid = .TRUE.

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
      CALL ppm_check_topoid(topoid,topo_ok,info)
      IF (.NOT. topo_ok) THEN
        valid = .FALSE.
        GOTO 9999
      ENDIF
      topo => ppm_topo(topoid)%t

      IF ((meshid .LT. 1).OR.(meshid .GT. topo%max_meshid)) THEN
          valid = .FALSE.
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_check_meshid',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_check_meshid',  &
     &            'Please call ppm_init first!',__LINE__,info)
              valid = .FALSE.
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,topo_ok,info)
          IF (.NOT. topo_ok) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_check_meshid',  &
     &            'Topoid out of range!',__LINE__,info)
              valid = .FALSE.
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_check_meshid
