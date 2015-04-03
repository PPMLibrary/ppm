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
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_topo_typedef
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: meshid
      !!! Mesh ID to be checked.
      INTEGER, INTENT(IN   ) :: topoid
      !!! Topology ID on which the mesh is defined.
      LOGICAL, INTENT(  OUT) :: valid
      !!! Returns `TRUE` if the given meshid is valid and defined,
      !!! `FALSE` otherwise.
      INTEGER, INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_), POINTER :: mesh

      REAL(ppm_kind_double) :: t0

      INTEGER :: meshid_

      CHARACTER(LEN=ppm_char) :: caller="ppm_check_meshid"

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      meshid_=-1
      mesh => ppm_mesh%begin()
      DO WHILE (ASSOCIATED(mesh))
         IF (mesh%ID.EQ.meshid) THEN
            IF (mesh%topoid.EQ.topoid) THEN
               meshid_=ppm_mesh%iter_id
               EXIT
            ENDIF
         ENDIF
         mesh => ppm_mesh%next()
      ENDDO

      IF (meshid_.LT.ppm_mesh%min_id.OR.meshid_.GT.ppm_mesh%max_id) THEN
         stdout("Input mesh ID ",meshid," Is Invalid!")
         valid=.FALSE.
      ELSE
         valid=.TRUE.
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT.ppm_initialized) THEN
             valid=.FALSE.
             fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_check_meshid
