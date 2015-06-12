      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_define
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

      SUBROUTINE ppm_mesh_define(topoid,meshid,Nm,info,Offset,ghostsize)
      !!! This routine defines a new mesh on an existing topology.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_topo_typedef
      USE ppm_module_mesh_typedef
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,                            INTENT(IN   ) :: topoid
      !!! Topology ID for which to create mesh
      INTEGER,                            INTENT(INOUT) :: meshid
      !!! Mesh ID of the new mesh. If .LE. 0 on input,
      !!! the routine will create an automatic one and return it here.
      !!! If .GT. 0 it will return an Error
      INTEGER,  DIMENSION(:  ),           INTENT(IN   ) :: Nm
      !!! Number of mesh POINTS in each dimension. Subs must be compatible
      !!! with this mesh, otherwise an error occurs.
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: Offset
      !!! Offset in each dimension
      INTEGER,  DIMENSION(:),   OPTIONAL, INTENT(IN   ) :: ghostsize
      !!! size of the ghost layer, in number of mesh nodes for each dimension

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_), POINTER :: mesh

      REAL(ppm_kind_double)                       :: t0
      REAL(ppm_kind_double), DIMENSION(1:ppm_dim) :: Offst

      INTEGER, DIMENSION(1:ppm_dim) :: ighostsize

      CHARACTER(LEN=ppm_char) :: caller='ppm_mesh_define'

      LOGICAL :: valid
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      IF (PRESENT(Offset)) THEN
         Offst=REAL(Offset,ppm_kind_double)
      ELSE
         Offst=0.0_ppm_kind_double
      ENDIF

      IF (PRESENT(ghostsize)) THEN
         ighostsize(1:ppm_dim)=ghostsize(1:ppm_dim)
      ELSE
         ighostsize=0
      ENDIF

      IF (meshid.LE.0) THEN
         !-----------------------------------------------------------------------
         !  Create a new mesh and mesh identifier if the user specified none
         !-----------------------------------------------------------------------
         ALLOCATE(ppm_t_equi_mesh::mesh,STAT=info)
         or_fail_alloc("mesh pointer")

         CALL mesh%create(topoid,Offst,info,Nm=Nm,ghostsize=ighostsize)
         or_fail("Failed to create mesh!")

         !-----------------------------------------------------------------------
         !  Created mesh and mesh identifier are stored in the collection
         !-----------------------------------------------------------------------
         meshid=mesh%ID

         mesh => NULL()
      ELSE
         stdout("meshid",meshid," is invalid! It should be .LE. 0 on input")
         fail(cbuf,ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
             fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
             fail('topoid is invalid!',ppm_err_argument,exit_point=8888)
          ENDIF
          IF (ppm_topo(topoid)%t%nsubs.LE.0) THEN
             fail('nsubs must be >0',exit_point=8888)
          ENDIF
          IF (ANY(Nm(1:ppm_dim).LT.2)) THEN
             fail('Nm must be >1 in all space dimensions',exit_point=8888)
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_define
