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

      SUBROUTINE ppm_mesh_define(topoid,meshid,Nm,istart,ndata,info,Offset,ghostsize)
      !!! This routine defines a new mesh on an existing topology.
      !!! The routine checks that the subdomains of the topology are
      !!! compatible with the specified mesh (i.e. are integer multiples
      !!! of the mesh spacing in extent). If not, an error is returned.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_mesh_on_subs
      USE ppm_module_mesh_store
      USE ppm_module_check_id
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
      INTEGER,  DIMENSION(:  ),           INTENT(IN   ) :: Nm
      !!! Number of mesh POINTS in each dimension. Subs must be compatible
      !!! with this mesh, otherwise an error occurs.
      INTEGER,  DIMENSION(:,:),           POINTER       :: istart
      !!! Start indices of all subs meshes in global mesh
      INTEGER,  DIMENSION(:,:),           POINTER       :: ndata
      !!! Number of mesh points in each direction on each sub mesh
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: Offset
      !!! Offset in each dimension
      INTEGER,  DIMENSION(:),   OPTIONAL, INTENT(IN   ) :: ghostsize
      !!! size of the ghost layer, in number of mesh nodes for each dimension

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      REAL(MK)                     :: t0
      REAL(MK), DIMENSION(ppm_dim) :: Offst

      INTEGER, DIMENSION(ppm_dim) :: ighostsize

      CHARACTER(LEN=ppm_char) :: caller='ppm_mesh_define'

      LOGICAL :: valid

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      IF (PRESENT(Offset)) THEN
         Offst=Offset
      ELSE
         Offst=0.0_MK
      ENDIF
      IF (PRESENT(ghostsize)) THEN
         ighostsize(1:ppm_dim)=ghostsize(1:ppm_dim)
      ELSE
         ighostsize=0
      ENDIF
      !-------------------------------------------------------------------------
      !  Create new mesh
      !-------------------------------------------------------------------------
      SELECT CASE (topo%prec)
      CASE (ppm_kind_single)
         CALL ppm_mesh_on_subs(Nm,topo%min_physs,topo%max_physs,topo%min_subs, &
         &    topo%max_subs,topo%nsubs,istart,ndata,info,Offset=Offst)

      CASE DEFAULT
         CALL ppm_mesh_on_subs(Nm,topo%min_physd,topo%max_physd,topo%min_subd, &
         &    topo%max_subd,topo%nsubs,istart,ndata,info,Offset=Offst)

      END SELECT
      or_fail('Defining meshes failed')

      !-------------------------------------------------------------------------
      !  Store new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info,&
      &    Offset=Offst,ghostsize=ighostsize)
      or_fail('Storing new mesh failed.')

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,caller,  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'topoid is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_define
