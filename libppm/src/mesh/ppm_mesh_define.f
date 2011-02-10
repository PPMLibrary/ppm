      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_define
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

      SUBROUTINE ppm_mesh_define(topoid,meshid,Nm,istart,ndata,info)
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

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID for which to create mesh
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Number of mesh POINTS in each dimension. Subs must be compatible
      !!! with this mesh, otherwise an error occurs.
      INTEGER                 , INTENT(INOUT) :: meshid
      !!! Mesh ID of the new mesh. If .LE. 0 on input,
      !!! the routine will create an automatic one and return it here.
      INTEGER , DIMENSION(:,:), POINTER       :: istart
      !!! Start indices of all subs meshes in global mesh
      INTEGER , DIMENSION(:,:), POINTER       :: ndata
      !!! Number of mesh points in each direction on each sub mesh
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)      :: t0
      LOGICAL                    :: valid
      TYPE(ppm_t_topo), POINTER  :: topo => NULL()

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_define',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Create new mesh
      !-------------------------------------------------------------------------
      IF (topo%prec .EQ. ppm_kind_single) THEN
        CALL ppm_mesh_on_subs(Nm,topo%min_physs,topo%max_physs,topo%min_subs, &
     &      topo%max_subs,topo%nsubs,istart,ndata,info)
      ELSE
        CALL ppm_mesh_on_subs(Nm,topo%min_physd,topo%max_physd,topo%min_subd, &
     &      topo%max_subd,topo%nsubs,istart,ndata,info)
      ENDIF

      IF (info .NE. 0) THEN
          info = ppm_error_error
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_define',   &
     &        'Storing new mesh failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_define',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_mesh_define',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_define',  &
     &            'topoid is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_define
