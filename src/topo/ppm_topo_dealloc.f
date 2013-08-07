      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_topo_dealloc
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

      SUBROUTINE ppm_topo_dealloc(topo,info)
      !!! This routine deallocates all members of a topology object.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), INTENT(INOUT) :: topo
      !!! Topology structure to be deallocated
      INTEGER,          INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3) :: ldc
      INTEGER                :: iopt
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_dealloc',t0,info)
      ldc  = 1
      iopt = ppm_param_dealloc

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Deallocate all array members
      !-------------------------------------------------------------------------
      CALL ppm_alloc(topo%min_subs,ldc,iopt,info)
      CALL ppm_alloc(topo%max_subs,ldc,iopt,info)
      CALL ppm_alloc(topo%min_subd,ldc,iopt,info)
      CALL ppm_alloc(topo%max_subd,ldc,iopt,info)
      CALL ppm_alloc(topo%sub2proc,ldc,iopt,info)
      CALL ppm_alloc(topo%isublist,ldc,iopt,info)
      CALL ppm_alloc(topo%bcdef,ldc,iopt,info)
      CALL ppm_alloc(topo%sub_costs,ldc,iopt,info)
      CALL ppm_alloc(topo%sub_costd,ldc,iopt,info)
      CALL ppm_alloc(topo%nneighsubs,ldc,iopt,info)
      CALL ppm_alloc(topo%ineighsubs,ldc,iopt,info)
      CALL ppm_alloc(topo%subs_bc,ldc,iopt,info)
      CALL ppm_alloc(topo%ineighproc,ldc,iopt,info)
      CALL ppm_alloc(topo%icommseq,ldc,iopt,info)
      CALL ppm_alloc(topo%min_physs,ldc,iopt,info)
      CALL ppm_alloc(topo%min_physd,ldc,iopt,info)
      CALL ppm_alloc(topo%max_physs,ldc,iopt,info)
      CALL ppm_alloc(topo%max_physd,ldc,iopt,info)

!      !-------------------------------------------------------------------------
!      !  Deallocate the meshes
!      !-------------------------------------------------------------------------
!      CALL ppm_mesh_alloc_equi(topo%mesh,ldc,iopt,info)
!      NULLIFY(topo%mesh)
!      topo%max_meshid = 0

      !-------------------------------------------------------------------------
      !  Mark this topology as not defined
      !-------------------------------------------------------------------------
      topo%isdefined = .FALSE.

      !-------------------------------------------------------------------------
      !  Allow this topo structure to be reused next time
      !-------------------------------------------------------------------------
      ppm_next_avail_topo = topo%ID

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_dealloc',t0,info)
      RETURN
      END SUBROUTINE ppm_topo_dealloc
