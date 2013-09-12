      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_store
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

      SUBROUTINE ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info)
      !!! This routine stores all relevant information about
      !!! a generated mesh on a certain topology.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      !USE ppm_module_data_mesh
      USE ppm_module_topo_typedef
      USE ppm_module_mesh_typedef
      !USE ppm_module_util_invert_list
      USE ppm_module_error
      !USE ppm_module_mesh_alloc
      !USE ppm_module_alloc
      !USE ppm_module_mesh_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_check_id

      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! Topology ID for which mesh has been created
      INTEGER                 , INTENT(INOUT) :: meshid
      !!! Mesh ID. If <= 0 on input, the
      !!! routine will create a new mesh and return the ID here.
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ndata
      !!! Number of mesh points in each direction on each sub
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: istart
      !!! Start of sub mesh in global mesh
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Global number of mesh points in the whole comput. domain
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo),        POINTER :: topo => NULL()
      CLASS(ppm_t_equi_mesh_), POINTER :: mesh => NULL()

      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(3) :: ldc
      INTEGER               :: iopt,ld,ud,kk,i,j,isub

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_mesh_store'

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

      IF (meshid .LE. 0) THEN
         !-----------------------------------------------------------------------
         !  Create a mesh identifier if the user specified none
         !-----------------------------------------------------------------------
         meshid          = ppm_mesh%max_id + 1
         ppm_mesh%max_id = meshid
         ppm_mesh%nb     = ppm_mesh%nb + 1

         IF (.NOT.ppm_mesh%exists(meshid)) THEN

            CALL ppm_mesh%grow_size(info)
            or_fail("ppm_mesh%grow_size()")

         ENDIF

      ELSE

         IF (meshid .GT. ppm_mesh%max_id) THEN
            fail('meshid is invalid',ppm_err_argument,exit_point=9999)
         ENDIF

      ENDIF

      SELECT CASE(ASSOCIATED(ppm_mesh%vec(meshid)%t))
      CASE (.FALSE.)
         ALLOCATE(ppm_t_equi_mesh::ppm_mesh%vec(meshid)%t)
      END SELECT

      mesh => ppm_mesh%vec(meshid)%t

      iopt   = ppm_param_alloc_fit

      ldc(1) = SIZE(ndata,1)
      ldc(2) = SIZE(ndata,2)
      CALL ppm_alloc(mesh%nnodes,ldc,iopt,info)
      or_fail_alloc('mesh%nnodes')

      ldc(1) = SIZE(istart,1)
      ldc(2) = SIZE(istart,2)
      CALL ppm_alloc(mesh%istart,ldc,iopt,info)
      or_fail_alloc('mesh%istart')

      CALL ppm_alloc(mesh%iend,ldc,iopt,info)
      or_fail_alloc('mesh%iend')

      ldc(1) = SIZE(Nm,1)
      CALL ppm_alloc(mesh%Nm,ldc,iopt,info)
      or_fail_alloc('mesh%Nm')

      !-------------------------------------------------------------------------
      !  Store the user-provided ID of this mesh
      !-------------------------------------------------------------------------
      mesh%ID=meshid
      mesh%topoid=topoid

      !-------------------------------------------------------------------------
      !  Store the mesh information
      !-------------------------------------------------------------------------
      DO isub = 1,topo%nsubs
         mesh%nnodes(1:ppm_dim,isub) =  ndata(1:ppm_dim,isub)
         mesh%istart(1:ppm_dim,isub) = istart(1:ppm_dim,isub)
      ENDDO
      mesh%Nm(1:ppm_dim) = Nm(1:ppm_dim)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
              &   'topoid not valid',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (ppm_topo(topoid)%t%nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
              &   'nsubs must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (Nm(i) .LT. 2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
                  &   'Nm must be >1 in all space dimensions',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_store
