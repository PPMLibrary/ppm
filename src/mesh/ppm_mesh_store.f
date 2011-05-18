      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_store
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

      SUBROUTINE ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info)
      !!! This routine stores all relevant information about
      !!! a generated mesh on a certain topology.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_util_invert_list
      USE ppm_module_error
      USE ppm_module_mesh_alloc
      USE ppm_module_alloc
      USE ppm_module_mesh_alloc
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
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ndata
      !!! Number of mesh points in each direction on each sub
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: istart
      !!! Start of sub mesh in global mesh
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Global number of mesh points in the whole comput. domain
      INTEGER                 , INTENT(INOUT) :: meshid
      !!! Mesh ID. If <= 0 on input, the
      !!! routine will create a new mesh and return the ID here.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3)    :: ldc
      INTEGER                   :: iopt,ld,ud,kk,i,j,isub
      REAL(ppm_kind_double)     :: t0
      LOGICAL                   :: valid
      TYPE(ppm_t_topo), POINTER :: topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_store',t0,info)

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
        meshid = topo%max_meshid + 1
        topo%max_meshid = meshid
      ELSE
        IF (meshid .GT. topo%max_meshid) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_mesh_store', &
     &                     'meshid is invalid',__LINE__,info)
            GOTO 9999
        ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the internal mesh list and Arrays at meshid
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldc(1) = topo%max_meshid
      CALL ppm_mesh_alloc_equi(topo%mesh,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'user meshid list ppm_t_topo%mesh',__LINE__,info)
         GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldc(1) = SIZE(ndata,1)
      ldc(2) = SIZE(ndata,2)
      CALL ppm_alloc(topo%mesh(meshid)%nnodes,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'ppm_t_equi_mesh%nnodes',__LINE__,info)
         GOTO 9999
      ENDIF

      ldc(1) = SIZE(istart,1)
      ldc(2) = SIZE(istart,2)
      CALL ppm_alloc(topo%mesh(meshid)%istart,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'ppm_t_equi_mesh%istart',__LINE__,info)
         GOTO 9999
      ENDIF

      ldc(1) = SIZE(Nm,1)
      CALL ppm_alloc(topo%mesh(meshid)%Nm,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_mesh_store',  &
     &       'ppm_t_equi_mesh%Nm',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the user-provided ID of this mesh
      !-------------------------------------------------------------------------
      topo%mesh(meshid)%ID = meshid

      !-------------------------------------------------------------------------
      !  Store the mesh information
      !-------------------------------------------------------------------------

      IF (ppm_dim .EQ. 3) THEN
          DO isub = 1,topo%nsubs
              topo%mesh(meshid)%nnodes(1,isub) = ndata(1,isub)
              topo%mesh(meshid)%nnodes(2,isub) = ndata(2,isub)
              topo%mesh(meshid)%nnodes(3,isub) = ndata(3,isub)
          ENDDO
          DO isub = 1,topo%nsubs
              topo%mesh(meshid)%istart(1,isub) = istart(1,isub)
              topo%mesh(meshid)%istart(2,isub) = istart(2,isub)
              topo%mesh(meshid)%istart(3,isub) = istart(3,isub)
          ENDDO
          topo%mesh(meshid)%Nm(1) = Nm(1)
          topo%mesh(meshid)%Nm(2) = Nm(2)
          topo%mesh(meshid)%Nm(3) = Nm(3)

      ELSE
          DO isub = 1,topo%nsubs
              topo%mesh(meshid)%nnodes(1,isub) = ndata(1,isub)
              topo%mesh(meshid)%nnodes(2,isub) = ndata(2,isub)
          ENDDO
          DO isub = 1,topo%nsubs
              topo%mesh(meshid)%istart(1,isub) = istart(1,isub)
              topo%mesh(meshid)%istart(2,isub) = istart(2,isub)
          ENDDO
          topo%mesh(meshid)%Nm(1) = Nm(1)
          topo%mesh(meshid)%Nm(2) = Nm(2)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_store',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &            'topoid not valid',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (ppm_topo(topoid)%t%nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &            'nsubs must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (Nm(i) .LT. 2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mesh_store',  &
     &                'Nm must be >1 in all space dimensions',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_store
