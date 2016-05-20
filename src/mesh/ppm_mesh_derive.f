      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_mesh_derive
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

      SUBROUTINE ppm_mesh_derive(topoid,template_meshid,new_meshid,act,factor, &
     &                           info)
      !!! This routine derives a new mesh from an existing one by refining
      !!! or coarsening the grid cells by a certain factor. For coarsening
      !!! the number of cells in the original mesh needs to be divisible by
      !!! that factor on every sub in every direction.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_mesh_store
      USE ppm_module_alloc
      USE ppm_module_error
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
      !!! Topology ID for which mesh has been created,
      INTEGER                 , INTENT(IN   ) :: template_meshid
      !!! Mesh id of the template mesh. It will NOT be overwritten.
      INTEGER                 , INTENT(IN   ) :: act
      !!! Action. One of:
      !!!
      !!! * ppm_param_mesh_refine
      !!! * ppm_param_mesh_coarsen
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: factor
      !!! Factor by which to refine/coarsen the mesh in each direction.
      !!! For coarsening the number of mesh points on the old mesh must
      !!! be divisible by this factor in every direction.
      INTEGER                 , INTENT(INOUT) :: new_meshid
      !!! Out mesh ID. If <= 0 on input, the
      !!! routine will create a new mesh and return the ID here.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)            :: t0
      INTEGER                          :: nsubs,i,j,iopt
      INTEGER, DIMENSION(2)            :: ldc
      INTEGER, DIMENSION(ppm_dim)      :: Nm
      INTEGER, DIMENSION(:,:), POINTER :: nno => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: ist => NULL()
      CHARACTER(LEN=ppm_char)          :: mesg
      LOGICAL                          :: valid
      TYPE(ppm_t_equi_mesh), POINTER   :: p_mesh
      TYPE(ppm_t_topo)     , POINTER   :: topo
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_mesh_derive',t0,info)



      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
      topo => ppm_topo(topoid)%t
      p_mesh => topo%mesh(template_meshid)
      !-------------------------------------------------------------------------
      !  Allocate memory for existing mesh data
      !-------------------------------------------------------------------------
      nsubs  = topo%nsubs

      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      CALL ppm_alloc(nno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_mesh_derive',   &
     &        'old number of grid points NNO',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ist,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_mesh_derive',   &
     &        'old starting points IST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Read existing mesh data and double-check topoid
      !-------------------------------------------------------------------------
      nno = p_mesh%nnodes(1:ppm_dim,1:nsubs)
      ist = p_mesh%istart(1:ppm_dim,1:nsubs)
      Nm  = p_mesh%Nm(1:ppm_dim)

      !-------------------------------------------------------------------------
      !  Coarsen mesh
      !-------------------------------------------------------------------------
      IF (act .EQ. ppm_param_mesh_coarsen) THEN
          DO j=1,ppm_dim
              IF (MOD(Nm(j)-1,factor(j)) .NE. 0) THEN
                  info = ppm_error_error
                  WRITE(mesg,'(A,I1,A,I5,A,I6)') 'Global mesh in dimension', &
     &                j,' has ',Nm(j)-1,        &
     &                ' mesh cells and is not divisible by ',factor(j)
                  CALL ppm_error(ppm_err_bad_meshop,'ppm_mesh_derive',  &
     &                mesg,__LINE__,info)
                  GOTO 9999
              ENDIF
              Nm(j) = ((Nm(j)-1)/factor(j))+1
          ENDDO
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF (MOD(nno(j,i)-1,factor(j)) .NE. 0) THEN
                      info = ppm_error_error
                      WRITE(mesg,'(A,I5,A,I1,A,I6,A,I3)') 'Mesh on sub ', &
     &                    i,' in dimension ',j,' has ',nno(j,i)-1,        &
     &                    ' mesh cells and is not divisible by ',factor(j)
                      CALL ppm_error(ppm_err_bad_meshop,'ppm_mesh_derive',  &
     &                    mesg,__LINE__,info)
                      GOTO 9999
                  ENDIF
                  nno(j,i) = ((nno(j,i)-1)/factor(j))+1
                  ist(j,i) = ((ist(j,i)-1)/factor(j))+1
              ENDDO
          ENDDO

      !-------------------------------------------------------------------------
      !  Refine mesh
      !-------------------------------------------------------------------------
      ELSEIF (act .EQ. ppm_param_mesh_refine) THEN
          DO j=1,ppm_dim
              Nm(j) = ((Nm(j)-1)*factor(j))+1
          ENDDO
          DO i=1,nsubs
              DO j=1,ppm_dim
                  nno(j,i) = ((nno(j,i)-1)*factor(j))+1
                  ist(j,i) = ((nno(j,i)-1)*factor(j))+1
              ENDDO
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(topoid,new_meshid,nno,ist,Nm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_mesh_derive',   &
     &        'Storing new mesh failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_derive',   &
     &        'old number of grid points NNO',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ist,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mesh_derive',   &
     &        'old starting points IST',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mesh_derive',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_mesh_derive',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'topoid invalid',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(topoid,template_meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'meshid invalid',__LINE__,info)
              GOTO 8888
          ENDIF
          IF ((act .NE. ppm_param_mesh_refine) .AND.              &
     &        (act .NE. ppm_param_mesh_coarsen)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &            'invalid mesh operation specifiec',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (factor(i) .LE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mesh_derive',  &
     &                'factor must be > 0 in all directions',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_derive
