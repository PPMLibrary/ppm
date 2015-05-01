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

      SUBROUTINE ppm_mesh_derive(topoid,template_meshid,new_meshid,act,factor,info)
      !!! This routine derives a new mesh from an existing one by refining
      !!! or coarsening the grid cells by a certain factor. For coarsening
      !!! the number of cells in the original mesh needs to be divisible by
      !!! that factor on every sub in every direction.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_topo_typedef
      USE ppm_module_mesh_typedef
      USE ppm_module_mesh_define
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
      TYPE(ppm_t_topo)     , POINTER   :: topo
      CLASS(ppm_t_equi_mesh_), POINTER :: mesh

      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(:,:), POINTER :: nno => NULL()
      INTEGER, DIMENSION(1:ppm_dim)    :: Nm
      INTEGER                          :: nsubs,i,j,iopt,meshid_
      INTEGER, DIMENSION(2)            :: ldc

      CHARACTER(LEN=ppm_char) :: caller="ppm_mesh_derive"

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
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      meshid_=-1
      mesh => ppm_mesh%begin()
      DO WHILE (ASSOCIATED(mesh))
         IF (mesh%ID.EQ.template_meshid) THEN
            IF (mesh%topoid.EQ.topoid) THEN
               meshid_=ppm_mesh%iter_id
               EXIT
            ENDIF
         ENDIF
         mesh => ppm_mesh%next()
      ENDDO

      IF (meshid_.LT.ppm_mesh%min_id.OR.meshid_.GT.ppm_mesh%max_id) THEN
         stdout("Input mesh ID",template_meshid,"Is Invalid!")
         fail(cbuf,ppm_error=ppm_error_fatal)
      ENDIF

      !-------------------------------------------------------------------------
      !  Read existing mesh data
      !-------------------------------------------------------------------------
      Nm(1:ppm_dim) = mesh%Nm(1:ppm_dim)

      !-------------------------------------------------------------------------
      !  Coarsen mesh
      !-------------------------------------------------------------------------
      SELECT CASE (act)
      CASE (ppm_param_mesh_coarsen)
         DO j=1,ppm_dim
            IF (MOD(Nm(j)-1,factor(j)).NE.0) THEN
               WRITE(cbuf,'(A,I1,A,I5,A,I6)') 'Global mesh in dimension', &
               & j,' has ',Nm(j)-1,' mesh cells and is not divisible by ',factor(j)
               fail(cbuf,ppm_err_bad_meshop,ppm_error=ppm_error_fatal)
            ENDIF
            Nm(j)=((Nm(j)-1)/factor(j))+1
         ENDDO

         topo => ppm_topo(topoid)%t

         !-------------------------------------------------------------------------
         !  Allocate memory for existing mesh data
         !-------------------------------------------------------------------------
         nsubs = topo%nsubs

         iopt   = ppm_param_alloc_fit
         ldc(1) = ppm_dim
         ldc(2) = nsubs
         CALL ppm_alloc(nno,ldc,iopt,info)
         or_fail_alloc("old number of grid points NNO",ppm_error=ppm_error_fatal)

         !-------------------------------------------------------------------------
         !  Read existing mesh data and double-check topoid
         !-------------------------------------------------------------------------
         nno(1:ppm_dim,1:nsubs) = mesh%iend(1:ppm_dim,1:nsubs)-mesh%istart(1:ppm_dim,1:nsubs)+1

         DO i=1,nsubs
            DO j=1,ppm_dim
               IF (MOD(nno(j,i)-1,factor(j)).NE.0) THEN
                  WRITE(cbuf,'(A,I5,A,I1,A,I6,A,I3)') 'Mesh on sub ', &
                  & i,' in dimension ',j,' has ',nno(j,i)-1,        &
                  & ' mesh cells and is not divisible by ',factor(j)
                  fail(cbuf,ppm_err_bad_meshop,ppm_error=ppm_error_fatal)
               ENDIF
               nno(j,i) = ((nno(j,i)-1)/factor(j))+1
            ENDDO
         ENDDO

      !-------------------------------------------------------------------------
      !  Refine mesh
      !-------------------------------------------------------------------------
      CASE (ppm_param_mesh_refine)
         DO j=1,ppm_dim
            Nm(j) = ((Nm(j)-1)*factor(j))+1
         ENDDO

      END SELECT

      !-------------------------------------------------------------------------
      !  Create and store a new mesh
      !-------------------------------------------------------------------------
      CALL ppm_mesh_define(topoid,new_meshid,Nm,info,Offset=mesh%Offset,ghostsize=mesh%ghostsize)
      or_fail("ppm_mesh_define Failed!")

      !-------------------------------------------------------------------------
      !  Deallocate work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(nno,ldc,iopt,info)
      or_fail_dealloc("old number of grid points NNO")

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
             fail('Topoid invalid!',ppm_err_argument,exit_point=8888)
          ENDIF
          CALL ppm_check_meshid(topoid,template_meshid,valid,info)
          IF (.NOT. valid) THEN
             fail('meshid invalid!',ppm_err_argument,exit_point=8888)
          ENDIF
          IF ((act .NE. ppm_param_mesh_refine) .AND. &
          &   (act .NE. ppm_param_mesh_coarsen)) THEN
              fail('Invalid mesh operation specifiec!',ppm_err_argument,exit_point=8888)
          ENDIF
          DO i=1,ppm_dim
             IF (factor(i) .LE. 0) THEN
                fail('factor must be > 0 in all directions!',ppm_err_argument,exit_point=8888)
             ENDIF
          ENDDO
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_mesh_derive
