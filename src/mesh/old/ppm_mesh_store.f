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

      SUBROUTINE ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info,Offset,ghostsize)

      !!! This routine stores all relevant information about
      !!! a generated mesh on a certain topology.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_error
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_alloc
      USE ppm_module_check_id
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
      !!! Topology ID for which mesh has been created
      INTEGER,                            INTENT(INOUT) :: meshid
      !!! Mesh ID. If <= 0 on input, the
      !!! routine will create a new mesh and return the ID here.
      INTEGER,            DIMENSION(:,:), INTENT(IN   ) :: ndata
      !!! Number of mesh points in each direction on each sub
      INTEGER,            DIMENSION(:,:), INTENT(IN   ) :: istart
      !!! Start of sub mesh in global mesh
      !Cautious! istart and ndata in this routine only kept for compatibility
      !with the old routines
      INTEGER,            DIMENSION(:),   INTENT(IN   ) :: Nm
      !!! Global number of mesh points in the whole comput. domain
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      REAL(MK), OPTIONAL, DIMENSION(:),   INTENT(IN   ) :: Offset

      !!! Offset in each dimension
      INTEGER,  OPTIONAL, DIMENSION(:),   INTENT(IN   ) :: ghostsize
      !!! size of the ghost layer, in number of mesh nodes for each dimension

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_), POINTER :: mesh

      REAL(MK)                        :: t0
      REAL(MK),  DIMENSION(1:ppm_dim) :: Offst

      INTEGER, DIMENSION(ppm_dim) :: ighostsize
      INTEGER                     :: iopt,ld,ud,kk,i,j,isub

      CHARACTER(LEN=ppm_char) :: msg
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

      IF (PRESENT(Offset)) THEN
         Offst(1:ppm_dim)=Offset(1:ppm_dim)
      ELSE
         Offst=0.0_MK
      ENDIF

      IF (PRESENT(ghostsize)) THEN
         ighostsize(1:ppm_dim)=ghostsize(1:ppm_dim)
      ELSE
         ighostsize=0
      ENDIF

      IF (meshid .LE. 0) THEN

         !-----------------------------------------------------------------------
         !  Create a new mesh and mesh identifier if the user specified none
         !-----------------------------------------------------------------------
         ALLOCATE(ppm_t_equi_mesh::mesh, STAT=info)
         or_fail_alloc("mesh pointer")

         CALL mesh%create(topoid,Offst,info,Nm=Nm,ghostsize=ighostsize)
         or_fail("Failed to create mesh!")

         !-----------------------------------------------------------------------
         !  Push the created mesh and mesh identifier into the collection
         !-----------------------------------------------------------------------
!          CALL ppm_mesh%push(mesh,info,meshid)
         meshid=ppm_mesh%get_id(mesh)
         mesh => NULL()

         ppm_mesh%vec(meshid)%t%ID=meshid

         WRITE(msg,'(a,a,i0)')"default_mesh_name","_",meshid
         ppm_mesh%vec(meshid)%t%name=TRIM(msg)

      ELSE

         IF (meshid .GT. ppm_mesh%max_id) THEN
            fail('meshid is invalid',ppm_err_argument,exit_point=9999)
         ENDIF

         SELECT CASE (ASSOCIATED(ppm_mesh%vec(meshid)%t))
         CASE (.FALSE.)
            ALLOCATE(ppm_t_equi_mesh::ppm_mesh%vec(meshid)%t,STAT=info)
            or_fail_alloc("ppm_mesh%vec(meshid)%t pointer")
         END SELECT

         mesh => ppm_mesh%vec(meshid)%t

         WRITE(msg,'(a,a,i0)')"default_mesh_name","_",meshid

         CALL mesh%create(topoid,Offst,info,Nm=Nm,ghostsize=ighostsize,name=TRIM(msg))

         !-------------------------------------------------------------------------
         !  Store the user-provided ID of this mesh
         !-------------------------------------------------------------------------
         mesh%ID=meshid

      ENDIF

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
             fail('topoid not valid',exit_point=8888)
          ENDIF
          IF (ppm_topo(topoid)%t%nsubs .LE. 0) THEN
             fail('nsubs must be >0',exit_point=8888)
          ENDIF
          DO i=1,ppm_dim
             IF (Nm(i) .LT. 2) THEN
                fail('Nm must be >1 in all space dimensions',exit_point=8888)
             ENDIF
          ENDDO
      8888 CONTINUE
      END SUBROUTINE check

      END SUBROUTINE ppm_mesh_store

