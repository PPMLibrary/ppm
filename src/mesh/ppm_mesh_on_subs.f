      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_mesh_on_subs
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mesh_on_subs_s(Nm,min_phys,max_phys, &
      &          min_sub,max_sub,nsubs,istart,ndata,info,Offset)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mesh_on_subs_d(Nm,min_phys,max_phys, &
      &          min_sub,max_sub,nsubs,istart,ndata,info,Offset)
#endif
      !!! This routine defines meshes on a given collection of subs.
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(:  ),           INTENT(IN   ) :: Nm
      !!! The number of mesh points (not cells) in each direction of the
      !!! global comput. domain. (including those ON the boundaries)
      REAL(MK), DIMENSION(:  ),           INTENT(IN   ) :: min_phys
      !!! The minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:  ),           INTENT(IN   ) :: max_phys
      !!! The maximum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:,:),           INTENT(IN   ) :: min_sub
      !!! Min. extent of the subdomains
      REAL(MK), DIMENSION(:,:),           INTENT(IN   ) :: max_sub
      !!! Max. extent of the subdomains
      INTEGER,                            INTENT(IN   ) :: nsubs
      !!! Total number of subdomains
      INTEGER , DIMENSION(:,:),           POINTER       :: istart
      !!! Start indices (i,j,k) (first index) of mesh in sub isub
      !!! (second index) in global mesh.
      INTEGER , DIMENSION(:,:),           POINTER       :: ndata
      !!! Number of grid points in x,y[,z] (first index) of mesh on sub
      !!! isub (second index).
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
#if   __KIND == __SINGLE_PRECISION
      REAL(ppm_kind_double), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: Offset
      !!! Offset in each dimension
#elif __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: Offset
      !!! Offset in each dimension
#endif
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)        :: t0
      REAL(MK), DIMENSION(ppm_dim) :: len_phys,dx
      REAL(MK), DIMENSION(ppm_dim) :: Offst

      INTEGER, DIMENSION(ppm_dim) :: iend
      !!! Upper-right coordinates on the sub mesh
      INTEGER, DIMENSION(ppm_dim) :: ldu,Nc
      INTEGER                     :: iopt,i,j,k

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_mesh_on_subs'
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
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      IF (PRESENT(Offset)) THEN
#if   __KIND == __SINGLE_PRECISION
         Offst(1:ppm_dim) = REAL(Offset(1:ppm_dim),MK)
#elif __KIND == __DOUBLE_PRECISION
         Offst(1:ppm_dim) = Offset(1:ppm_dim)
#endif
      ELSE
         Offst = 0.0_MK
      ENDIF

      !-------------------------------------------------------------------------
      !  Mesh spacing
      !-------------------------------------------------------------------------
      Nc(1:ppm_dim)       = Nm(1:ppm_dim)-1
      len_phys(1:ppm_dim) = max_phys(1:ppm_dim) - min_phys(1:ppm_dim)
      dx(1:ppm_dim)       = len_phys(1:ppm_dim)/REAL(Nc(1:ppm_dim),MK)

      !check for round-off problems and fix them if necessary
      DO k=1,ppm_dim
         DO WHILE (min_phys(k)+Nc(k)*dx(k).LT.max_phys(k))
            dx(k)=dx(k)+EPSILON(dx(k))
         ENDDO
      ENDDO
      check_true(<#ALL(min_phys(1:ppm_dim)+Nc(1:ppm_dim)*dx(1:ppm_dim).GE.max_phys(1:ppm_dim))#>,"round-off problem in mesh creation")

      !-------------------------------------------------------------------------
      !  Allocate memory for the meshes
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu(1) = ppm_dim
      ldu(2) = nsubs
      CALL ppm_alloc(istart,ldu,iopt,info)
      or_fail_alloc('sub mesh start indices ISTART')

      CALL ppm_alloc(ndata,ldu,iopt,info)
      or_fail_alloc('sub mesh sizes NDATA')

      !-------------------------------------------------------------------------
      !  Determine number of mesh points.
      !-------------------------------------------------------------------------
      DO i=1,nsubs
         istart(1:ppm_dim,i)=1+CEILING((min_sub(1:ppm_dim,i)-Offst(1:ppm_dim))/dx(1:ppm_dim))
         iend(1:ppm_dim)    =1+FLOOR((  max_sub(1:ppm_dim,i)-Offst(1:ppm_dim))/(dx(1:ppm_dim)-EPSILON(dx(1:ppm_dim))))
         ndata(1:ppm_dim,i) =1+iend(1:ppm_dim)-istart(1:ppm_dim,i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Some diagnostics
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         stdout_f('(A,I5)',"number of meshes on subs created: ",nsubs)
      ENDIF
      IF (ppm_debug .GT. 1) THEN
         IF (ppm_dim .LT. 3) THEN
            DO i=1,nsubs
               stdout_f('(I4,A,2I8)',i," istart: ",'istart(1:2,i)')
               stdout_f('(I4,A,2I8)',i," ndata : ",'ndata(1:2,i)')
               stdout_f('(A)',"------------------------------------")
            ENDDO
         ELSE
            DO i=1,nsubs
               stdout_f('(I4,A,3I8)',i," istart: ",'istart(1:3,i)')
               stdout_f('(I4,A,3I8)',i," ndata : ",'ndata(1:3,i)')
               stdout_f('(A)',"------------------------------------")
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (nsubs .LE. 0) THEN
             fail('nsubs must be > 0',exit_point=8888)
          ENDIF
          IF (SIZE(min_sub,2) .LT. nsubs .OR. SIZE(max_sub,2) .LT. nsubs) THEN
             fail('not enough subs specified in min_sub, max_sub',exit_point=8888)
          ENDIF
          IF (SIZE(min_sub,1).LT.ppm_dim.OR.SIZE(max_sub,1).LT.ppm_dim) THEN
             fail('leading dimension in min_sub, max_sub wrong',exit_point=8888)
          ENDIF
          DO i=1,nsubs
             DO j=1,ppm_dim
                IF (max_sub(j,i) .LT. min_sub(j,i)) THEN
                   fail('max_sub must always be >= min_sub',exit_point=8888)
                ENDIF
             ENDDO
          ENDDO
          DO i=1,ppm_dim
             IF (Nm(i) .LT. 2) THEN
                fail('Nm must be >1 in all space dimensions',exit_point=8888)
             ENDIF
             IF (max_phys(i) .LE. min_phys(i)) THEN
                fail('max_phys must be > min_phys',exit_point=8888)
             ENDIF
          ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mesh_on_subs_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mesh_on_subs_d
#endif
