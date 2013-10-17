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
      &          min_sub,max_sub,nsubs,istart,ndata,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mesh_on_subs_d(Nm,min_phys,max_phys, &
      &          min_sub,max_sub,nsubs,istart,ndata,info)
#endif
      !!! This routine defines meshes on a given collection of subs. The
      !!! subs must have been defined such that their extent is an integer
      !!! multiple of the mesh spacing in all directions. If this is not
      !!! the case this routine will fail.
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
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! The number of mesh points (not cells) in each direction of the
      !!! global comput. domain. (including those ON the boundaries)
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys
      !!! The minimum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: max_phys
      !!! The maximum coordinate of the physical/computational domain
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub
      !!! Min. extent of the subdomains
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_sub
      !!! Max. extent of the subdomains
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! Total number of subdomains
      INTEGER , DIMENSION(:,:), POINTER       :: istart
      !!! Start indices (i,j,k) (first index) of mesh in sub isub
      !!! (second index) in global mesh.
      INTEGER , DIMENSION(:,:), POINTER       :: ndata
      !!! Number of grid points in x,y[,z] (first index) of mesh on sub
      !!! isub (second index).
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                     :: t0,lmyeps
      REAL(MK), DIMENSION(ppm_dim) :: len_phys,dx,rat

      INTEGER, DIMENSION(ppm_dim) :: ldu,Nc
      INTEGER                     :: iopt,i,j,k

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller = 'ppm_mesh_on_subs'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
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
      !  Check that the subs align with the mesh points and determine
      !  number of mesh points. 2D and 3D case have separate loops for
      !  vectorization.
      !-------------------------------------------------------------------------
      SELECT CASE (ppm_dim)
      CASE (2)
          DO i=1,nsubs
             len_phys(1) = max_sub(1,i)-min_sub(1,i)
             len_phys(2) = max_sub(2,i)-min_sub(2,i)
             rat(1) = len_phys(1)/dx(1)
             rat(2) = len_phys(2)/dx(2)
             Nc(1)  = NINT(rat(1))
             Nc(2)  = NINT(rat(2))
             IF (ABS(rat(1)-REAL(Nc(1),MK)) .GT. lmyeps*rat(1)) THEN
                WRITE(mesg,'(2(A,F12.6))') 'in dimension 1: sub_length=',  &
                & len_phys(1),' mesh_spacing=',dx(1)
                info = ppm_error_error
                CALL ppm_error(ppm_err_subs_incomp,caller,     &
                &    mesg,__LINE__,info)
                GOTO 9999
             ENDIF
             IF (ABS(rat(2)-REAL(Nc(2),MK)) .GT. lmyeps*rat(2)) THEN
                  WRITE(mesg,'(2(A,F12.6))') 'in dimension 2: sub_length=',  &
                  & len_phys(2),' mesh_spacing=',dx(2)
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_subs_incomp,caller,     &
                  &    mesg,__LINE__,info)
                  GOTO 9999
              ENDIF
              ndata(1,i) = Nc(1) + 1
              ndata(2,i) = Nc(2) + 1
              !-----------------------------------------------------------------
              !  Determine the start indices in the global mesh
              !-----------------------------------------------------------------
              istart(1,i) = NINT((min_sub(1,i)-min_phys(1))/dx(1)) + 1
              istart(2,i) = NINT((min_sub(2,i)-min_phys(2))/dx(2)) + 1
          ENDDO

      CASE DEFAULT
         DO i=1,nsubs
            len_phys(1:3) = max_sub(1:3,i)-min_sub(1:3,i)
            rat(1:3)      = len_phys(1:3)/dx(1:3)
            Nc(1:3)       = NINT(rat(1:3))
            IF (ABS(rat(1)-REAL(Nc(1),MK)) .GT. lmyeps*rat(1)) THEN
               WRITE(mesg,'(2(A,F12.6))') 'in dimension 1: sub_length=',  &
               & len_phys(1),' mesh_spacing=',dx(1)
               info = ppm_error_error
               CALL ppm_error(ppm_err_subs_incomp,caller,     &
               &    mesg,__LINE__,info)
               GOTO 9999
            ENDIF
            IF (ABS(rat(2)-REAL(Nc(2),MK)) .GT. lmyeps*rat(2)) THEN
               WRITE(mesg,'(2(A,F12.6))') 'in dimension 2: sub_length=',  &
               & len_phys(2),' mesh_spacing=',dx(2)
               info = ppm_error_error
               CALL ppm_error(ppm_err_subs_incomp,caller,     &
               &    mesg,__LINE__,info)
               GOTO 9999
            ENDIF
            IF (ABS(rat(3)-REAL(Nc(3),MK)) .GT. lmyeps*rat(3)) THEN
               WRITE(mesg,'(2(A,F12.6))') 'in dimension 3: sub_length=',  &
               & len_phys(3),' mesh_spacing=',dx(3)
               info = ppm_error_error
               CALL ppm_error(ppm_err_subs_incomp,caller,     &
               &    mesg,__LINE__,info)
               GOTO 9999
            ENDIF
            ndata(1:3,i) = Nc(1:3) + 1
            !-----------------------------------------------------------------
            !  Determine the start indices in the global mesh
            !-----------------------------------------------------------------
            istart(1:3,i) = NINT((min_sub(1:3,i)-min_phys(1:3))/dx(1:3)) + 1
         ENDDO

      END SELECT

      !-------------------------------------------------------------------------
      !  Some diagnostics
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          WRITE(mesg,'(A,I5)') 'number of meshes on subs created: ',nsubs
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
      IF (ppm_debug .GT. 1) THEN
          IF (ppm_dim .LT. 3) THEN
              DO i=1,nsubs
                  WRITE(mesg,'(I4,A,2I8)') i,' istart: ',istart(1:2,i)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    mesg,info)
                  WRITE(mesg,'(I4,A,2I8)') i,' ndata : ',ndata(1:2,i)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    mesg,info)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    '------------------------------------',info)
              ENDDO
          ELSE
              DO i=1,nsubs
                  WRITE(mesg,'(I4,A,3I8)') i,' istart: ',istart(1:3,i)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    mesg,info)
                  WRITE(mesg,'(I4,A,3I8)') i,' ndata : ',ndata(1:3,i)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    mesg,info)
                  CALL ppm_write(ppm_rank,caller,  &
                  &    '------------------------------------',info)
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
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'nsubs must be > 0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(min_sub,2) .LT. nsubs .OR. SIZE(max_sub,2) .LT. nsubs) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'not enough subs specified in min_sub, max_sub',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (SIZE(min_sub,1).LT.ppm_dim.OR.SIZE(max_sub,1).LT.ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,caller,  &
     &            'leading dimension in min_sub, max_sub wrong',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF (max_sub(j,i) .LT. min_sub(j,i)) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_argument,caller,  &
     &                    'max_sub must always be >= min_sub',__LINE__,info)
                      GOTO 8888
                  ENDIF
              ENDDO
          ENDDO
          DO i=1,ppm_dim
              IF (Nm(i) .LT. 2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,caller,  &
     &                'Nm must be >1 in all space dimensions',__LINE__,info)
                  GOTO 8888
              ENDIF
              IF (max_phys(i) .LE. min_phys(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,caller,  &
     &                'max_phys must be > min_phys',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mesh_on_subs_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mesh_on_subs_d
#endif
