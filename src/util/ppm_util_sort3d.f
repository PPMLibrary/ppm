      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_util_sort3d
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
      SUBROUTINE ppm_util_sort3d_s(xp,Np,xmin,xmax,Nm,npbx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_sort3d_d(xp,Np,xmin,xmax,Nm,npbx,info)
#endif
      !!! Re-orders the particle locations such that
      !!! subsequent particles are in the same box (cell).
      !!!
      !!! NOTE: Two do loops do not vectorize. If particles are
      !!! outside the mesh the code will fail.
      !!!
      !!! NOTE: The routine uses no (0) automatic arrays! since
      !!! they silently fail when resources are exhausted.
      !!!
      !!! [WARNING]
      !!! This routine may loose real particles if used together with
      !!! ghost layers
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_util_rank
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      !!! Particle coordinates, sorted upon output
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of particles
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmin
      !!! Minimum extent of mesh
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmax
      !!! Maximum extent of mesh
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm
      !!! Number of cells in each direction
      INTEGER , DIMENSION(:)  , INTENT(INOUT) :: npbx
      !!! Number of particles in each box
      INTEGER                 , INTENT(INOUT) :: info
      !!! Return status. If info = 1 on input
      !!! the extent of the particles will be checked.
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                 :: t0
      REAL(MK), DIMENSION(:,:), ALLOCATABLE :: work

      !  counters
      INTEGER                        :: ipart,ibox,i
      ! number of cells
      INTEGER                        :: nbox
      ! dimensions for alloc
      INTEGER, DIMENSION(2)          :: lda
      INTEGER                        :: iopt
      ! index list of particles in cells (allocated within rank3d)
      INTEGER, DIMENSION(:), POINTER :: lpdx
      ! pointer to first particle in each cell (allocated within rank3d)
      INTEGER, DIMENSION(:), POINTER :: lhbx
      ! dummy array to store the number of ghost layers (0)
      INTEGER, DIMENSION(6)          :: Ngl
      ! local info level
      INTEGER                        :: info2

      CHARACTER(LEN=ppm_char) :: caller = 'ppm_util_sort'
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      ! save input info level (substart will reset info to 0)
      info2 = info
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Call ppm_util_rank3d to get the particle index arrays
      !-------------------------------------------------------------------------
      NULLIFY(lpdx,lhbx)

      Ngl=0

      CALL ppm_util_rank3d(xp,Np,xmin,xmax,Nm,Ngl,lpdx,lhbx,info2)

      ! check if all particles have been ranked
      IF (info2 .GT. 0) THEN
         fail('Not all particles have been ranked',ppm_error=ppm_error_error)
      ENDIF

      !-------------------------------------------------------------------------
      !  Total number of cells
      !-------------------------------------------------------------------------
      nbox = Nm(1)*Nm(2)*Nm(3)   ! total number of boxes

      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
      ALLOCATE(work(ppm_dim,Np), STAT=info)
      or_fail_alloc('work array WORK')

      !-------------------------------------------------------------------------
      !  Re-arrange the particles in the correct order
      !-------------------------------------------------------------------------
      DO ibox=1,nbox
         ! compute number of particles per box (NEEDED by the calling
         ! routine ppm_decomp_boxsplit)
         npbx(ibox) = lhbx(ibox+1) - lhbx(ibox)
         DO ipart=lhbx(ibox),(lhbx(ibox+1)-1)
            work(1,ipart) = xp(1,lpdx(ipart))
            work(2,ipart) = xp(2,lpdx(ipart))
            work(3,ipart) = xp(3,lpdx(ipart))
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy back
      !-------------------------------------------------------------------------
      DO ipart=1,Np
         xp(1,ipart) = work(1,ipart)
         xp(2,ipart) = work(2,ipart)
         xp(3,ipart) = work(3,ipart)
      ENDDO

      !-------------------------------------------------------------------------
      !  Free work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(lpdx,lda,iopt,info)
      or_fail_dealloc('particle index list LPDX')

      CALL ppm_alloc(lhbx,lda,iopt,info)
      or_fail_dealloc('box header pointers LHBX')

      DEALLOCATE(work, STAT=info)
      or_fail_dealloc('work array WORK')

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Np .LT. 0) THEN
           fail('Np must be >0',exit_point=8888,ppm_error=ppm_error_error)
        ENDIF
        IF (Nm(1) .LE. 0) THEN
           fail('Nm(1) must be >0',exit_point=8888,ppm_error=ppm_error_error)
        ENDIF
        IF (Nm(2) .LE. 0) THEN
           fail('Nm(2) must be >0',exit_point=8888,ppm_error=ppm_error_error)
        ENDIF
        IF (Nm(3) .LE. 0) THEN
           fail('Nm(3) must be >0',exit_point=8888,ppm_error=ppm_error_error)
        ENDIF
        DO i=1,ppm_dim
           IF (xmax(i) .LE. xmin(i)) THEN
              fail('xmin must be < xmax',exit_point=8888,ppm_error=ppm_error_error)
           ENDIF
        ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_sort3d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_sort3d_d
#endif
