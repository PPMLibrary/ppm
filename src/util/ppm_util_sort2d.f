      !  Subroutine   :                   ppm_util_sort2d
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
      SUBROUTINE ppm_util_sort2d_s(xp,Np,xmin,xmax,Nm,npbx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_sort2d_d(xp,Np,xmin,xmax,Nm,npbx,info)
#endif
      !!! Re-orders the particle locations such that
      !!! subsequent particles are in the same box (cell)
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmin
      !!! Minimum extent of cell mesh
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmax
      !!! Minimum extent of cell mesh
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm     
      !!! Number of cells in each direction
      INTEGER , DIMENSION(:)  , POINTER       :: npbx
      !!! Number of particles in each box
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of particles
      INTEGER                 , INTENT(INOUT) :: info
      !!! If info = 1 on input the extent of the particles will be checked.
      !!!
      !!! return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), POINTER      :: work => NULL()
      ! timer
      REAL(MK)                               :: t0
      !  counters
      INTEGER                                :: ipart,ibox,i
      ! number of cells
      INTEGER                                :: nbox
      ! dimensions for alloc
      INTEGER, DIMENSION(2)                  :: lda
      INTEGER                                :: iopt
      ! index list of particles in cells (allocated within rank2d)
      INTEGER, DIMENSION(:), POINTER         :: lpdx => NULL()
      ! pointer to first particle in each cell (allocated within rank2d)
      INTEGER, DIMENSION(:), POINTER         :: lhbx => NULL()
      ! dummy array to store the number of ghost layers (0)
      INTEGER, DIMENSION(4)                  :: Ngl
      ! local info level
      INTEGER                                :: info2
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      ! save input info (substart will reset info to 0)
      info2 = info
      CALL substart('ppm_util_sort2d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Total number of cells
      !-------------------------------------------------------------------------
      nbox = Nm(1)*Nm(2)  

      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      lda(1) = 2
      lda(2) = Np
      CALL ppm_alloc(work,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_sort2d',     &
     &        'work array WORK',__LINE__,info)
          GOTO 9999
      ENDIF
      lda(1) = nbox
      CALL ppm_alloc(npbx,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_sort2d',     &
     &        'number of particles per box NPBX',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Call ppm_util_rank2d to get the particle index arrays
      !-------------------------------------------------------------------------
      DO i=1,4
         Ngl(i) = 0
      ENDDO
      CALL ppm_util_rank2d(xp,Np,xmin,xmax,Nm,Ngl,lpdx,lhbx,info2)

      ! check if all particles have been ranked
      IF (info2 .GT. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_part_lost,'ppm_util_sort2d',     &
     &        'Not all particles have been ranked',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Re-arrange the particles in the correct order
      !-------------------------------------------------------------------------
      DO ibox=1,nbox
          ! compute the number of particles per box (NEEDED by the calling
          ! routine ppm_decomp_boxsplit)
          npbx(ibox) = lhbx(ibox+1) - lhbx(ibox)
          DO ipart=lhbx(ibox),(lhbx(ibox+1)-1)
              work(1,ipart)      = xp(1,lpdx(ipart))
              work(2,ipart)      = xp(2,lpdx(ipart))
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Copy back
      !-------------------------------------------------------------------------
      DO ipart=1,Np
         xp(1,ipart) = work(1,ipart)
         xp(2,ipart) = work(2,ipart)
      ENDDO

      !-------------------------------------------------------------------------
      !  Free work memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(lpdx,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort2d',     &
     &        'particle index list LPDX',__LINE__,info) 
      ENDIF
      CALL ppm_alloc(lhbx,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort2d',     &
     &        'box header pointers LHBX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(work,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort2d',     &
     &        'work array WORK',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_sort2d',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Np must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Nm(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Nm(1) must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Nm(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Nm(2) must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (xmax(i) .LE. xmin(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &                'xmin must be < xmax',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_sort2d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_sort2d_d
#endif
