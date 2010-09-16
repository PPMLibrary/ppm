      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_util_sort3d
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
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
#include "ppm_define.h"
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
      !!! Minimum extent of mesh
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmax
      !!! Maximum extent of mesh
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm
      !!! Number of cells in each direction
      INTEGER , DIMENSION(:)  , POINTER       :: npbx
      !!! Number of particles in each box
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of particles
      INTEGER                 , INTENT(INOUT) :: info
      !!! Return status. If info = 1 on input
      !!! the extent of the particles will be checked.
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), POINTER      :: work
      ! timer
      REAL(MK)                               :: t0
      !  counters
      INTEGER                                :: ipart,ibox,i
      ! number of cells
      INTEGER                                :: nbox
      ! dimensions for alloc
      INTEGER, DIMENSION(2)                  :: lda
      INTEGER                                :: iopt
      ! index list of particles in cells (allocated within rank3d)
      INTEGER, DIMENSION(:), POINTER         :: lpdx
      ! pointer to first particle in each cell (allocated within rank3d)
      INTEGER, DIMENSION(:), POINTER         :: lhbx
      ! dummy array to store the number of ghost layers (0)
      INTEGER, DIMENSION(6)                  :: Ngl
      ! local info level
      INTEGER                                :: info2
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      ! save input info level (substart will reset info to 0)
      info2 = info
      CALL substart('ppm_util_sort3d',t0,info)

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
      nbox = Nm(1)*Nm(2)*Nm(3)   ! total number of boxes

      !-------------------------------------------------------------------------
      !  Allocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      lda(1) = 3
      lda(2) = Np
      CALL ppm_alloc(work,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_sort3d',     &
     &        'work array WORK',__LINE__,info)
          GOTO 9999
      ENDIF
      lda(1) = nbox
      CALL ppm_alloc(npbx,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_sort3d',     &
     &        'number of particles per box NPBX',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Call ppm_util_rank3d to get the particle index arrays
      !-------------------------------------------------------------------------
      DO i=1,6
         Ngl(i) = 0
      ENDDO
      CALL ppm_util_rank3d(xp,Np,xmin,xmax,Nm,Ngl,lpdx,lhbx,info2)

      ! check if all particles have been ranked
      IF (info2 .GT. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_part_lost,'ppm_util_sort3d',     &
     &        'Not all particles have been ranked',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Re-arrange the particles in the correct order
      !-------------------------------------------------------------------------
      DO ibox=1,nbox
          ! compute number of particles per box (NEEDED by the calling
          ! routine ppm_decomp_boxsplit)
          npbx(ibox) = lhbx(ibox+1) - lhbx(ibox)
          DO ipart=lhbx(ibox),(lhbx(ibox+1)-1)
              work(1,ipart)      = xp(1,lpdx(ipart))
              work(2,ipart)      = xp(2,lpdx(ipart))
              work(3,ipart)      = xp(3,lpdx(ipart))
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
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort3d',     &
     &        'particle index list LPDX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(lhbx,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort3d',     &
     &        'box header pointers LHBX',__LINE__,info)
      ENDIF
      CALL ppm_alloc(work,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_sort3d',     &
     &        'work array WORK',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_sort3d',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort3d',  &
     &            'Np must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Nm(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort3d',  &
     &            'Nm(1) must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Nm(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort3d',  &
     &            'Nm(2) must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (Nm(3) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort3d',  &
     &            'Nm(3) must be >0',__LINE__,info)
              GOTO 8888
          ENDIF
          DO i=1,ppm_dim
              IF (xmax(i) .LE. xmin(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_util_sort3d',  &
     &                'xmin must be < xmax',__LINE__,info)
                  GOTO 8888
              ENDIF
          ENDDO
 8888     CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_sort3d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_sort3d_d
#endif
