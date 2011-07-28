      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_util_sort2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Re-orders the particle locations such that
      !                 subsequent particles are in the same box (cell).
      !
      !  Input        : Np         (I) number of particles
      !                 xmin(:)    (F) minimum extent of cell mesh
      !                 xmax(:)    (F) maximum extent of cell mesh
      !                 Nm(:)      (I) number of cells in each direct.
      !
      !  Input/output : xp(:,:)    (F) particle co-ordinates. Sorted upon
      !                                output.
      !                 info       (I) return status. If info = 1 on input
      !                                the extent of the particles
      !                                will be checked.
      !
      !  Output       : npbx(nbox) (I) number of particles in a box
      !                 
      !  Remarks      : Two do loops do not vectorize. If particles are
      !                 outside the mesh the code will fail.
      !
      !                 The routine uses no (0) automatic arrays! since
      !                 they silently fail when resources are exhausted.
      !
      !                 Warning! This routine may loose real particles if 
      !                 used together with ghost layers
      !
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_sort2d.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.16  2006/09/04 18:35:00  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.14  2004/10/28 12:38:18  davidch
      !  Fixed numerical bug in cell lists that resulted in real particles being treated as ghosts
      !  and vice versa. The new ranking and cell list routines are supposed to be exact. All
      !  epsilons that were added to the domains in order to prevent the mentioned problems were
      !  removed since they are no longer needed.
      !  Modified Files:
      !      ppm_util_rank2d.f ppm_util_rank3d.f ppm_util_sort2d.f
      !      ppm_util_sort3d.f ppm_find_duplicates.f ppm_neighlist_clist.f
      !      ppm_error.h
      !
      !  Revision 1.13  2004/10/01 16:09:14  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.12  2004/07/26 07:42:34  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/06/10 16:20:05  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.10  2004/03/22 09:54:17  walther
      !  Bug fix: Np can be zero - the building of the tree relies on this.
      !
      !  Revision 1.9  2004/01/23 17:22:13  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.8  2004/01/08 17:49:07  ivos
      !  Added dealloc for lpdx and lhbx at the end.
      !
      !  Revision 1.7  2004/01/08 12:54:04  ivos
      !  Now checks if all particles have been ranked and terminates if not.
      !  This is needed since ppm_util_rank no longer performs this check (it
      !  would be harmful for the cell lists).
      !
      !  Revision 1.6  2004/01/06 13:42:41  ivos
      !  Now uses ppm_util_rank instead of ppm_util_cellsort.
      !
      !  Revision 1.5  2003/12/05 10:10:07  ivos
      !  Bugfix. It now compiles.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_sort2d_s(xp,Np,xmin,xmax,Nm,npbx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_sort2d_d(xp,Np,xmin,xmax,Nm,npbx,info)
#endif

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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: xmin,xmax
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: Nm     
      INTEGER , DIMENSION(:)  , POINTER       :: npbx
      INTEGER                 , INTENT(IN   ) :: Np
      INTEGER                 , INTENT(INOUT) :: info
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
      ! index list of particles in cells (allocated within rank2d)
      INTEGER, DIMENSION(:), POINTER         :: lpdx
      ! pointer to first particle in each cell (allocated within rank2d)
      INTEGER, DIMENSION(:), POINTER         :: lhbx
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
          IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Np must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Nm(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Nm(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Nm(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &            'Nm(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (xmax(i) .LE. xmin(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_util_sort2d',  &
     &                'xmin must be < xmax',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
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
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_util_sort2d_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_util_sort2d_d
#endif
