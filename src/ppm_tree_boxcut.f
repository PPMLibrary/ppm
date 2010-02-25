      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_boxcut
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine cuts a box into 2, 4 or 8 pieces using
      !                 the given cut directions and cut positions. 
      !
      !  Input        : xp(:,:)      (F) The positions of the points.
      !                 cutbox       (I) ID of the box to be cut.
      !                 min_box(:)   (F) the minimum coordinates of the 
      !                                  box to be cut.
      !                 max_box(:)   (F) the maximum coordinates of the 
      !                                  box to be cut.
      !                 ncut         (I) number of cuts to apply
      !                 cutdir(:)    (I) directions of the cuts.
      !                                  cutdir(1)=2 means that the first
      !                                  cut is perpendicular to the second
      !                                  coordinate axis (y axis).
      !                 cutpos(:)    (F) positions of the cuts along the
      !                                  given axes.
      !
      !  Input/output :  
      !
      !  Output       : mincut(:,:)  (F) minimum coordinates of the cut
      !                                  boxes. 1st: x,y[,z], 2nd: ibox.
      !                 maxcut(:,:)  (F) maximum coordinates of the cut
      !                                  boxes. 1st: x,y[,z], 2nd: ibox.
      !                 lhbx(:)      (I) pointer to the first point
      !                                  (in lpdx) in each sub-box. This is 
      !                                  only allocated and returned if there 
      !                                  are particles at all.
      !                 lpdx(:)      (I) index of particles in sub-boxes. 
      !                                  This is only allocated and returned 
      !                                  if there are particles at all.
      !                 info         (I) return status
      !
      !  Remarks      : The routines only works for ncut=1, 2, or 3.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_boxcut.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.12  2006/10/24 09:48:18  pchatela
      !  Rollback to 1.08 after confusion and merge with wrong version of file
      !
      !  Revision 1.11  2006/09/04 18:34:57  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.9  2006/08/09 16:12:02  ivos
      !  Changed to code that determines the box index for each particle
      !  to a more efficient version that only considers the cut directions.
      !
      !  Revision 1.8  2005/08/31 12:43:45  ivos
      !  Shark optimizations.
      !
      !  Revision 1.7  2005/08/31 11:24:31  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.6  2005/08/30 13:17:26  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.5  2005/02/01 13:22:58  ivos
      !  Moved declarations of lhbx_cut and lpdx_cut to module_data_tree.
      !
      !  Revision 1.4  2004/12/03 17:43:18  ivos
      !  Removed debug output.
      !
      !  Revision 1.3  2004/12/03 17:19:19  ivos
      !  Now determines the new particle ranking lists and returns them.
      !
      !  Revision 1.2  2004/09/22 17:26:30  ivos
      !  bugfix: fixed array boundary overflow in cutpos.
      !
      !  Revision 1.1  2004/09/22 10:32:04  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_boxcut_s(xp,cutbox,min_box,max_box,ncut,cutdir,  &
     &    cutpos,mincut,maxcut,lhbx,lpdx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_boxcut_d(xp,cutbox,min_box,max_box,ncut,cutdir,  &
     &    cutpos,mincut,maxcut,lhbx,lpdx,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_box,max_box,cutpos
      INTEGER                 , INTENT(IN   ) :: ncut,cutbox
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: cutdir
      REAL(MK), DIMENSION(:,:), POINTER       :: mincut,maxcut
      INTEGER , DIMENSION(:  ), POINTER       :: lhbx,lpdx
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER                                 :: nnbox,cd,i,j,iopt,Np,ip
      INTEGER , DIMENSION(2)                  :: ldc
      INTEGER , DIMENSION(:), POINTER         :: pbox
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_boxcut',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (cutbox .LT. 1) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &          'cutbox must be >0 !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF ((ppm_dim .LT. 3) .AND. (ncut .GT. 2)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &          'Cannot cut more than 2 times in 2D!',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            IF (min_box(i) .GT. max_box(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &             'min_box must be <= max_box !',__LINE__,info)
               GOTO 9999
            ENDIF 
         ENDDO
         DO i=1,ncut-1
            DO j=i+1,ncut
               IF (cutdir(i) .EQ. cutdir(j)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &                'cut directions must be distinct !',__LINE__,info)
                  GOTO 9999
               ENDIF 
            ENDDO
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  Exit if there is nothing to be cut
      !-------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_boxcut',    &
     &            'Nothing to cut. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Number of new boxes to be created
      !-------------------------------------------------------------------------
      nnbox = 2**ncut
      iopt  = ppm_param_alloc_grow
      IF (have_particles) THEN
          Np = - tree_lhbx(1,cutbox) + tree_lhbx(2,cutbox) + 1
          ldc(1) = Np
          CALL ppm_alloc(lpdx,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_boxcut',          &
     &            'particle list pointers LPDX',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Split box into 2 pieces
      !-------------------------------------------------------------------------
      IF (ncut .EQ. 1) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd  = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = max_box(cd)
          
      !-------------------------------------------------------------------------
      !  Split box into 4 pieces
      !-------------------------------------------------------------------------
      ELSEIF (ncut .EQ. 2) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = cutpos(1)
          mincut(cd,4) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = cutpos(1)
          maxcut(cd,3) = max_box(cd)
          maxcut(cd,4) = max_box(cd)

          cd = cutdir(2)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(2)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = cutpos(2)

          maxcut(cd,1) = cutpos(2)
          maxcut(cd,2) = max_box(cd)
          maxcut(cd,3) = cutpos(2)
          maxcut(cd,4) = max_box(cd)
      
      !-------------------------------------------------------------------------
      !  Split box into 8 pieces
      !-------------------------------------------------------------------------
      ELSEIF (ncut .EQ. 3) THEN
          IF (ppm_dim .GT. 2) THEN
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  mincut(3,i) = min_box(3)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
                  maxcut(3,i) = max_box(3)
              ENDDO
          ELSE
              DO i=1,nnbox
                  mincut(1,i) = min_box(1)
                  mincut(2,i) = min_box(2)
                  maxcut(1,i) = max_box(1)
                  maxcut(2,i) = max_box(2)
              ENDDO
          ENDIF

          cd = cutdir(1)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = cutpos(1)
          mincut(cd,4) = cutpos(1)

          maxcut(cd,1) = cutpos(1)
          maxcut(cd,2) = cutpos(1)
          maxcut(cd,3) = max_box(cd)
          maxcut(cd,4) = max_box(cd)

          mincut(cd,5) = min_box(cd)
          mincut(cd,6) = min_box(cd)
          mincut(cd,7) = cutpos(1)
          mincut(cd,8) = cutpos(1)

          maxcut(cd,5) = cutpos(1)
          maxcut(cd,6) = cutpos(1)
          maxcut(cd,7) = max_box(cd)
          maxcut(cd,8) = max_box(cd)

          cd = cutdir(2)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = cutpos(2)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = cutpos(2)

          maxcut(cd,1) = cutpos(2)
          maxcut(cd,2) = max_box(cd)
          maxcut(cd,3) = cutpos(2)
          maxcut(cd,4) = max_box(cd)

          mincut(cd,5) = min_box(cd)
          mincut(cd,6) = cutpos(2)
          mincut(cd,7) = min_box(cd)
          mincut(cd,8) = cutpos(2)

          maxcut(cd,5) = cutpos(2)
          maxcut(cd,6) = max_box(cd)
          maxcut(cd,7) = cutpos(2)
          maxcut(cd,8) = max_box(cd)

          cd = cutdir(3)
          mincut(cd,1) = min_box(cd)
          mincut(cd,2) = min_box(cd)
          mincut(cd,3) = min_box(cd)
          mincut(cd,4) = min_box(cd)

          maxcut(cd,1) = cutpos(3)
          maxcut(cd,2) = cutpos(3)
          maxcut(cd,3) = cutpos(3)
          maxcut(cd,4) = cutpos(3)

          mincut(cd,5) = cutpos(3)
          mincut(cd,6) = cutpos(3)
          mincut(cd,7) = cutpos(3)
          mincut(cd,8) = cutpos(3)

          maxcut(cd,5) = max_box(cd)
          maxcut(cd,6) = max_box(cd)
          maxcut(cd,7) = max_box(cd)
          maxcut(cd,8) = max_box(cd)
      
      !-------------------------------------------------------------------------
      !  Unknown cutting mode
      !-------------------------------------------------------------------------
      ELSE
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_tree_boxcut',     &
     &        'ncut must be < 4.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Rank particles into new boxes. We cannot use ppm_util_rank
      !  here, since the child boxes need not be all of the same size!
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          !---------------------------------------------------------------------
          !  Allocate particle box ID list
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldc(1) = Np
          CALL ppm_alloc(pbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_boxcut',   &
     &            'list of particle box IDs PBOX',__LINE__,info)
              GOTO 9999
          ENDIF
          pbox = -1
          npbx = 0

          !---------------------------------------------------------------------
          !  Determine sub-box for each particle. TODO: use the snipplet
          !  here... but avoid alloc!
          !---------------------------------------------------------------------
          IF (ppm_dim .EQ. 3) THEN
              DO i=1,nnbox
                  DO j=1,Np
                      ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
                      IF (.NOT.(xp(1,ip) .LT. mincut(1,i)) .AND.    &
     &                    (xp(1,ip) .LT. maxcut(1,i)) .AND.    &
     &                    .NOT.(xp(2,ip) .LT. mincut(2,i)) .AND.    &
     &                    (xp(2,ip) .LT. maxcut(2,i)) .AND.    &
     &                    .NOT.(xp(3,ip) .LT. mincut(3,i)) .AND.    &
     &                    (xp(3,ip) .LT. maxcut(3,i))) THEN
                          pbox(j) = i
                      ENDIF
                  ENDDO
              ENDDO
          ELSE
              DO i=1,nnbox
                  DO j=1,Np
                      ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
                      IF (.NOT.(xp(1,ip) .LT. mincut(1,i)) .AND.    &
     &                    (xp(1,ip) .LT. maxcut(1,i)) .AND.    &
     &                    .NOT.(xp(2,ip) .LT. mincut(2,i)) .AND.    &
     &                    (xp(2,ip) .LT. maxcut(2,i))) THEN
                          pbox(j) = i
                      ENDIF
                  ENDDO
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Determine number of particles in each sub-box
          !---------------------------------------------------------------------
          DO j=1,Np
              i = pbox(j)
              IF (i .GT. 0) npbx(i) = npbx(i) + 1
          ENDDO

          !---------------------------------------------------------------------
          !  Build the box particle pointers
          !---------------------------------------------------------------------
          cbox(1) = 1
          lhbx(1) = 1
          DO i=2,nnbox
              cbox(i) = cbox(i-1) + npbx(i-1)
              lhbx(i) = cbox(i)
          ENDDO
          lhbx(nnbox+1) = lhbx(nnbox) + npbx(nnbox)

          !---------------------------------------------------------------------
          !  Rank the particles in the sub-boxes
          !---------------------------------------------------------------------
          DO j=1,Np
              ip = tree_lpdx(tree_lhbx(1,cutbox)+j-1)
              i = pbox(j)
              IF (i .GT. 0) THEN
                  lpdx(cbox(i)) = ip
                  cbox(i)       = cbox(i) + 1
              ENDIF
          ENDDO

          !---------------------------------------------------------------------
          !  Deallocate memory
          !---------------------------------------------------------------------
          iopt = ppm_param_dealloc
          CALL ppm_alloc(pbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree_boxcut',     &
     &            'list of particle box IDs PBOX',__LINE__,info)
          ENDIF
      ENDIF           ! have_particles

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_boxcut',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_boxcut_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_boxcut_d
#endif
