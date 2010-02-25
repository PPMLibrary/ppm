      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_tree_divcheck
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine checks how many dimensions of a box
      !                 are divisible.
      !
      !  Input        : min_box(:,:) (F) lower coordinates of the boxes.
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 max_box(:,:) (F) upper coordinates of the boxes.
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 nbox         (I) number of boxes to check
      !                 minboxsize(:)(F) minimum box size in all
      !                                  dimensions.
      !                 fixed(:)     (L) Flags telling which dimensions are
      !                                  fixed (i.e. must not be divided).
      !                 boxcost(:)   (F) costs associated with boxes
      !
      !  Input/output :                                            
      !
      !  Output       : ndiv(:)      (I) number of divisible directions of
      !                                  each box (1..nbox).
      !                 info         (I) return status.
      !
      !  Remarks      : A box with cost 0 is counted as non-divisible.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_divcheck.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.8  2006/09/05 08:01:28  pchatela
      !  Proper scaling for REAL comparisons
      !  Added module_alloc to ppm_decomp_boxsplit
      !
      !  Revision 1.7  2005/08/31 11:24:32  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.6  2005/08/30 13:17:27  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.5  2005/01/27 09:24:24  ivos
      !  minboxsize .EQ. 0 is now allowed.
      !
      !  Revision 1.4  2004/09/30 15:53:51  ivos
      !  fix: replaces array expression with explicit DO loop, so pgf90 does
      !  not crash any more when compiling with -O2 -fast.
      !
      !  Revision 1.3  2004/09/23 09:48:35  ivos
      !  introduced lmyeps for real comparisons.
      !
      !  Revision 1.2  2004/09/22 17:28:36  ivos
      !  fixed wrong error message strings.
      !
      !  Revision 1.1  2004/09/22 10:32:05  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_divcheck_s(min_box,max_box,nbox,minboxsize,   &
     &    fixed,boxcost,ndiv,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_divcheck_d(min_box,max_box,nbox,minboxsize,   &
     &    fixed,boxcost,ndiv,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box,max_box
      INTEGER                 , INTENT(IN   ) :: nbox
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: minboxsize,boxcost
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      INTEGER , DIMENSION(:  ), POINTER       :: ndiv
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,lmyeps,boxlen
      REAL(MK), DIMENSION(ppm_dim)            :: ms2
      INTEGER                                 :: iopt,i,j
      INTEGER, DIMENSION(2)                   :: ldc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_divcheck',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (nbox .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',     &
     &          'Number of boxes must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            IF (minboxsize(i) .LT. 0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',     &
     &             'the minimum box size must be > 0 !',__LINE__,info)
               GOTO 9999
            ENDIF 
            DO j=1,nbox
               IF (min_box(i,j) .GT. max_box(i,j)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_tree_divcheck',   &
     &                'min_box must be <= max_box !',__LINE__,info)
                  GOTO 9999
               ENDIF 
            ENDDO
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  If there are no boxes to check, we quit
      !-------------------------------------------------------------------------
      IF (nbox .LT. 1) THEN
          CALL ppm_write(ppm_rank,'ppm_tree_divcheck',   &
     &        'No boxes to be checked. Exiting.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Count and determine divisible dimensions
      !-------------------------------------------------------------------------
      DO i=1,ppm_dim
          ms2(i) = 2.0_MK*minboxsize(i)
      ENDDO

      DO i=1,nbox
          ndiv(i) = 0
          IF (boxcost(i) .GT. lmyeps) THEN
              DO j=1,ppm_dim
                  boxlen = max_box(j,i)-min_box(j,i)
                  IF (((boxlen-ms2(j)).GT.lmyeps*boxlen).AND.(.NOT.fixed(j))) THEN
                      ndiv(i) = ndiv(i) + 1
                  ENDIF
              ENDDO
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_divcheck',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_divcheck_d
#endif
