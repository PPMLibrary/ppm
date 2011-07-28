      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_done
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine decides if a decomposition is done or
      !                 one or more boxes need further refinement.  
      !
      !  Input        : minboxes     (I) minimum number of childless boxes 
      !                                  of non-zero cost to create.
      !                 nsubs        (I) number of childless boxes with
      !                                  non-zero cost. This is not the
      !                                  same as number of divisible boxes
      !                                  as they need not be larger than
      !                                  2*minboxsize in this case.
      !                 boxcost(:)   (F) computational cost associated with
      !                                  each box
      !                 iboxlist(:)  (I) list of divisible boxes.
      !                 nboxlist     (I) Number of boxes which could be
      !                                  divided further (i.e. length of
      !                                  iboxlist).
      !                 nlevel       (I) Number if tree levels (tree
      !                                  depth) so far.
      !                 maxvariance  (F) Maximum variance of cost allowed 
      !                                  between boxes. Set to .LE. 0 to
      !                                  disable this.
      !                 maxboxcost   (F) Maximum allowed cost of a box. Tree 
      !                                  will stop if all boxes are
      !                                  below this cost. If .LE. 0, cost
      !                                  is unlimited.
      !                 maxlevels    (I) Maximum number of levels to be
      !                                  created. Tree stops as soon as
      !                                  this is reached. If .LE. 0,
      !                                  levels are unlimited.
      !
      !  Input/output :                                            
      !
      !  Output       : lcontinue    (L) .FALSE. if no further subdivision
      !                                  is needed, .TRUE. otherwise.
      !                 info         (I) return status.
      !
      !  Remarks      : Decomposition is considered done if more boxes than
      !                 minboxes are present and the variance of the costs
      !                 of all divisible (i.e. childless and non-empty)
      !                 boxes is below maxvariance. Alternatively, a
      !                 decomposition is always done when no more divisible 
      !                 boxes exist or the specified maximum number of levels
      !                 has been reached.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_done.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2005/08/31 12:43:45  ivos
      !  Shark optimizations.
      !
      !  Revision 1.8  2005/08/31 11:24:32  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.7  2004/12/03 17:16:56  ivos
      !  Changed from ldone to lcontinue. Formulate positive sentences!
      !
      !  Revision 1.6  2004/12/02 10:00:04  ivos
      !  maxvariance is now allowed to be negative.
      !
      !  Revision 1.5  2004/11/30 15:05:20  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.4  2004/11/04 12:57:03  ivos
      !  Added maxlevels to args of ppm_tree in order to be able to specify
      !  the maximum number of tree levels before done.
      !
      !  Revision 1.3  2004/09/23 09:47:54  ivos
      !  bugfix: removed extra argument from ppm_write.
      !
      !  Revision 1.2  2004/09/22 17:29:04  ivos
      !  bugfix: added USE ppm_module_write.
      !
      !  Revision 1.1  2004/09/22 10:32:03  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_done_s(minboxes,nsubs,boxcost,iboxlist,nboxlist,   &
     &    nlevel,maxvariance,maxboxcost,maxlevels,lcontinue,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_done_d(minboxes,nsubs,boxcost,iboxlist,nboxlist,   &
     &    nlevel,maxvariance,maxboxcost,maxlevels,lcontinue,info)
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: boxcost
      INTEGER                 , INTENT(IN   ) :: nboxlist,minboxes,nsubs
      INTEGER                 , INTENT(IN   ) :: nlevel,maxlevels
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: iboxlist
      REAL(MK)                , INTENT(IN   ) :: maxvariance,maxboxcost
      LOGICAL                 , INTENT(  OUT) :: lcontinue
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0,meancost,diffcost,varcost
      REAL(MK)                                :: maxcost,dm
      INTEGER                                 :: i,j
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_done',t0,info)
      lcontinue = .TRUE.

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (nsubs .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_done',     &
     &          'Number of non-empty boxes must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (nlevel .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_done',     &
     &          'Number of levels must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (nboxlist .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_done',     &
     &          'Number of boxes in list must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  If there are no boxes, we quit
      !-------------------------------------------------------------------------
      IF (nsubs .LT. 1) THEN
          lcontinue = .FALSE. ! no boxes = nothing to subdivide !
          CALL ppm_write(ppm_rank,'ppm_tree_done',   &
     &        'No non-empty boxes present. Done.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are no more boxes that could potentially be divided, we
      !  have to stop.
      !-------------------------------------------------------------------------
      IF (nboxlist .LT. 1) THEN
          lcontinue = .FALSE.
          !---------------------------------------------------------------------
          !  If there are less boxes than processors, we have a problem
          !  THIS SHOULD GO TO THE DECOMP ROUTINE AND NOT INTO THE GENERIC
          !  TREE !!!
          !---------------------------------------------------------------------
          IF (nsubs .LT. minboxes) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_few_subs,'ppm_tree_done',     &
     &            'Could not create the minimum number of non-empty boxes!', &
     &            __LINE__,info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If the max tree depth has been reached, we stop
      !-------------------------------------------------------------------------
      IF (maxlevels .GT. 0) THEN
          IF (nlevel .GE. maxlevels) THEN
              lcontinue = .FALSE.
              IF (ppm_debug .GT. 1) THEN
                  CALL ppm_write(ppm_rank,'ppm_tree_done',     &
     &                'Max number of levels reached. Done.',info)
              ENDIF
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  If there are less boxes than processors, decomposition is not OK
      !-------------------------------------------------------------------------
      IF (nsubs .LT. minboxes) GOTO 9999

      !-------------------------------------------------------------------------
      !  Compute the variance of divisible box costs
      !-------------------------------------------------------------------------
      varcost  = 0.0_MK
      meancost = 0.0_MK
      maxcost  = -HUGE(maxcost)
      DO i=1,nboxlist
          j = iboxlist(i)
          dm = boxcost(j)
          meancost = meancost + dm
          IF (dm .GT. maxcost) THEN
              maxcost = dm
          ENDIF
      ENDDO
      IF (nboxlist .GT. 1) THEN
          meancost = meancost/REAL(nboxlist,MK)
          DO i=1,nboxlist
              j = iboxlist(i)
              diffcost = boxcost(j) - meancost
              varcost  = varcost + (diffcost*diffcost)
          ENDDO
          varcost  = varcost/REAL(nboxlist-1,MK)
      ENDIF
          
      !-------------------------------------------------------------------------
      !  If variance of costs is below threshold, decomposition is OK
      !-------------------------------------------------------------------------
      IF (varcost .LT. maxvariance) lcontinue = .FALSE.

      !-------------------------------------------------------------------------
      !  If all boxes are below the maxcost we are done.
      !-------------------------------------------------------------------------
      IF (maxcost .LT. maxboxcost) lcontinue = .FALSE. 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_done',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_done_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_done_d
#endif
