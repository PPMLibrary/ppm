      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_topo_box2subs
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine converts boxes from a ppm tree
      !                 (created by ppm_tree) to a set of subdomains which
      !                 define a valid domain decomposition. If no tree level
      !                 is specified, this is done by taking each childless
      !                 tree box as a sub. Otherwise the boxes of the
      !                 specific level are taken. If a certain branch of
      !                 the tree does not extend to that level, the next
      !                 higher existing box can be taken to fill the gap
      !                 in order to return a decomposition with no holes 
      !                 in space.
      !
      !  Input        : min_box(:,:) (F) lower coordinates of the boxes.
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 max_box(:,:) (F) upper coordinates of the boxes.
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 nchld(:)     (I) number of children of each box.
      !                 nbox         (I) total number of boxes 
      !                 blevel(:)    (I) tree level of each box as
      !                                  returned by ppm_tree. OPTIONAL.
      !                                  Needs to be present if level is
      !                                  present.
      !                 child(:,:)   (I) children (1st index) of all 
      !                                  boxes (2nd index). OPTIONAL.
      !                                  Needs to be present if level is
      !                                  present.
      !                 level        (I) OPTIONAL argument to specify
      !                                  tree level from which boxes are
      !                                  to be taken. If .GT. 0, only
      !                                  the boxes which are precisely
      !                                  on that level are returned,
      !                                  possibly leading to holes in
      !                                  the decomposition. If .LT. 0,
      !                                  holes are filled with the next
      !                                  higher-level existing box.
      !
      !  Input/output :                                            
      !
      !  Output       : min_sub(:,:) (F) lower coordinates of the subs.   
      !                 max_sub(:,:) (F) upper coordinates of the subs.
      !                 nsubs        (I) number of subs.
      !                 info         (I) return status.
      !                 boxid(:)     (I) to which of the original boxes
      !                                  does each sub correspond?
      !                                  OPTIONAL only allocated and
      !                                  returned of present.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_box2subs.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2005/02/01 15:30:10  ivos
      !  bugfix: size test for stack reallocation was wrong.
      !
      !  Revision 1.3  2005/01/24 08:43:43  ivos
      !  Added level,blevel,child OPTIONAL arguments to allow tree level to
      !  be specified instead of just leafs.
      !
      !  Revision 1.2  2004/11/30 15:06:29  ivos
      !  Added optional argument boxid and improved vectorization.
      !
      !  Revision 1.1  2004/09/23 10:00:04  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_box2subs_s(min_box,max_box,nchld,nbox,   &
     &    min_sub,max_sub,nsubs,info,boxid,level,blevel,child)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_box2subs_d(min_box,max_box,nchld,nbox,   &
     &    min_sub,max_sub,nsubs,info,boxid,level,blevel,child)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: nchld
      INTEGER , DIMENSION(:  ), INTENT(IN   ), OPTIONAL :: blevel
      INTEGER , DIMENSION(:,:), INTENT(IN   ), OPTIONAL :: child
      INTEGER                 , INTENT(IN   ), OPTIONAL :: level
      INTEGER , DIMENSION(:  ), POINTER, OPTIONAL :: boxid
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER                 , INTENT(  OUT) :: nsubs,info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER                                 :: iopt,i,j,istack
      INTEGER, DIMENSION(2)                   :: ldc
      INTEGER, DIMENSION(:), POINTER          :: subbox,boxstack
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_box2subs',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (nbox .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_topo_box2subs',     &
     &          'Number of boxes must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            DO j=1,nbox
               IF (min_box(i,j) .GT. max_box(i,j)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_box2subs',   &
     &                'min_box must be <= max_box !',__LINE__,info)
                  GOTO 9999
               ENDIF 
            ENDDO
         ENDDO
         IF (PRESENT(level)) THEN
             IF (.NOT.PRESENT(blevel)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_box2subs',     &
     &              'blevel must be present if level is.',__LINE__,info)
                GOTO 9999
             ENDIF
             IF (.NOT.PRESENT(child)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_box2subs',     &
     &              'child must be present if level is.',__LINE__,info)
                GOTO 9999
             ENDIF
         ENDIF
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate memory for list
      !-------------------------------------------------------------------------
      iopt     = ppm_param_alloc_fit
      ldc(1)   = nbox
      CALL ppm_alloc(subbox,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &        'list of box IDs SUBBOX',__LINE__,info)
          GOTO 9999
      ENDIF 
      subbox   = ppm_param_undefined

      !-------------------------------------------------------------------------
      !  Count the number of subs and build a list
      !-------------------------------------------------------------------------
      nsubs = 0
      IF (PRESENT(level)) THEN
          IF (level .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Count boxes on specified level
              !-----------------------------------------------------------------
              DO i=1,nbox
                  IF (blevel(i) .EQ. level) THEN
                      nsubs = nsubs + 1
                      subbox(nsubs) = i
                  ENDIF
              ENDDO
          ELSE
              !-----------------------------------------------------------------
              !  Traverse tree at most down to the specified level
              !-----------------------------------------------------------------
              iopt     = ppm_param_alloc_fit
              ldc(1)   = MIN((-level-2)*(nchld(1)-1) + nchld(1),1)
              CALL ppm_alloc(boxstack,ldc,iopt,info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &                'stack of boxes to traverse BOXSTACK',__LINE__,info)
                  GOTO 9999
              ENDIF 

              ! push root box
              istack = 1
              boxstack(istack) = 1

              ! traverse the tree to find boxes
              DO WHILE (istack .GT. 0) 
                  i = boxstack(istack)
                  istack = istack - 1
                  ! add box if it is on desired level or above and has
                  ! no children
                  IF ((blevel(i).EQ.-level).OR.(nchld(i).EQ.0)) THEN
                      nsubs = nsubs + 1
                      subbox(nsubs) = i
                  ELSE
                      ! grow stack if needed
                      IF ((istack+nchld(i)) .GT. ldc(1)) THEN
                          ldc(1) = ldc(1) + nchld(i)
                          iopt     = ppm_param_alloc_grow_preserve
                          CALL ppm_alloc(boxstack,ldc,iopt,info)
                          IF (info.NE.0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',&
     &                            'stack of boxes to traverse BOXSTACK',   &
     &                            __LINE__,info)
                              GOTO 9999
                          ENDIF 
                      ENDIF 
                      ! push all its children to the stack
                      DO j=1,nchld(i)
                          istack = istack + 1
                          boxstack(istack) = child(j,i)
                      ENDDO
                  ENDIF
              ENDDO
              
              ! deallocate stack
              iopt     = ppm_param_dealloc
              CALL ppm_alloc(boxstack,ldc,iopt,info)
              IF (info.NE.0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_dealloc,'ppm_topo_box2subs',        &
     &                'stack of boxes to traverse BOXSTACK',__LINE__,info)
              ENDIF 
          ENDIF
      ELSE
          !---------------------------------------------------------------------
          !  Count childless boxes
          !---------------------------------------------------------------------
          DO i=1,nbox
              IF (nchld(i) .LT. 1) THEN
                  nsubs = nsubs + 1
                  subbox(nsubs) = i
              ENDIF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for subs
      !-------------------------------------------------------------------------
      iopt     = ppm_param_alloc_fit
      ldc(1)   = ppm_dim
      ldc(2)   = nsubs
      CALL ppm_alloc(min_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &        'minimum extent of subs MIN_SUB',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(max_sub,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &        'maximum extent of subs MAX_SUB',__LINE__,info)
          GOTO 9999
      ENDIF 
      IF (PRESENT(boxid)) THEN
          ldc(1) = nsubs
          CALL ppm_alloc(boxid,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &            'box ID of subs BOXID',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Store subs
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO i=1,nsubs
              j = subbox(i)
              min_sub(1,i) = min_box(1,j)
              min_sub(2,i) = min_box(2,j)
              max_sub(1,i) = max_box(1,j)
              max_sub(2,i) = max_box(2,j)
          ENDDO
      ELSE
          DO i=1,nsubs
              j = subbox(i)
              min_sub(1,i) = min_box(1,j)
              min_sub(2,i) = min_box(2,j)
              min_sub(3,i) = min_box(3,j)
              max_sub(1,i) = max_box(1,j)
              max_sub(2,i) = max_box(2,j)
              max_sub(3,i) = max_box(3,j)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Store box IDs if needed
      !-------------------------------------------------------------------------
      IF (PRESENT(boxid)) THEN
          DO i=1,nsubs
              boxid(i) = subbox(i)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate list
      !-------------------------------------------------------------------------
      iopt     = ppm_param_dealloc
      CALL ppm_alloc(subbox,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_alloc,'ppm_topo_box2subs',          &
     &        'list of box IDs SUBBOX',__LINE__,info)
      ENDIF 

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_box2subs',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_box2subs_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_box2subs_d
#endif
