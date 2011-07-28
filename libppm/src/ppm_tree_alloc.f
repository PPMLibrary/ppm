      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_tree_alloc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine (re)allocates the tree data
      !                 structures.
      !
      !  Input        : iopt         (I) Allocation mode (passed on to
      !                                  ppm_alloc)
      !                 nbox         (I) new number of boxes to allocate
      !                 nbpd         (I) number of children per parent
      !
      !                 The following argument is only present in the
      !                 __TREE case.
      !
      !                 nlevel       (I) number of levels in the tree
      !
      !  Input/output : min_box(:,:) (F) lower coordinates of the box.  
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 max_box(:,:) (F) upper coordinates of the box.
      !                                  1st index: x,y[,z], 2nd: box ID
      !                 boxcost(:)   (F) cost of all the boxes.
      !                 nchld(:)     (I) number of children of each box.
      !                 blevel(:)    (I) tree level of each box.
      !
      !                 The following arguments are optional, but all have
      !                 to be present if nlevel is present:
      !
      !                 parent(:)    (I) index of the parent box of each
      !                                  box. ppm_param_undefined if no 
      !                                  parent (i.e. root box)
      !                 child(:,:)   (I) indices of all children of a box.
      !                                  1st index: child ID, 2nd: box ID.
      !                 nbpl(:)      (I) the number of boxes per level
      !
      !  Output       : info         (I) return status.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_alloc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.9  2005/09/01 11:36:01  ivos
      !  Uses more accurate estimate for number of boxes. fixed alloc bug in the
      !  case where nbox.LT.nbpd.
      !
      !  Revision 1.8  2005/08/31 13:33:49  ivos
      !  bugfix: removed doubly-declared variables and unused arguments.
      !
      !  Revision 1.7  2005/08/31 11:24:30  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.6  2005/05/24 23:28:17  ivos
      !  Removed unused local variables.
      !
      !  Revision 1.5  2004/12/03 17:17:17  ivos
      !  Now only allocates what is needed.
      !
      !  Revision 1.4  2004/11/25 13:01:59  ivos
      !  bugfix: type in cpp directive fixed.
      !
      !  Revision 1.3  2004/11/04 12:57:02  ivos
      !  Added maxlevels to args of ppm_tree in order to be able to specify
      !  the maximum number of tree levels before done.
      !
      !  Revision 1.2  2004/09/22 17:24:48  ivos
      !  Added nchld also for TYPE==DECOMP since box2subs needs it.
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

#if   __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_alloc_ts(iopt,nbox,nbpd,nlevel,min_box,     &
     &    max_box,boxcost,parent,nchld,child,blevel,nbpl,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_alloc_td(iopt,nbox,nbpd,nlevel,min_box,     &
     &    max_box,boxcost,parent,nchld,child,blevel,nbpl,info)
#endif
#elif __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_alloc_ds(iopt,nbox,nbpd,min_box,max_box,    &
     &    boxcost,nchld,blevel,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_alloc_dd(iopt,nbox,nbpd,min_box,max_box,    &
     &    boxcost,nchld,blevel,info)
#endif
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
      REAL(MK), DIMENSION(:,:), POINTER       :: min_box,max_box
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      INTEGER , DIMENSION(:  ), POINTER       :: nchld
      INTEGER                 , INTENT(IN   ) :: iopt,nbox,nbpd
      INTEGER                 , INTENT(  OUT) :: info
      INTEGER , DIMENSION(:  ), POINTER       :: blevel
#if   __TYPE == __TREE
      INTEGER                 , INTENT(IN   ) :: nlevel
      INTEGER , DIMENSION(:  ), POINTER       :: parent,nbpl
      INTEGER , DIMENSION(:,:), POINTER       :: child
#endif
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK)                                :: t0
      INTEGER, DIMENSION(2)                   :: ldc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_alloc',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (nbox .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_alloc',     &
     &          'Number of boxes must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (nbpd .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_alloc',     &
     &          'Number of boxes per step must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
#if   __TYPE == __TREE
         IF (nlevel .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_alloc',     &
     &          'Number of levels must be >= 0',__LINE__,info)
            GOTO 9999
         ENDIF 
#endif
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          ldc(1) = 2
          ldc(2) = nbox
          CALL ppm_alloc(tree_lhbx,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &            'pointer to headers TREE_LHBX',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 
      ldc(1)   = ppm_dim
      ldc(2)   = nbox
      CALL ppm_alloc(min_box,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'lower box boundaries MIN_BOX',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(max_box,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'upper box boundaries MAX_BOX',__LINE__,info)
          GOTO 9999
      ENDIF 
      IF (have_mesh) THEN
          CALL ppm_alloc(Nm_box,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &            'box grid size NM_BOX',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 
      ldc(1) = nbox
      CALL ppm_alloc(ndiv,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'number of divisible directions NDIV',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(blevel,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'tree levels of boxes BLEVEL',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(boxcost,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'box costs BOXCOST',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(nchld,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'number of children NCHLD',__LINE__,info)
          GOTO 9999
      ENDIF 
#if   __TYPE == __TREE
      ldc(1) = nbpd
      ldc(2) = nbox
      CALL ppm_alloc(child,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'list of children CHILD',__LINE__,info)
          GOTO 9999
      ENDIF 
      ldc(1) = nbox
      CALL ppm_alloc(parent,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'parent pointer PARENT',__LINE__,info)
          GOTO 9999
      ENDIF 
      ldc(1) = nlevel
      CALL ppm_alloc(nbpl,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree_alloc',          &
     &        'number of boxes per level NBPL',__LINE__,info)
          GOTO 9999
      ENDIF 
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_alloc',t0,info)
      RETURN
#if   __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_ts
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_td
#endif
#elif __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_alloc_dd
#endif
#endif
