      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_tree
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs a generic tree decomposition
      !                 using recursive orthogonal multisection (ROM) in either
      !                 2, 4, or 8 subunits (binary tree, quadtree,
      !                 octtree). The tree can be based on scattered data
      !                 points (if Np.GT.0), on regular mesh points (if
      !                 Nm.GT.0) or purely on the geometry of the domain
      !                 (if both Np and Nm are 0).
      !
      !                 If preprocessed with __TYPE == __TREE, the full
      !                 tree information (connectivity, levels, etc.) is
      !                 computed and returned. If __TYPE == __DECOMP, only
      !                 the space decomposition is done, but no tree
      !                 information returned.
      !
      !  Input        : xp(:,:)      (F) the data points
      !                 Np           (I) the number of data points. If .LE.
      !                                  0, decomposition is based on
      !                                  geometry and mesh only.
      !                 Nm(:)        (I) number of grid points in the
      !                                  global mesh. (0,0,0) if there is
      !                                  no mesh. If a mesh is present, the
      !                                  box boundaries will be aligned
      !                                  with mesh planes.
      !                 min_dom(:)   (F) the minimum coordinate of the 
      !                                  domain 
      !                 max_dom(:)   (F) the maximum coordinate of the 
      !                                  domain 
      !                 treetype     (I) type of multisection tree. One of:
      !                                     ppm_param_tree_bin
      !                                     ppm_param_tree_quad
      !                                     ppm_param_tree_oct (3D only)
      !                                  For binary, quad- or oct-tree. 
      !                 minboxes     (I) minimum number of childless boxes 
      !                                  (leaves) of non-zero cost to be 
      !                                  created. Set this to -1 if there 
      !                                  is no minimum requirement.
      !                 pruneboxes   (L) .TRUE. to prune the tree to only
      !                                  contain boxes of non-zero cost.
      !                                  .FALSE. to get a tree with all
      !                                  boxes (also the empty ones).
      !                 minboxsize(:)(F) the miminum box size in all
      !                                  directions
      !                 maxvariance  (F) maximum variance of cost allowed
      !                                  between boxes. The tree stops as
      !                                  soon as the variance of costs of
      !                                  all boxes is below this max.
      !                                  Set this to -1 to disable this
      !                                  criterion.
      !                 maxboxcost   (F) Maximum cost per
      !                                  box. Subdivision will stop
      !                                  as soon as all boxes have costs 
      !                                  below this value. Set this to -1
      !                                  to not impose any limit. 
      !                 maxlevels    (I) Maximum number of levels to
      !                                  create. Tree stops as soon as
      !                                  this is reached. The root box
      !                                  is counted as level 1. Set to .LE. 0
      !                                  for unlimited levels. This
      !                                  input is only present in the
      !                                  __TREE version, not in the
      !                                  __DECOMP version.
      !                 fixed(1:dim) (L) Flag which tells for each spatial
      !                                  dimension (1..ppm_dim) if it is
      !                                  fixed (i.e. no cuts perpendicular
      !                                  to it are allowed). fixed=(F,F,T)
      !                                  will e.g. enforce z-pencil
      !                                  decompositions.
      !                 weights(3,2) (F) weights for the three cost
      !                                  contributions (particles, mesh
      !                                  points, volume) (1st index) for 
      !                                  box cost (weights(:,1)) and the 
      !                                  determination of the cut planes
      !                                  (weights(:,2)).
      !                 pcost(:)     (F) OPTIONAL argument of length
      !                                  Np, specifying the
      !                                  cost of each data point.
      !
      !  Input/output :                                            
      !
      !  Output       : min_box(:,:) (F) the min. extents of the boxes
      !                 max_box(:,:) (F) the max. extents of the boxes
      !                 nbox         (I) the total number of boxes
      !                 nchld(:)     (I) number of children of each box.
      !                 info         (I) return status
      !      
      !                 The following outputs are optional. If these
      !                 arguments are present, they will be computed and
      !                 returned, otherwise not:
      !
      !                 lhbx(:,:)    (I) pointer to first (1,:) and last
      !                                  (2,:) data point (in
      !                                  lpdx) in each tree box.
      !                                  This is only present if
      !                                  Np.GT.0. Only points on the
      !                                  local processor are considered.
      !                                  Entries for non-leaf boxes are
      !                                  only true if pruneboxes=.FALSE.
      !                 lpdx(:)      (I) pointers to the data points (in
      !                                  xp) in each tree box. Box ib
      !                                  conatins points
      !                                    lpdx(lhbx(1,ib):(lhbx(2,ib)))
      !                                  This is only present if Np.GT.0. 
      !                                  Only points on the local
      !                                  processor are considered.
      !                                  Entries for non-leaf boxes are
      !                                  only true if pruneboxes=.FALSE.
      !                 boxcost(:)   (F) costs of all boxes 1..nbox.
      !                 parent(:)    (I) index of the parent box of each
      !                                  box. ppm_param_undefined if no 
      !                                  parent (i.e. root box)
      !                 child(:,:)   (I) indices of all children of a box.
      !                                  1st index: child ID, 2nd: box ID.
      !                 blevel(:)    (I) tree level of each box. 1..nbox.
      !                                  Level 1 is the root box.
      !                 nbpl(:)      (I) the number of boxes per level.
      !                                  Level 1 is the root box.
      !                 nlevel       (I) the number of levels. Level 1 is
      !                                  the root box.
      !
      !  Remarks      : If Np.GT.0, the particle-caused cost per sub is
      !                 given by number(particles) [or sum(pcost)] in that
      !                 sub.
      !                 If Nm(1:ppm_dim).GT.1, the mesh-causes cost is
      !                 given by the number of mesh points in each sub.
      !                 The geometry-based cost is always given by the
      !                 volume of each sub.
      !
      !                 These costs can be freely composed using the
      !                 weights argument.
      !
      !                 ppm_decomp_cartesian is not made obsolete by this
      !                 routine since the latter can gracefully produce ANY
      !                 (even prime) number of equisized subdomains,
      !                 whereas this routine is restricted to powers of
      !                 2, 4, or 8 (for equisized subs!).
      !        
      !                 directly compute box costs in tree_cutpos
      !                 (since we do the allreduce there anyway) and return
      !                 them. only recompute (with tree_boxcost) if the cut
      !                 planes were shifted.
      !
      !                 speed it up by mapping the particles once there are
      !                 more boxes than processors. each processor then
      !                 computes its own subtree. The trees are reduced and
      !                 merged at the end.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.29  2006/11/15 13:50:27  pchatela
      !  Re-disabled the ppm_decomp_tree, which does not have index arrays yet
      !
      !  Revision 1.28  2006/09/04 18:34:56  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.27  2005/09/06 14:47:37  ivos
      !  bugfix: for Np.LT.3 the alloc failed because nlevelalloc was
      !  estimated to 0. Now forced to a minimum of 1.
      !
      !  Revision 1.26  2005/09/01 11:36:00  ivos
      !  Uses more accurate estimate for number of boxes. fixed alloc bug in the
      !  case where nbox.LT.nbpd.
      !
      !  Revision 1.25  2005/08/31 13:33:48  ivos
      !  bugfix: removed doubly-declared variables and unused arguments.
      !
      !  Revision 1.24  2005/08/31 12:43:44  ivos
      !  Shark optimizations.
      !
      !  Revision 1.23  2005/08/31 11:24:30  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.22  2005/08/30 13:17:25  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.21  2005/08/30 12:29:44  ivos
      !  Removed debug PRINT.
      !
      !  Revision 1.20  2005/08/30 12:27:54  ivos
      !  Added estimate for nlevels and nboxes for allocation so that the
      !  lists are not grown at every step. Speedup of a factor of 1000!
      !
      !  Revision 1.19  2005/02/01 13:22:57  ivos
      !  Moved declarations of lhbx_cut and lpdx_cut to module_data_tree.
      !
      !  Revision 1.18  2005/01/31 08:23:36  ivos
      !  Fixed declaration of mhbx_cut
      !
      !  Revision 1.17  2005/01/27 09:24:24  ivos
      !  minboxsize .EQ. 0 is now allowed.
      !
      !  Revision 1.16  2004/12/09 20:23:29  ivos
      !  bugfix: forgot to initialize nbpl after (re)allocation.
      !
      !  Revision 1.15  2004/12/03 17:43:18  ivos
      !  Removed debug output.
      !
      !  Revision 1.14  2004/12/03 17:16:03  ivos
      !  Now uses and returns particle-in-box lists lhbx and lpdx.
      !  Weights was changed to be 2d array to allow different cost
      !  definitions for refinementd and cutting.
      !
      !  Revision 1.13  2004/12/02 10:01:14  ivos
      !  Added check of info after ppm_tree_done calls.
      !
      !  Revision 1.12  2004/11/30 15:26:52  ivos
      !  Added TODO comment.
      !
      !  Revision 1.11  2004/11/30 15:05:20  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.10  2004/11/04 12:57:02  ivos
      !  Added maxlevels to args of ppm_tree in order to be able to specify
      !  the maximum number of tree levels before done.
      !
      !  Revision 1.9  2004/11/03 16:26:26  ivos
      !  bugfix: info is now correctly propagared out if errors occur.
      !
      !  Revision 1.8  2004/10/27 11:48:31  ivos
      !  bugfix: wrong index in boxcost fixed for the pruneboxes case.
      !
      !  Revision 1.7  2004/09/30 10:49:42  ivos
      !  bugfix: SIZE of Nm is now checked before accessing it when checking
      !  if there is a mesh at all. This allows the user to pass dummy pointers.
      !
      !  Revision 1.6  2004/09/29 11:05:42  ivos
      !  bugfix: in the case of a mesh cut failure, it was not checked if
      !  tree_done. This could result in a length 0 boxlist being attempted
      !  to be used. Fixed by checking after the failure.
      !
      !  Revision 1.5  2004/09/24 14:59:27  ivos
      !  Added initialization for have_mesh and have_particles.
      !
      !  Revision 1.4  2004/09/24 08:02:15  ivos
      !  bugfix: nchld and child of newly created boxes are now properly
      !  initialized. Mesh cuts cannot create boxes smaller than the
      !  required minimum any more.
      !
      !  Revision 1.3  2004/09/23 09:51:30  ivos
      !  bugfix: introduces inextboxlist since this is the box which needs
      !  to be replaced and NOT inext.
      !
      !  Revision 1.2  2004/09/22 17:25:24  ivos
      !  Added nchld also for TYPE==DECOMP since box2subs needs it.
      !
      !  Revision 1.1  2004/09/22 10:32:02  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_ds(xp,Np,Nm,min_dom,max_dom,treetype,     &
     &   minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,     &
     &   fixed,weights,min_box,max_box,nbox,nchld,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_dd(xp,Np,Nm,min_dom,max_dom,treetype,     &
     &   minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,     &
     &   fixed,weights,min_box,max_box,nbox,nchld,info,pcost)
#endif
#elif __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_ts(xp,Np,Nm,min_dom,max_dom,treetype,            &
     &   minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,  &
     &   fixed,weights,min_box,max_box,lhbx,lpdx,boxcost,        &
     &   parent,nchld,child,blevel,nbox,nbpl,nlevel,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_td(xp,Np,Nm,min_dom,max_dom,treetype,            &
     &   minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,  &
     &   fixed,weights,min_box,max_box,lhbx,lpdx,boxcost,        &
     &   parent,nchld,child,blevel,nbox,nbpl,nlevel,info,pcost)
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
      USE ppm_module_write
      USE ppm_module_tree_alloc
      USE ppm_module_tree_divcheck
      USE ppm_module_tree_cutdir
      USE ppm_module_tree_cutpos
      USE ppm_module_tree_boxcut
      USE ppm_module_tree_boxcost
      USE ppm_module_tree_done
      USE ppm_module_util_rank
      USE ppm_module_decomp_tree
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
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_dom,max_dom,minboxsize
      REAL(MK), DIMENSION(3,2), INTENT(IN   ) :: weights
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      LOGICAL                 , INTENT(IN   ) :: pruneboxes
      REAL(MK), DIMENSION(:,:), POINTER       :: min_box,max_box
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(IN   ) :: Np,treetype,minboxes
      REAL(MK)                , INTENT(IN   ) :: maxvariance
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      INTEGER                 , INTENT(  OUT) :: nbox,info
      INTEGER , DIMENSION(:  ), POINTER       :: nchld
#if   __TYPE == __TREE
      INTEGER                 , INTENT(IN   ) :: maxlevels
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      INTEGER , DIMENSION(:  ), POINTER       :: parent,nbpl,blevel
      INTEGER , DIMENSION(:,:), POINTER       :: lhbx
      INTEGER , DIMENSION(:  ), POINTER       :: lpdx
      INTEGER , DIMENSION(:,:), POINTER       :: child
      INTEGER                 , INTENT(  OUT) :: nlevel
#endif
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: mins,maxs,meshdx,meshdxinv
      INTEGER , DIMENSION(ppm_dim)            :: thisNm
      INTEGER , DIMENSION(2*ppm_dim)          :: ghostNm
      INTEGER , DIMENSION(2)                  :: ldc
      INTEGER , DIMENSION(1)                  :: ldc1
      INTEGER , DIMENSION(3)                  :: cutdir
      REAL(MK), DIMENSION(3)                  :: cutpos,boxlen
      REAL(MK), DIMENSION(:  ), POINTER       :: cpos,costc
      INTEGER , DIMENSION(:  ), POINTER       :: icut
      REAL(MK), DIMENSION(:,:), POINTER       :: minc,maxc
#if   __TYPE == __DECOMP
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      INTEGER , DIMENSION(:  ), POINTER       :: blevel
      INTEGER                                 :: nlevel
#endif
      INTEGER                                 :: nboxlist,nadd,k2,itype
      INTEGER                                 :: nboxalloc,nlevelalloc
      REAL(MK)                                :: t0,lmyeps,r0,r1,maxcost
      LOGICAL                                 :: up,nofixed,simpleweights
      LOGICAL                                 :: lcontinue
      INTEGER                                 :: i,j,k,l,iopt,inext,ncut,nbpd
      INTEGER                                 :: ibox,nboxold,nsubs
      INTEGER                                 :: inextboxlist,lctr
      INTEGER                                 :: nboxlistalloc
      INTEGER                                 :: info2,mxlev,bpc,istart,iend
      CHARACTER(LEN=ppm_char)                 :: mesg
#ifdef __MPI
      INTEGER                                 :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      lcontinue = .TRUE.
      itype = treetype
#if   __TYPE == __TREE
      mxlev = maxlevels
#else
      mxlev = -1
#endif

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
        IF (weights(1,1).EQ.0.0_MK.AND.weights(2,1).EQ.0.0_MK.AND.   &
     &      weights(3,1).EQ.0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree',    &
     &          'At least one weights(:,1) must be non-zero!',__LINE__,info)
            GOTO 9999
         ENDIF 
        IF (weights(1,2).EQ.0.0_MK.AND.weights(2,2).EQ.0.0_MK.AND.   &
     &      weights(3,2).EQ.0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree',    &
     &          'At least one weights(:,2) must be non-zero!',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (treetype.EQ.ppm_param_tree_oct.AND.ppm_dim.EQ.2) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_argument,'ppm_tree',    &
     &          'Octtree is not possible in 2d. Reverting to quadtree.',  &
     &          __LINE__,info)
            itype = ppm_param_tree_quad
         ENDIF 
         DO i=1,ppm_dim
            IF (minboxsize(i) .LT. 0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree',     &
     &             'the minimum box size must be > 0 !',__LINE__,info)
               GOTO 9999
            ENDIF 
            IF (min_dom(i) .GT. max_dom(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree',   &
     &             'min_dom must be <= max_dom !',__LINE__,info)
               GOTO 9999
            ENDIF 
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  Check what kind of input data is given
      !-------------------------------------------------------------------------
      have_particles = .FALSE.
      have_mesh      = .FALSE.
      IF (Np .GT. 0) have_particles = .TRUE.
      IF (SIZE(Nm,1) .GE. ppm_dim) THEN
          IF (ppm_dim .GT. 2) THEN
              IF ((Nm(1).GT.1).AND.(Nm(2).GT.1).AND.(Nm(3).GT.1))      &
     &            have_mesh = .TRUE.
          ELSE
              IF ((Nm(1).GT.1).AND.(Nm(2).GT.1)) have_mesh = .TRUE.
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Revert to the more efficient (specialized) decomposition trees where 
      !  possible
      !  PC: Disabled it: not fully compatible yet.
      !  Need to generate index lists etc...
      !-------------------------------------------------------------------------
      nofixed = .TRUE.
      DO i=1,ppm_dim
          IF (fixed(i)) nofixed = .FALSE.
      ENDDO
      simpleweights = .TRUE.
      IF (weights(1,1) .NE. 1.0_MK) simpleweights = .FALSE.
      IF (weights(2,1) .NE. 0.0_MK) simpleweights = .FALSE.
      IF (weights(3,1) .NE. 0.0_MK) simpleweights = .FALSE.
      IF (weights(1,2) .NE. 0.0_MK) simpleweights = .FALSE.
      IF (weights(2,2) .NE. 0.0_MK) simpleweights = .FALSE.
      IF (weights(3,2) .NE. 1.0_MK) simpleweights = .FALSE.

      !IF ((itype .EQ. ppm_param_tree_quad) .AND. (.NOT.have_mesh) .AND.  &
     !&    (have_particles) .AND. (.NOT.pruneboxes) .AND. (nofixed) .AND. &
     !&    (simpleweights) .AND. (ppm_dim .EQ. 2)) THEN
     !     IF (PRESENT(pcost)) THEN
     !         CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&          minboxsize(2)),maxvariance,min_box,max_box,nbox,info,pcost)
     !     ELSE
     !         CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&          minboxsize(2)),maxvariance,min_box,max_box,nbox,info)
     !     ENDIF
     !     GOTO 8000
     !ENDIF
     !IF ((itype .EQ. ppm_param_tree_oct) .AND. (.NOT.have_mesh) .AND.   &
     !&    (have_particles) .AND. (.NOT.pruneboxes) .AND. (nofixed) .AND. &
     !&    (simpleweights) .AND. (ppm_dim .EQ. 3)) THEN
     !     IF (PRESENT(pcost)) THEN
     !         CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&          minboxsize(2),minboxsize(3)),maxvariance,min_box,max_box,   &
     !&          nbox,info,pcost)
     !     ELSE
     !         CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&          minboxsize(2),minboxsize(3)),maxvariance,min_box,max_box,   &
     !&          nbox,info)
     !     ENDIF
     !     GOTO 8000
     !ENDIF

      !-------------------------------------------------------------------------
      !  Store the mesh spacings (if needed)
      !-------------------------------------------------------------------------
      IF (have_mesh) THEN
          meshdx(1) = (max_dom(1) - min_dom(1))/REAL(Nm(1)-1,MK)
          meshdxinv(1) = 1.0_MK/meshdx(1)
          meshdx(2) = (max_dom(2) - min_dom(2))/REAL(Nm(2)-1,MK)
          meshdxinv(2) = 1.0_MK/meshdx(2)
          IF (ppm_dim .GT. 2) THEN
              meshdx(3) = (max_dom(3) - min_dom(3))/REAL(Nm(3)-1,MK)
              meshdxinv(3) = 1.0_MK/meshdx(3)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the number of boxes and cuts per subdivision
      !-------------------------------------------------------------------------
      IF (itype .EQ. ppm_param_tree_bin) THEN
          nbpd = 2
          ncut = 1
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree','Creating binary tree.',info)
          ENDIF
      ELSEIF (itype .EQ. ppm_param_tree_quad) THEN
          nbpd = 4
          ncut = 2
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree','Creating quad-tree.',info)
          ENDIF
      ELSEIF (itype .EQ. ppm_param_tree_oct) THEN
          nbpd = 8
          ncut = 3
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree','Creating oct-tree.',info)
          ENDIF
      ELSE
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,'ppm_tree',     &
     &        'unknown tree type specified !',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Clear module pointers
      !-------------------------------------------------------------------------
      NULLIFY(tree_lhbx)
      NULLIFY(tree_lpdx)

      !-------------------------------------------------------------------------
      !  Allocate tree data structures
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      !-------------------------------------------------------------------------
      !  Guess tree size based on assumed uniform distribution
      !-------------------------------------------------------------------------
      nbox        = 1
      nlevel      = 1
      IF (Np .GT. 0) THEN
          IF (maxboxcost .GT. 0.0_MK) THEN
              ! the number of levels needed under uniform particle
              ! distribution
              nlevelalloc = CEILING(LOG(REAL(Np,MK)/maxboxcost)/   &
     &            LOG(REAL(nbpd,MK)))
#if __TYPE ==  __TREE
          ELSEIF (maxlevels .GT. 0) THEN
              ! assume we hit maxlevels
              nlevelalloc = maxlevels
#endif
          ELSE
              ! default assumption
              nlevelalloc = 3
          ENDIF
      ELSE
          ! Assume 3 levels if no particles are present
          nlevelalloc = 3
      ENDIF
      ! the number of boxes (not only leafs) is the geometric series sum
      nlevelalloc = nlevelalloc + 1   ! we start counting levels at 1
      IF (nlevelalloc .LT. 1) nlevelalloc = 1
      nboxalloc   = (1-(nbpd**nlevelalloc))/(1-nbpd)
      IF (nboxalloc .LT. 1) nboxalloc = 1
      IF (ppm_debug .GT. 0) THEN
#if   __TYPE == __TREE
          WRITE(mesg,'(A,I3,A,I6,A)') 'Allocating ',nlevelalloc,   &
     &        ' levels and ',nboxalloc,' boxes.'
#else
          WRITE(mesg,'(A,I3,A)') 'Allocating ',nboxalloc,' boxes.'
#endif
          CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
      ENDIF
#if   __TYPE == __TREE
      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,nlevelalloc,min_box,max_box,   &
     &    boxcost,parent,nchld,child,blevel,nbpl,info)
#elif __TYPE == __DECOMP
      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,min_box,max_box,   &
     &    boxcost,nchld,blevel,info)
#endif
      IF (info .NE. ppm_param_success) GOTO 9999

      !-------------------------------------------------------------------------
      !  Allocate memory for box cut
      !-------------------------------------------------------------------------
      ldc(1) = ppm_dim
      ldc(2) = 2**ncut
      CALL ppm_alloc(minc,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',          &
     &        'lower coordinates of new boxes MINC',__LINE__,info)
          GOTO 9999
      ENDIF 
      CALL ppm_alloc(maxc,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',          &
     &        'upper coordinates of new boxes MAXC',__LINE__,info)
          GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate local data structures
      !-------------------------------------------------------------------------
      nboxlist = 1
      nboxlistalloc = nbpd**(nlevelalloc-1)   ! list of leaves. guess...
      ldc(1) = nboxlistalloc
      CALL ppm_alloc(boxlist,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',          &
     &        'list of divisible boxes BOXLIST',__LINE__,info)
          GOTO 9999
      ENDIF 
      boxlist(1) = 1
      IF (have_mesh) THEN
          ldc(1) = ppm_dim
          ldc(2) = nbpd
          CALL ppm_alloc(Nmc,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree',          &
     &            'list of divisible boxes BOXLIST',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 
      IF (have_particles) THEN
          ldc(1) = 2**ncut
          CALL ppm_alloc(cbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree',   &
     &            'temporary box pointers CBOX',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(npbx,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree',   &
     &            'number of particles per box NPBX',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF 

      !-------------------------------------------------------------------------
      !  The domain itself is the root box. Get the tree started!
      !-------------------------------------------------------------------------
      IF (ppm_dim .GT. 2) THEN
          min_box(1,1) = min_dom(1)
          min_box(2,1) = min_dom(2)
          min_box(3,1) = min_dom(3)
          max_box(1,1) = max_dom(1)
          max_box(2,1) = max_dom(2)
          max_box(3,1) = max_dom(3)
          IF (have_mesh) THEN
              Nm_box(1,1) = Nm(1)
              Nm_box(2,1) = Nm(2)
              Nm_box(3,1) = Nm(3)
          ENDIF
      ELSE
          min_box(1,1) = min_dom(1)
          min_box(2,1) = min_dom(2)
          max_box(1,1) = max_dom(1)
          max_box(2,1) = max_dom(2)
          IF (have_mesh) THEN
              Nm_box(1,1) = Nm(1)
              Nm_box(2,1) = Nm(2)
          ENDIF
      ENDIF
      nsubs     = 1
      nchld     = 0
      blevel(1) = 1
#if   __TYPE == __TREE
      parent(1) = ppm_param_undefined
      child     = ppm_param_undefined
      nbpl(1)   = 1
#endif

      !-------------------------------------------------------------------------
      !  Rank the particles in the root box
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          IF (ppm_dim .EQ. 2) THEN
              thisNm(1)    = 1
              thisNm(2)    = 1
              ghostNm(1:4) = 0
              CALL ppm_util_rank2d(xp,Np,min_dom,max_dom,thisNm,ghostNm,     &
     &            tree_lpdx,lhbx_cut,info)
          ELSE
              thisNm(1)    = 1
              thisNm(2)    = 1
              thisNm(3)    = 1
              ghostNm(1:6) = 0
              CALL ppm_util_rank3d(xp,Np,min_dom,max_dom,thisNm,ghostNm,     &
     &            tree_lpdx,lhbx_cut,info)
          ENDIF
          IF (info .NE. ppm_param_success) GOTO 9999
          tree_lhbx(1,1) = lhbx_cut(1)
          tree_lhbx(2,1) = lhbx_cut(2) - 1
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the work memory for ppm_tree_boxcost. Do this once outside
      !  of the loop. This assumes that only the newly created boxes are
      !  tested, thus the size of nbpd.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nbpd
      CALL ppm_alloc(costc,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',     &
     &        'box costs COSTC',__LINE__,info)
          GOTO 9999
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcst_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcst_d,ldc,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',     &
     &        'particle cost part PCOST',__LINE__,info)
          GOTO 9999
      ENDIF
#ifdef __MPI
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcsum_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcsum_d,ldc,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',     &
     &        'particle cost sums PCSUM',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Compute cost of root box
      !-------------------------------------------------------------------------
      IF (PRESENT(pcost)) THEN
          CALL ppm_tree_boxcost(Nm_box,weights(:,1),min_box,max_box,  &
     &        1,lhbx_cut,tree_lpdx,boxcost,info,pcost)
      ELSE
          CALL ppm_tree_boxcost(Nm_box,weights(:,1),min_box,max_box,  &
     &        1,lhbx_cut,tree_lpdx,boxcost,info)
      ENDIF
      IF (info .NE. ppm_param_success) GOTO 9999
      
      !-------------------------------------------------------------------------
      !  Grow the list to the proper size as util_rank has only allocated
      !  it to length 2
      !-------------------------------------------------------------------------
      IF (have_particles) THEN
          iopt = ppm_param_alloc_grow
          ldc(1) = 2**ncut + 1
          CALL ppm_alloc(lhbx_cut,ldc,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_tree',          &
     &            'particle list header pointers LHBX_CUT',__LINE__,info)
              GOTO 9999
          ENDIF 
      ENDIF 

      !-------------------------------------------------------------------------
      !  Check if there is anything to be done at all
      !-------------------------------------------------------------------------
      CALL ppm_tree_done(minboxes,nsubs,boxcost,boxlist,nboxlist,   &
     &    nlevel,maxvariance,maxboxcost,mxlev,lcontinue,info)
      IF (info .NE. ppm_param_success) GOTO 9999
      IF ((.NOT.lcontinue) .AND. (ppm_debug .GT. 0)) THEN
          CALL ppm_write(ppm_rank,'ppm_tree',     &
     &        'Nothing to be done. Exiting.',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that root box is divisible
      !-------------------------------------------------------------------------
      CALL ppm_tree_divcheck(min_box,max_box,1,minboxsize,fixed,     &
     &    boxcost,ndiv,info)
      IF (info .NE. 0) GOTO 9999
      IF (ndiv(1) .LT. ncut) THEN
          lcontinue = .FALSE.
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree',     &
     &            'Initial domain is not divisible. Done.',info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for cut directions and positions
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = ncut
      CALL ppm_alloc(icut,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',     &
     &        'list of cut directions ICUT',__LINE__,info)
          GOTO 9999
      ENDIF
      icut = ppm_param_undefined
      CALL ppm_alloc(cpos,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_tree',     &
     &        'list of cut positions CPOS',__LINE__,info)
          GOTO 9999
      ENDIF
      cpos = -HUGE(cpos(1))

      !-------------------------------------------------------------------------
      !  Subdivide until done
      !-------------------------------------------------------------------------
      lctr = 0
      DO WHILE (lcontinue)
          lctr = lctr + 1
      
!         WRITE(mesg,'(a,i4.4)') 'boxes',lctr
!         OPEN(10,FILE=mesg)
!         DO i=1,nbox
! x-y plan
!           WRITE(10,'(2e12.4)') min_box(1,i),min_box(2,i)
!           WRITE(10,'(2e12.4)') max_box(1,i),min_box(2,i)
!           WRITE(10,'(2e12.4)') max_box(1,i),max_box(2,i)
!           WRITE(10,'(2e12.4)') min_box(1,i),max_box(2,i)
!           WRITE(10,'(2e12.4)') min_box(1,i),min_box(2,i)
!           WRITE(10,'(   a  )')
!         ENDDO
!         CLOSE(10)

          !---------------------------------------------------------------------
          !  Choose next subdomain to refine. This is always the one with
          !  maximum cost. In the case of particle, cost is number of
          !  particles, for meshes the number of mesh points and for
          !  geometric decompositions the sub volume.
          !---------------------------------------------------------------------
          maxcost = -HUGE(maxcost)
          inextboxlist = -1
          inext = -1
          DO i=1,nboxlist
              j = boxlist(i)
              r0 = boxcost(j)
              IF (r0 .GT. maxcost) THEN
                  maxcost = r0
                  inextboxlist = i
                  inext = j
              ENDIF
          ENDDO
          IF (ppm_dim .GT. 2) THEN
              mins(1)   = min_box(1,inext)
              mins(2)   = min_box(2,inext)
              mins(3)   = min_box(3,inext)
              maxs(1)   = max_box(1,inext)
              maxs(2)   = max_box(2,inext)
              maxs(3)   = max_box(3,inext)
              IF (have_mesh) THEN
                  thisNm(1) = Nm_box(1,inext)
                  thisNm(2) = Nm_box(2,inext)
                  thisNm(3) = Nm_box(3,inext)
              ENDIF
          ELSE
              mins(1)   = min_box(1,inext)
              mins(2)   = min_box(2,inext)
              maxs(1)   = max_box(1,inext)
              maxs(2)   = max_box(2,inext)
              IF (have_mesh) THEN
                  thisNm(1) = Nm_box(1,inext)
                  thisNm(2) = Nm_box(2,inext)
              ENDIF
          ENDIF

          !---------------------------------------------------------------------
          !  Determine best cut direction(s)
          !---------------------------------------------------------------------
          IF (PRESENT(pcost)) THEN 
              CALL ppm_tree_cutdir(xp,Np,weights(:,1),min_box,max_box, &
     &                             inext,ncut,fixed,minboxsize,icut,info,pcost)
          ELSE 
              CALL ppm_tree_cutdir(xp,Np,weights(:,1),min_box,max_box, &
     &                             inext,ncut,fixed,minboxsize,icut,info)
          ENDIF
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Determine best cut position(s)
          !---------------------------------------------------------------------
          IF (PRESENT(pcost)) THEN 
              CALL ppm_tree_cutpos(xp,Np,weights(:,2),min_box,max_box, &
     &                             inext,ncut,minboxsize,icut,cpos,info,pcost)
          ELSE 
              CALL ppm_tree_cutpos(xp,Np,weights(:,2),min_box,max_box, &
     &                             inext,ncut,minboxsize,icut,cpos,info)
          ENDIF
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Align positions with mesh planes if needed
          !---------------------------------------------------------------------
          IF (have_mesh) THEN
              DO i=1,ncut
                  j  = icut(i)
                  r0 = (cpos(i)-mins(j))*meshdxinv(j)
                  k  = NINT(r0)
                  r1 = REAL(k,MK)*meshdx(j)
                  cpos(i) = mins(j) + r1
                  up = .TRUE.
                  IF ((r1-r0) .LT. 0.0_MK) up = .FALSE.
                  !-------------------------------------------------------------
                  !  Check if minboxsizes are respected
                  !-------------------------------------------------------------
                  IF (((cpos(i)-mins(j)) .LT. minboxsize(j)) .OR.   &
     &                ((maxs(j)-cpos(i)) .LT. minboxsize(j))) THEN
                      IF (up) THEN
                          !-----------------------------------------------------
                          !  If we moved up, try down now
                          !-----------------------------------------------------
                          k = k - 1
                          cpos(i) = mins(j) + (REAL(k,MK)*meshdx(j))
                      ELSE
                          !-----------------------------------------------------
                          !  If we moved down, try up
                          !-----------------------------------------------------
                          k = k + 1
                          cpos(i) = mins(j) + (REAL(k,MK)*meshdx(j))
                      ENDIF
                      !---------------------------------------------------------
                      !  Check if minboxsizes are respected now
                      !---------------------------------------------------------
                      IF (((cpos(i)-mins(j)) .LT. minboxsize(j)) .OR.   &
     &                    ((maxs(j)-cpos(i)) .LT. minboxsize(j))) THEN
                          !-----------------------------------------------------
                          !  Cannot subdivide this box along grid lines.
                          !  Remove it from the list of divisible boxes and
                          !  loop.
                          !-----------------------------------------------------
                          DO l=inextboxlist,nboxlist-1
                              boxlist(l) = boxlist(l+1)
                          ENDDO
                          nboxlist = nboxlist - 1
                          GOTO 100
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Subdivide this box
          !---------------------------------------------------------------------
          CALL ppm_tree_boxcut(xp,inext,mins,maxs,ncut,icut,cpos,minc,maxc,  &
     &        lhbx_cut,lpdx_cut,info)
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Update the Nm of the sub-boxes. This needs to be done here for
          !  all sub-boxes since tree_boxcost needs it. 
          !---------------------------------------------------------------------
          IF (have_mesh) THEN
              DO i=1,nbpd
                  Nmc(1,i) = NINT((maxc(1,i)-minc(1,i))*meshdxinv(1))+1
                  Nmc(2,i) = NINT((maxc(2,i)-minc(2,i))*meshdxinv(2))+1
                  IF (ppm_dim .GT. 2) THEN
                      Nmc(3,i) = NINT((maxc(3,i)-minc(3,i))*meshdxinv(3))+1
                  ENDIF
              ENDDO
          ENDIF

          !---------------------------------------------------------------------
          !  Update the costs of the new boxes.
          !  This also grows boxcost.
          !---------------------------------------------------------------------
          IF (PRESENT(pcost)) THEN 
              CALL ppm_tree_boxcost(Nmc,weights(:,1),minc,maxc,   &
     &            nbpd,lhbx_cut,lpdx_cut,costc,info,pcost)
          ELSE
              CALL ppm_tree_boxcost(Nmc,weights(:,1),minc,maxc,   &
     &            nbpd,lhbx_cut,lpdx_cut,costc,info,pcost)
          ENDIF

          !---------------------------------------------------------------------
          !  Add the new boxes to the tree
          !---------------------------------------------------------------------
          IF (have_particles) istart = tree_lhbx(1,inext)
          nadd = 0
          DO i=1,nbpd
              !-----------------------------------------------------------------
              !  If pruneboxes is set, only add boxes of non-zero cost.
              !-----------------------------------------------------------------
              IF ((pruneboxes .AND. ABS(boxcost(i)) .GT. lmyeps) .OR.   &
     &            (.NOT. pruneboxes)) THEN
                  nbox = nbox + 1
                  k    = blevel(inext) + 1
                  IF (k .GT. nlevel) THEN
                      nlevel = k
                      up = .TRUE.
                  ELSE
                      up = .FALSE.
                  ENDIF

                  !-------------------------------------------------------------
                  !  Grow the lists if needed
                  !-------------------------------------------------------------
                  iopt = ppm_param_alloc_grow_preserve
#if   __TYPE == __TREE
                  IF (nbox .GT. nboxalloc .OR. nlevel .GT. nlevelalloc) THEN
                      IF(nbox.GT.nboxalloc) nboxalloc=nboxalloc+   &
     &                    (nbpd**(nlevel-1))
                      IF(nlevel.GT.nlevelalloc) nlevelalloc=nlevelalloc+1
                      IF (ppm_debug .GT. 0) THEN
                          WRITE(mesg,'(A,I3,A,I6,A)') 'Reallocating to ',   &
     &                        nlevelalloc,' levels and ',nboxalloc,' boxes.'
                          CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
                      ENDIF
                      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,nlevelalloc,  &
     &                    min_box,max_box,boxcost,parent,nchld,child,&
     &                    blevel,nbpl,info)
                  ENDIF
#elif __TYPE == __DECOMP
                  IF (nbox .GT. nboxalloc) THEN
                      nboxalloc = nboxalloc + (nbpd**(nlevel-1))
                      IF (ppm_debug .GT. 0) THEN
                          WRITE(mesg,'(A,I3,A)') 'Reallocating to ',   &
     &                        nboxalloc,' boxes.'
                          CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
                      ENDIF
                      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,min_box,max_box,&
     &                    boxcost,nchld,blevel,info)
                  ENDIF
#endif
                  IF (info .NE. ppm_param_success) GOTO 9999

                  !-------------------------------------------------------------
                  !  Store the new boxes
                  !-------------------------------------------------------------
                  IF (ppm_dim .GT. 2) THEN
                      min_box(1,nbox)   = minc(1,i)
                      min_box(2,nbox)   = minc(2,i)
                      min_box(3,nbox)   = minc(3,i)
                      max_box(1,nbox)   = maxc(1,i)
                      max_box(2,nbox)   = maxc(2,i)
                      max_box(3,nbox)   = maxc(3,i)
                      IF (have_mesh) THEN
                          Nm_box(1,nbox)= Nmc(1,i)
                          Nm_box(2,nbox)= Nmc(2,i)
                          Nm_box(3,nbox)= Nmc(3,i)
                      ENDIF
                  ELSE
                      min_box(1,nbox)   = minc(1,i)
                      min_box(2,nbox)   = minc(2,i)
                      max_box(1,nbox)   = maxc(1,i)
                      max_box(2,nbox)   = maxc(2,i)
                      IF (have_mesh) THEN
                          Nm_box(1,nbox)= Nmc(1,i)
                          Nm_box(2,nbox)= Nmc(2,i)
                      ENDIF
                  ENDIF
                  boxcost(nbox)             = costc(i)
                  nchld(nbox)               = 0
                  nchld(inext)              = nchld(inext) + 1
                  blevel(nbox)              = k
#if   __TYPE == __TREE
                  child(1:nbpd,nbox)        = ppm_param_undefined
                  child(nchld(inext),inext) = nbox
                  parent(nbox)              = inext
                  IF (up) nbpl(k)           = 0
                  nbpl(k)                   = nbpl(k) + 1
#endif
                  !-------------------------------------------------------------
                  !  Update the particle index lists
                  !-------------------------------------------------------------
                  IF (have_particles) THEN
                      bpc                   = lhbx_cut(i+1)-lhbx_cut(i)
                      iend                  = istart+bpc-1
                      tree_lpdx(istart:iend)= lpdx_cut(lhbx_cut(i):   &
     &                                                (lhbx_cut(i+1)-1))
                      tree_lhbx(1,nbox)     = istart
                      tree_lhbx(2,nbox)     = iend
                      istart                = iend+1
                  ENDIF
                  nadd                      = nadd + 1
              ENDIF
          ENDDO

          !---------------------------------------------------------------------
          !  Update the list of divisible boxes
          !---------------------------------------------------------------------
          IF (nadd .GT. 0) THEN
              !-----------------------------------------------------------------
              !  Only do this if we added new boxes to the tree
              !-----------------------------------------------------------------
              ibox = nbox - nadd + 1
              CALL ppm_tree_divcheck(min_box(1:ppm_dim,ibox:nbox),    &
     &            max_box(1:ppm_dim,ibox:nbox),nbpd,minboxsize,fixed, &
     &            boxcost(ibox:nbox),ndiv,info)
              IF (info .NE. 0) GOTO 9999
              k  = 0    ! number of added boxes
              k2 = 0    ! number of boxes of non-zero cost
              DO i=1,nadd
                  IF (boxcost(ibox+i-1) .GT. lmyeps) k2 = k2 + 1
                  IF (ndiv(i) .GE. ncut) THEN
                      !---------------------------------------------------------
                      !  If yes, add them to the list of divisible boxes
                      !---------------------------------------------------------
                      IF (k .EQ. 0) THEN
                          !-----------------------------------------------------
                          !  First box to add replaces its parent
                          !-----------------------------------------------------
                          j = inextboxlist
                      ELSE
                          !-----------------------------------------------------
                          !  Following ones are appended at the end of the list
                          !-----------------------------------------------------
                          nboxlist = nboxlist + 1
                          j = nboxlist
                      ENDIF
                      IF (nboxlist .GT. nboxlistalloc) THEN
                          iopt   = ppm_param_alloc_grow_preserve
                          nboxlistalloc = nboxlistalloc + (nbpd**(nlevel-1))
                          ldc(1) = nboxlistalloc
                          CALL ppm_alloc(boxlist,ldc,iopt,info)
                          IF (info.NE.0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_tree',  &
     &                            'list of divisible boxes BOXLIST',    &
     &                            __LINE__,info)
                              GOTO 9999
                          ENDIF 
                      ENDIF 
                      boxlist(j) = ibox+i-1
                      k = k + 1
                  ENDIF
              ENDDO

              !-----------------------------------------------------------------
              !  Compress list if no replacement of parent could be made
              !-----------------------------------------------------------------
              IF (k .EQ. 0) THEN
                  DO i=inextboxlist,nboxlist-1
                      boxlist(i) = boxlist(i+1)
                  ENDDO
                  nboxlist = nboxlist - 1
              ENDIF

              !-----------------------------------------------------------------
              !  Update the number of childless boxes of non-empty cost
              !-----------------------------------------------------------------
              nsubs = nsubs - 1 + k2
          ENDIF           ! nadd.GT.0

          !---------------------------------------------------------------------
          !  Determine if tree is finished
          !---------------------------------------------------------------------
 100      CALL ppm_tree_done(minboxes,nsubs,boxcost,boxlist,nboxlist,   &
     &        nlevel,maxvariance,maxboxcost,mxlev,lcontinue,info)
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Debug output
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I8)') 'Completed iteration ',lctr
              CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
              WRITE(mesg,'(A,I8)') 'Total number of boxes: ',nbox
              CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
              WRITE(mesg,'(A,I8)') 'Number of further divisible boxes: ', &
     &            nboxlist
              CALL ppm_write(ppm_rank,'ppm_tree',mesg,info)
          ENDIF
      ENDDO             ! while lcontinue

#if   __TYPE == __TREE
      !-------------------------------------------------------------------------
      !  Return particle lists
      !-------------------------------------------------------------------------
      NULLIFY(lhbx)
      NULLIFY(lpdx)
      IF (have_particles) THEN
          lhbx => tree_lhbx
          lpdx => tree_lpdx
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Free memory
      !-------------------------------------------------------------------------
 9999 CONTINUE
      iopt = ppm_param_dealloc
      CALL ppm_alloc(boxlist,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'list of divisible boxes BOXLIST',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ndiv,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'number of divisible dimensions NDIV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(icut,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'cut directions ICUT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(cpos,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'cut positions CPOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(minc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'minimum positions of newly cut boxes MINC',__LINE__,info)
      ENDIF
      CALL ppm_alloc(maxc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'maximum positions of newly cut boxes MAXC',__LINE__,info)
      ENDIF
      CALL ppm_alloc(costc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'costs of new boxes COSTC',__LINE__,info)
      ENDIF
      IF (have_mesh) THEN
          CALL ppm_alloc(Nm_box,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'number of mesh points per box NM_BOX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(Nmc,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'number of mesh points of newly cut boxes NMC',__LINE__,info)
          ENDIF
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcst_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcst_d,ldc,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'particle cost part PCOST',__LINE__,info)
      ENDIF
#ifdef __MPI
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcsum_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcsum_d,ldc,iopt,info)
#endif
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'particle cost sums PCSUM',__LINE__,info)
      ENDIF
#endif
      IF (have_particles) THEN
          CALL ppm_alloc(lpdx_cut,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'particle index pointers LPDX_CUT',__LINE__,info)
          ENDIF
          CALL ppm_alloc(lhbx_cut,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'pointers to first particle in box LHBX_CUT',__LINE__,info)
          ENDIF
      ENDIF
#if   __TYPE == __DECOMP
      CALL ppm_alloc(boxcost,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'costs of all boxes BOXCOST',__LINE__,info)
      ENDIF
      CALL ppm_alloc(blevel,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &        'tree levels of boxes BLEVEL',__LINE__,info)
      ENDIF
      IF (have_particles) THEN
          CALL ppm_alloc(npbx,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'number of particles per box NPBX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(cbox,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'temporary box pointers CBOX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(tree_lhbx,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'list of divisible boxes BOXLIST',__LINE__,info)
          ENDIF
          CALL ppm_alloc(tree_lpdx,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_tree',     &
     &            'list of divisible boxes BOXLIST',__LINE__,info)
          ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 8000 CALL substop('ppm_tree',t0,info)
      RETURN
#if   __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_dd
#endif
#elif __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_ts
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_td
#endif
#endif
