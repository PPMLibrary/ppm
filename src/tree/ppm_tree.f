      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_tree
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

#if   __TYPE == __DECOMP
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_ds(xp,Np,Nm,min_dom,max_dom,treetype,         &
      &          minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost, &
      &          fixed,weights,min_box,max_box,nbox,nchld,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_dd(xp,Np,Nm,min_dom,max_dom,treetype,         &
      &          minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost, &
      &          fixed,weights,min_box,max_box,nbox,nchld,info,pcost)
#endif
#elif __TYPE == __TREE
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_ts(xp,Np,Nm,min_dom,max_dom,treetype,             &
      &          minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,     &
      &          maxlevels,fixed,weights,min_box,max_box,lhbx,lpdx,boxcost, &
      &          parent,nchld,child,blevel,nbox,nbpl,nlevel,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_td(xp,Np,Nm,min_dom,max_dom,treetype,             &
      &          minboxes,pruneboxes,minboxsize,maxvariance,maxboxcost,     &
      &          maxlevels,fixed,weights,min_box,max_box,lhbx,lpdx,boxcost, &
      &          parent,nchld,child,blevel,nbox,nbpl,nlevel,info,pcost)
#endif
#endif
      !!! This routine performs a generic tree decomposition
      !!! using recursive orthogonal multisection (ROM) in either
      !!! 2, 4, or 8 subunits (binary tree, quadtree,
      !!! octtree). The tree can be based on scattered data
      !!! points (if Np > 0), on regular mesh points (if
      !!! Nm > 0) or purely on the geometry of the domain
      !!! (if both Np and Nm are 0).
#if __TYPE == __TREE
      !!!
      !!! In this version the full tree information (connectivity,
      !!! levels, etc.) is computed and returned.
      !!! If __TYPE == __DECOMP, only the space decomposition is
      !!! done, but no tree information returned.
#endif
      !!!
      !!! [NOTE]
      !!! .Remarks
      !!! ==============================================================
      !!! If Np > 0, the particle-caused cost per sub is
      !!! given by number(particles) [or sum(pcost)] in that sub.
      !!! If Nm(1:ppm_dim) > 1, the mesh-caused cost is
      !!! given by the number of mesh points in each sub.
      !!! The geometry-based cost is always given by the
      !!! volume of each sub.
      !!!
      !!! These costs can be freely composed using the
      !!! weights argument.
      !!!
      !!! ppm_decomp_cartesian is not made obsolete by this
      !!! routine since the latter can gracefully produce ANY
      !!! (even prime) number of equisized subdomains,
      !!! whereas this routine is restricted to powers of
      !!! 2, 4, or 8 (for equisized subs!).
      !!!
      !!! directly compute box costs in tree_cutpos
      !!! (since we do the allreduce there anyway) and return
      !!! them. only recompute (with tree_boxcost) if the cut
      !!! planes were shifted.
      !!! ==============================================================

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
      USE ppm_module_decomp
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! The data points
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_dom
      !!! Minimum coordinate of the domain
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: max_dom
      !!! Maximum coordinate of the domain
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: minboxsize
      !!! Miminum box size in all directions
      REAL(MK), DIMENSION(3,2), INTENT(IN   ) :: weights
      !!! Weights for the three cost contributions (particles,
      !!! mesh points, volume) (1st index) for box cost (weights(:,1))
      !!! and the determination of the cut planes (weights(:,2)).
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Np, specifying the
      !!! cost of each data point.
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      !!! Flag which tells for each spatial dimension (1..ppm_dim) if it is
      !!! fixed (i.e. no cuts perpendicular to it are allowed). fixed=(F,F,T)
      !!! will e.g. enforce z-pencil decompositions.
      LOGICAL                 , INTENT(IN   ) :: pruneboxes
      !!! `TRUE` to prune the tree to only contain boxes of non-zero cost.
      !!! `FALSE` to get a tree with all boxes (also the empty ones).
      REAL(MK), DIMENSION(:,:), POINTER       :: min_box
      !!! Min. extents of the boxes
      REAL(MK), DIMENSION(:,:), POINTER       :: max_box
      !!! Max. extents of the boxes
      INTEGER                 , INTENT(IN   ) :: Np
      !!! Number of data points.
      !!! If <= 0, decomposition is based on geometry and mesh only.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Number of grid points in the global mesh. (0,0,0) if there is
      !!! no mesh. If a mesh is present, the box boundaries will be aligned
      !!! with mesh planes.
      INTEGER                 , INTENT(IN   ) :: treetype
      !!! Type of multisection tree. One of:
      !!!
      !!! *  ppm_param_tree_bin
      !!! *  ppm_param_tree_quad
      !!! *  ppm_param_tree_oct (3D only)
      !!!
      !!! For binary, quad- or oct-tree.
      INTEGER                 , INTENT(IN   ) :: minboxes
      !!! Minimum number of childless boxes (leaves) of non-zero cost to be
      !!! created. Set this to -1 if there is no minimum requirement.
      REAL(MK)                , INTENT(IN   ) :: maxvariance
      !!! Maximum variance of cost allowed between boxes. The tree stops as
      !!! soon as the variance of costs of all boxes is below this max.
      !!! Set this to -1 to disable this criterion.
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      !!! Maximum cost per box. Subdivision will stop
      !!! as soon as all boxes have costs below this value. Set this to -1
      !!! to not impose any limit.
      INTEGER                 , INTENT(  OUT) :: nbox
      !!! The total number of boxes
      INTEGER , DIMENSION(:  ), POINTER       :: nchld
      !!! Number of children of each box.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 upon success

#if   __TYPE == __TREE
      INTEGER                 , INTENT(IN   ) :: maxlevels
      !!! Maximum number of levels to create. Tree stops as soon as
      !!! this is reached. The root box is counted as level 1. Set to <= 0
      !!! for unlimited levels. This input is only present in the
      !!! __TREE version, not in the __DECOMP version.
      REAL(MK), DIMENSION(:  ), POINTER       :: boxcost
      !!! Costs of all boxes 1..nbox.
      INTEGER , DIMENSION(:  ), POINTER       :: parent
      !!! Index of the parent box of each box. ppm_param_undefined if no
      !!! parent (i.e. root box)
      INTEGER , DIMENSION(:  ), POINTER       :: nbpl
      !!! the number of boxes per level. Level 1 is the root box.
      INTEGER , DIMENSION(:  ), POINTER       :: blevel
      !!! tree level of each box. 1..nbox. Level 1 is the root box.
      INTEGER , DIMENSION(:,:), POINTER       :: lhbx
      !!! pointer to first (1,:) and last (2,:) data point (in
      !!! lpdx) in each tree box. This is only present if
      !!! Np > 0. Only points on the local processor are considered.
      !!! Entries for non-leaf boxes are only true if `pruneboxes=FALSE`
      INTEGER , DIMENSION(:  ), POINTER       :: lpdx
      !!! Pointers to the data points (in xp) in each tree box. Box ib
      !!! conatins points lpdx(lhbx(1,ib):(lhbx(2,ib)))
      !!! This is only present if Np > 0. Only points on the local
      !!! processor are considered. Entries for non-leaf boxes are
      !!! only true if `pruneboxes=FALSE`
      INTEGER , DIMENSION(:,:), POINTER       :: child
      !!! Indices of all children of a box.
      !!!
      !!! 1st index: child ID                                                  +
      !!! 2nd: box ID.
      INTEGER                 , INTENT(  OUT) :: nlevel
      !!! The number of levels. Level 1 is the root box.
#endif
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(1:ppm_dim)    :: mins,maxs
      REAL(MK), DIMENSION(1:ppm_dim)    :: meshdx,meshdxinv
      REAL(MK), DIMENSION(:  ), POINTER :: cpos  => NULL()
      REAL(MK), DIMENSION(:  ), POINTER :: costc => NULL()

      REAL(MK), DIMENSION(:,:), POINTER :: minc  => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: maxc  => NULL()
      REAL(MK)                          :: t0,lmyeps,r0,r1,maxcost
#if   __TYPE == __DECOMP
      REAL(MK), DIMENSION(:  ), POINTER :: boxcost => NULL()

      INTEGER , DIMENSION(:  ), POINTER :: blevel  => NULL()
      INTEGER                           :: nlevel
#endif
      INTEGER , DIMENSION(:  ), POINTER :: icut  => NULL()
      INTEGER , DIMENSION(ppm_dim)      :: thisNm
      INTEGER , DIMENSION(2*ppm_dim)    :: ghostNm
      INTEGER , DIMENSION(2)            :: ldc
      INTEGER                           :: nboxlist,nadd,k2,itype
      INTEGER                           :: nboxalloc,nlevelalloc
      INTEGER                           :: i,j,k,l,iopt,inext,ncut,nbpd
      INTEGER                           :: ibox,nsubs
      INTEGER                           :: inextboxlist,lctr
      INTEGER                           :: nboxlistalloc
      INTEGER                           :: info2,mxlev,bpc,istart,iend

      CHARACTER(LEN=ppm_char) :: mesg
      CHARACTER(LEN=ppm_char) :: caller = 'ppm_tree'

      LOGICAL :: up,nofixed,simpleweights
      LOGICAL :: lcontinue

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)
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
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check what kind of input data is given
      !-------------------------------------------------------------------------
      have_particles = .FALSE.
      have_mesh      = .FALSE.
      IF (Np .GT. 0) have_particles = .TRUE.
      IF (SIZE(Nm,1) .GE. ppm_dim) THEN
         IF (ppm_dim .GT. 2) THEN
            IF ((Nm(1).GT.1).AND.(Nm(2).GT.1).AND.(Nm(3).GT.1)) &
            &  have_mesh = .TRUE.
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

     ! IF ((itype .EQ. ppm_param_tree_quad) .AND. (.NOT.have_mesh) .AND.  &
     !&    (have_particles) .AND. (.NOT.pruneboxes) .AND. (nofixed) .AND. &
     !&    (simpleweights) .AND. (ppm_dim .EQ. 2)) THEN
     !   IF (PRESENT(pcost)) THEN
     !       CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&           minboxsize(2)),maxvariance,min_box,max_box,nbox,info,pcost)
     !   ELSE
     !       CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&           minboxsize(2)),maxvariance,min_box,max_box,nbox,info)
     !   ENDIF
     !   GOTO 8000
     ! ENDIF
     ! IF ((itype .EQ. ppm_param_tree_oct) .AND. (.NOT.have_mesh) .AND.   &
     !&    (have_particles) .AND. (.NOT.pruneboxes) .AND. (nofixed) .AND. &
     !&    (simpleweights) .AND. (ppm_dim .EQ. 3)) THEN
     !   IF (PRESENT(pcost)) THEN
     !       CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&           minboxsize(2),minboxsize(3)),maxvariance,min_box,max_box,   &
     !&           nbox,info,pcost)
     !   ELSE
     !       CALL ppm_decomp_tree(xp,Np,min_dom,max_dom,MAX(minboxsize(1), &
     !&           minboxsize(2),minboxsize(3)),maxvariance,min_box,max_box,   &
     !&           nbox,info)
     !   ENDIF
     !   GOTO 8000
     ! ENDIF

      !-------------------------------------------------------------------------
      !  Store the mesh spacings (if needed)
      !-------------------------------------------------------------------------
      IF (have_mesh) THEN
         meshdx(1:ppm_dim) = (max_dom(1:ppm_dim) - &
         &                    min_dom(1:ppm_dim))/ &
         &                    REAL(Nm(1:ppm_dim)-1,MK)
         !check for round-off problems and fix them if necessary
         DO k=1,ppm_dim
            DO WHILE (min_dom(k)+(Nm(k)-1)*meshdx(k).LT.max_dom(k))
               meshdx(k)=meshdx(k)+EPSILON(meshdx(k))
            ENDDO
         ENDDO
         check_true(<#ALL(min_dom(1:ppm_dim)+(Nm(1:ppm_dim)-1)*meshdx(1:ppm_dim).GE.max_dom(1:ppm_dim))#>,"round-off problem in mesh creation")

         meshdxinv(1:ppm_dim) = 1.0_MK/meshdx(1:ppm_dim)
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the number of boxes and cuts per subdivision
      !-------------------------------------------------------------------------
      SELECT CASE (itype)
      CASE (ppm_param_tree_bin)
          nbpd = 2
          ncut = 1
          IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,caller,'Creating binary tree.',info)
          ENDIF
      CASE (ppm_param_tree_quad)
          nbpd = 4
          ncut = 2
          IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,caller,'Creating quad-tree.',info)
          ENDIF
      CASE (ppm_param_tree_oct)
          nbpd = 8
          ncut = 3
          IF (ppm_debug .GT. 0) THEN
             CALL ppm_write(ppm_rank,caller,'Creating oct-tree.',info)
          ENDIF
      CASE DEFAULT
          info = ppm_error_error
          CALL ppm_error(ppm_err_argument,caller,     &
          &   'unknown tree type specified !',__LINE__,info)
          GOTO 9999
      END SELECT

      !-------------------------------------------------------------------------
      !  Clear module pointers
      !-------------------------------------------------------------------------
      NULLIFY(tree_lhbx)
      NULLIFY(tree_lpdx)

      !-------------------------------------------------------------------------
      !  Allocate tree data structures
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      !-------------------------------------------------------------------------
      !  Guess tree size based on assumed uniform distribution
      !-------------------------------------------------------------------------
      nbox   = 1
      nlevel = 1
      IF (Np .GT. 0) THEN
          IF (maxboxcost .GT. 0.0_MK) THEN
              ! the number of levels needed under uniform particle
              ! distribution
              nlevelalloc = CEILING(LOG(REAL(Np,MK)/maxboxcost)/   &
              &             LOG(REAL(nbpd,MK)))
#if __TYPE ==  __TREE
          ELSEIF (maxlevels .GT. 0) THEN
              ! assume we hit maxlevels
              nlevelalloc = maxlevels
#endif
          ELSE
              ! default assumption
              nlevelalloc = 3
              ! nlevelalloc = MAX(NINT(LOG(REAL(Np,MK))/LOG(2.0_MK)),3)
              ! JHW 20061108
          ENDIF
      ELSE
          ! Assume 3 levels if no particles are present
          nlevelalloc = 3
          ! nlevelalloc = MAX(NINT(LOG(REAL(Np,MK))/LOG(2.0_MK)),3)
          ! JHW 20061108
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
          CALL ppm_write(ppm_rank,caller,mesg,info)
      ENDIF
#if   __TYPE == __TREE
      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,nlevelalloc,min_box,max_box, &
      &    boxcost,parent,nchld,child,blevel,nbpl,info)
#elif __TYPE == __DECOMP
      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,min_box,max_box, &
      &    boxcost,nchld,blevel,info)
#endif
      IF (info .NE. ppm_param_success) GOTO 9999

      !-------------------------------------------------------------------------
      !  Allocate memory for box cut
      !-------------------------------------------------------------------------
      ldc(1) = ppm_dim
      ldc(2) = 2**ncut
      CALL ppm_alloc(minc,ldc,iopt,info)
      or_fail_alloc('lower coordinates of new boxes MINC')

      CALL ppm_alloc(maxc,ldc,iopt,info)
      or_fail_alloc('upper coordinates of new boxes MAXC')

      !-------------------------------------------------------------------------
      !  Allocate local data structures
      !-------------------------------------------------------------------------
      nboxlist = 1
      nboxlistalloc = 10*nbpd**(nlevelalloc-1)   ! list of leaves. guess...
      ! alocate 10 times more space for leaves... JHW 20061108
      ldc(1) = nboxlistalloc
      CALL ppm_alloc(boxlist,ldc,iopt,info)
      or_fail_alloc('list of divisible boxes BOXLIST')

      boxlist(1) = 1
      IF (have_mesh) THEN
         ldc(1) = ppm_dim
         ldc(2) = nbpd
         CALL ppm_alloc(Nmc,ldc,iopt,info)
         or_fail_alloc('list of divisible boxes BOXLIST')
      ENDIF
      IF (have_particles) THEN
         ldc(1) = 2**ncut
         CALL ppm_alloc(cbox,ldc,iopt,info)
         or_fail_alloc('temporary box pointers CBOX')

         CALL ppm_alloc(npbx,ldc,iopt,info)
         or_fail_alloc('number of particles per box NPBX')
      ENDIF

      !-------------------------------------------------------------------------
      !  The domain itself is the root box. Get the tree started!
      !-------------------------------------------------------------------------
      IF (ppm_dim .GT. 2) THEN
         min_box(1:3,1) = min_dom(1:3)
         max_box(1:3,1) = max_dom(1:3)
         IF (have_mesh) THEN
            Nm_box(1:3,1) = Nm(1:3)
         ENDIF
      ELSE
         min_box(1:2,1) = min_dom(1:2)
         max_box(1:2,1) = max_dom(1:2)
         IF (have_mesh) THEN
            Nm_box(1:2,1) = Nm(1:2)
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
         thisNm  = 1
         ghostNm = 0
         IF (ppm_dim .EQ. 2) THEN
            CALL ppm_util_rank2d(xp,Np,min_dom,max_dom,thisNm,ghostNm,     &
            &    tree_lpdx,lhbx_cut,info)
         ELSE
            CALL ppm_util_rank3d(xp,Np,min_dom,max_dom,thisNm,ghostNm,     &
            &    tree_lpdx,lhbx_cut,info)
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
      or_fail_alloc('box costs COSTC')

#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcst_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcst_d,ldc,iopt,info)
#endif
      or_fail_alloc('particle cost part PCOST')

#ifdef __MPI
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcsum_s,ldc,iopt,info)
#else
      CALL ppm_alloc(pcsum_d,ldc,iopt,info)
#endif
      or_fail_alloc('particle cost sums PCSUM')
#endif

      !-------------------------------------------------------------------------
      !  Compute cost of root box
      !-------------------------------------------------------------------------
      IF (PRESENT(pcost)) THEN
         CALL ppm_tree_boxcost(Nm_box,weights(:,1),min_box,max_box,  &
         &    1,lhbx_cut,tree_lpdx,boxcost,info,pcost)
      ELSE
         CALL ppm_tree_boxcost(Nm_box,weights(:,1),min_box,max_box,  &
         &    1,lhbx_cut,tree_lpdx,boxcost,info)
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
         or_fail_alloc('particle list header pointers LHBX_CUT')
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if there is anything to be done at all
      !-------------------------------------------------------------------------
      CALL ppm_tree_done(minboxes,nsubs,boxcost,boxlist,nboxlist,   &
      &    nlevel,maxvariance,maxboxcost,mxlev,lcontinue,info)
      IF (info .NE. ppm_param_success) GOTO 9999
      IF ((.NOT.lcontinue) .AND. (ppm_debug .GT. 0)) THEN
         CALL ppm_write(ppm_rank,caller,     &
         &    'Nothing to be done. Exiting.',info)
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
            CALL ppm_write(ppm_rank,caller,     &
            &    'Initial domain is not divisible. Done.',info)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for cut directions and positions
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldc(1) = ncut
      CALL ppm_alloc(icut,ldc,iopt,info)
      or_fail_alloc('list of cut directions ICUT')

      icut = ppm_param_undefined
      CALL ppm_alloc(cpos,ldc,iopt,info)
      or_fail_alloc('list of cut positions CPOS')

      cpos = -HUGE(cpos(1))

      !-------------------------------------------------------------------------
      !  Subdivide until done
      !-------------------------------------------------------------------------
      lctr = 0
      lcontinue_do: DO WHILE (lcontinue)
          lctr = lctr + 1

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
             & inext,ncut,fixed,minboxsize,icut,info,pcost)
          ELSE
             CALL ppm_tree_cutdir(xp,Np,weights(:,1),min_box,max_box, &
             & inext,ncut,fixed,minboxsize,icut,info)
          ENDIF
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Determine best cut position(s)
          !---------------------------------------------------------------------
          IF (PRESENT(pcost)) THEN
             CALL ppm_tree_cutpos(xp,Np,weights(:,2),min_box,max_box, &
             & inext,ncut,minboxsize,icut,cpos,info,pcost)
          ELSE
             CALL ppm_tree_cutpos(xp,Np,weights(:,2),min_box,max_box, &
             & inext,ncut,minboxsize,icut,cpos,info)
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
                &   ((maxs(j)-cpos(i)) .LT. minboxsize(j))) THEN
                   IF (up) THEN
                      !-----------------------------------------------------
                      !  If we moved up, try down now
                      !-----------------------------------------------------
                      k = k - 1
                      cpos(i) = mins(j) + REAL(k,MK)*meshdx(j)
                   ELSE
                      !-----------------------------------------------------
                      !  If we moved down, try up
                      !-----------------------------------------------------
                      k = k + 1
                      cpos(i) = mins(j) + REAL(k,MK)*meshdx(j)
                   ENDIF
                   !---------------------------------------------------------
                   !  Check if minboxsizes are respected now
                   !---------------------------------------------------------
                   IF (((cpos(i)-mins(j)) .LT. minboxsize(j)) .OR.   &
                   &   ((maxs(j)-cpos(i)) .LT. minboxsize(j))) THEN
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
             ENDDO !i=1,ncut
          ENDIF !(have_mesh)

          !---------------------------------------------------------------------
          !  Subdivide this box
          !---------------------------------------------------------------------
          CALL ppm_tree_boxcut(xp,inext,mins,maxs,ncut,icut,cpos,minc,maxc,  &
          &    lhbx_cut,lpdx_cut,info)
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
             & nbpd,lhbx_cut,lpdx_cut,costc,info,pcost)
          ELSE
             CALL ppm_tree_boxcost(Nmc,weights(:,1),minc,maxc,   &
             & nbpd,lhbx_cut,lpdx_cut,costc,info,pcost)
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
              &   (.NOT. pruneboxes)) THEN
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
                      &  (nbpd**(nlevel-1))
                      IF(nlevel.GT.nlevelalloc) nlevelalloc=nlevelalloc+1
                      IF (ppm_debug .GT. 0) THEN
                          WRITE(mesg,'(A,I3,A,I6,A)') 'Reallocating to ',   &
                          & nlevelalloc,' levels and ',nboxalloc,' boxes.'
                          CALL ppm_write(ppm_rank,caller,mesg,info)
                      ENDIF
                      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,nlevelalloc,  &
                      &    min_box,max_box,boxcost,parent,nchld,child,&
                      &    blevel,nbpl,info)
                  ENDIF
#elif __TYPE == __DECOMP
                  IF (nbox .GT. nboxalloc) THEN
                      nboxalloc = nboxalloc + (nbpd**(nlevel-1))
                      IF (ppm_debug .GT. 0) THEN
                          WRITE(mesg,'(A,I3,A)') 'Reallocating to ',   &
                          & nboxalloc,' boxes.'
                          CALL ppm_write(ppm_rank,caller,mesg,info)
                      ENDIF
                      CALL ppm_tree_alloc(iopt,nboxalloc,nbpd,min_box,max_box,&
                      &    boxcost,nchld,blevel,info)
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
                      &                       (lhbx_cut(i+1)-1))
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
              &    max_box(1:ppm_dim,ibox:nbox),nbpd,minboxsize,fixed, &
              &    boxcost(ibox:nbox),ndiv,info)
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
                          or_fail_alloc('list of divisible boxes BOXLIST')
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
          &    nlevel,maxvariance,maxboxcost,mxlev,lcontinue,info)
          IF (info .NE. ppm_param_success) GOTO 9999

          !---------------------------------------------------------------------
          !  Debug output
          !---------------------------------------------------------------------
          IF (ppm_debug .GT. 1) THEN
              WRITE(mesg,'(A,I8)') 'Completed iteration ',lctr
              CALL ppm_write(ppm_rank,caller,mesg,info)
              WRITE(mesg,'(A,I8)') 'Total number of boxes: ',nbox
              CALL ppm_write(ppm_rank,caller,mesg,info)
              WRITE(mesg,'(A,I8)') 'Number of further divisible boxes: ', &
              & nboxlist
              CALL ppm_write(ppm_rank,caller,mesg,info)
          ENDIF
      ENDDO lcontinue_do ! while lcontinue

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
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'list of divisible boxes BOXLIST',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ndiv,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'number of divisible dimensions NDIV',__LINE__,info)
      ENDIF
      CALL ppm_alloc(icut,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'cut directions ICUT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(cpos,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'cut positions CPOS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(minc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'minimum positions of newly cut boxes MINC',__LINE__,info)
      ENDIF
      CALL ppm_alloc(maxc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'maximum positions of newly cut boxes MAXC',__LINE__,info)
      ENDIF
      CALL ppm_alloc(costc,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'costs of new boxes COSTC',__LINE__,info)
      ENDIF
      IF (have_mesh) THEN
          CALL ppm_alloc(Nm_box,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'number of mesh points per box NM_BOX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(Nmc,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'number of mesh points of newly cut boxes NMC',__LINE__,info)
          ENDIF
      ENDIF
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcst_s,ldc,iopt,info2)
#else
      CALL ppm_alloc(pcst_d,ldc,iopt,info2)
#endif
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'particle cost part PCOST',__LINE__,info)
      ENDIF
#ifdef __MPI
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(pcsum_s,ldc,iopt,info2)
#else
      CALL ppm_alloc(pcsum_d,ldc,iopt,info2)
#endif
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'particle cost sums PCSUM',__LINE__,info)
      ENDIF
#endif
      IF (have_particles) THEN
          CALL ppm_alloc(lpdx_cut,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'particle index pointers LPDX_CUT',__LINE__,info)
          ENDIF
          CALL ppm_alloc(lhbx_cut,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'pointers to first particle in box LHBX_CUT',__LINE__,info)
          ENDIF
      ENDIF
#if   __TYPE == __DECOMP
      CALL ppm_alloc(boxcost,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'costs of all boxes BOXCOST',__LINE__,info)
      ENDIF
      CALL ppm_alloc(blevel,ldc,iopt,info2)
      IF (info2 .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,caller,     &
     &        'tree levels of boxes BLEVEL',__LINE__,info)
      ENDIF
      IF (have_particles) THEN
          CALL ppm_alloc(npbx,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'number of particles per box NPBX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(cbox,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'temporary box pointers CBOX',__LINE__,info)
          ENDIF
          CALL ppm_alloc(tree_lhbx,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'list of divisible boxes BOXLIST',__LINE__,info)
          ENDIF
          CALL ppm_alloc(tree_lpdx,ldc,iopt,info2)
          IF (info2 .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,     &
     &            'list of divisible boxes BOXLIST',__LINE__,info)
          ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 8000 CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (weights(1,1).EQ.0.0_MK.AND.weights(2,1).EQ.0.0_MK.AND.   &
     &      weights(3,1).EQ.0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,    &
     &          'At least one weights(:,1) must be non-zero!',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (weights(1,2).EQ.0.0_MK.AND.weights(2,2).EQ.0.0_MK.AND.   &
     &      weights(3,2).EQ.0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,    &
     &          'At least one weights(:,2) must be non-zero!',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (treetype.EQ.ppm_param_tree_oct.AND.ppm_dim.EQ.2) THEN
            info = ppm_error_warning
            CALL ppm_error(ppm_err_argument,caller,    &
     &          'Octtree is not possible in 2d. Reverting to quadtree.',  &
     &          __LINE__,info)
            itype = ppm_param_tree_quad
         ENDIF
         DO i=1,ppm_dim
            IF (minboxsize(i) .LT. 0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,caller,     &
     &             'the minimum box size must be > 0 !',__LINE__,info)
               GOTO 8888
            ENDIF
            IF (min_dom(i) .GT. max_dom(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,caller,   &
     &             'min_dom must be <= max_dom !',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
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
