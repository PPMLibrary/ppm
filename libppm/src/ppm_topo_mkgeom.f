      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mkgeom
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is the topology routine for purely
      !                 geometry-based decompositions, i.e. without
      !                 particles and without meshes.
      !                     ppm_param_decomp_bisection
      !                     ppm_param_decomp_xpencil 
      !                     ppm_param_decomp_ypencil 
      !                     ppm_param_decomp_zpencil 
      !                     ppm_param_decomp_xy_slab
      !                     ppm_param_decomp_xz_slab
      !                     ppm_param_decomp_yz_slab
      !                     ppm_param_decomp_cuboid
      !                     ppm_param_decomp_user_defined
      !
      !                 In the user_defined case, the user must submit 
      !                 existing subdomains and all of min_sub, max_sub, 
      !                 cost, nsubs must be provided. 
      !                 Subs are then mapped onto processors and
      !                 the topology is stored.
      !
      !  Input        : decomp       (I) type of decomposition (see above)
      !                 assig        (I) the type of subdomain-to-processor
      !                                  assignment. One of:
      !                                     ppm_param_assign_internal
      !                                     ppm_param_assign_nodal_cut
      !                                     ppm_param_assign_nodal_comm
      !                                     ppm_param_assign_dual_cut
      !                                     ppm_param_assign_dual_comm
      !                                     ppm_param_assign_user_defined
      !                                  The latter uses the external
      !                                  library METIS and is only
      !                                  available if ppm was compiled with
      !                                  METIS support.
      !                 min_phys(:)  (F) the minimum coordinate of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinate of the 
      !                                  physical/computational domain 
      !                 bcdef(:)     (I) the definition of the BC 
      !                 ghostsize    (F) the size (width) of the GL 
      !
      !  Input/output : topo_id      (I) topology identifier (user
      !                                  numbering) for which to create a
      !                                  mesh. Can be an existing one (in
      !                                  which case it is REPLACED), a new
      !                                  one or <0 (in which case ppm
      !                                  internally generates the next
      !                                  available topology ID). Topology 0
      !                                  is reserved for the ring topology
      !                                  (null decomposition) and may only
      !                                  be used with particles.
      !                                  On return: actually used or
      !                                  created topology ID.
      !                 min_sub(:,:) (F) the min. extent of the subdomains.
      !                                  (either user-specified on input or
      !                                  decomposition result on output)
      !                 max_sub(:,:) (F) the max. extent of the subdomains
      !                 sub2proc(:)  (I) processor affiliation of each
      !                                  subdomain. User-defined on
      !                                  input or library result on
      !                                  output.
      !                 nsubs        (I) the total number of subdomains (on
      !                                  all processors in total). Either
      !                                  user-specified on input or
      !                                  decomposoiton result on output.
      !                 cost(:)      (F) estimated cost associated with
      !                                  subdomains. Either user-specified
      !                                  on input or decomp. result. The
      !                                  cost of a subdomain is given by
      !                                  its volume.
      !
      !  Output       : isublist(:)  (I) list of subdomains handled by the
      !                                  local processor
      !                 nsublist     (I) number of subdomains handled by
      !                                  the local processor
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      : In the case of user-defined subs, a lot is
      !                 currently trusted. We may want to include further
      !                 checks in that case.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_mkgeom.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.11  2006/09/04 18:34:56  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.10  2006/02/03 09:34:01  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.9  2005/02/01 16:39:27  ivos
      !  bugfix: fixed wrong symbolic constant in argument check.
      !
      !  Revision 1.8  2005/02/01 15:29:43  ivos
      !  bugfix: info test in assign_user_defined was wrong.
      !
      !  Revision 1.7  2005/01/16 14:44:18  michaebe
      !  added use ppm_module_alloc
      !
      !  Revision 1.6  2005/01/13 13:20:05  ivos
      !  Added possiblity of user-defined sub2proc assignment.
      !
      !  Revision 1.5  2004/12/03 17:18:50  ivos
      !  Changed to use 2d array weights since ppm_tree has changed.
      !
      !  Revision 1.4  2004/11/30 15:05:20  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.3  2004/11/05 14:15:17  ivos
      !  Commented file output of subs.
      !
      !  Revision 1.2  2004/10/05 09:04:13  ivos
      !  security fix: topoid=0 is reserved for the ring. ppm_topo_store relies
      !  on this when determining the internal topoid. Added check that topoid 0
      !  is not used for non-ring topologies. If it is, it will be reset to -1
      !  and a new topology will be created. Updated the comment header
      !  accordingly.
      !
      !  Revision 1.1  2004/09/24 14:58:53  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_mkgeom_s(decomp,assig,min_phys,max_phys,bcdef,  &
     &    ghostsize,topo_id,min_sub,max_sub,cost,sub2proc,nsubs,isublist, &
     &    nsublist,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkgeom_d(decomp,assig,min_phys,max_phys,bcdef,  &
     &    ghostsize,topo_id,min_sub,max_sub,cost,sub2proc,nsubs,isublist, &
     &    nsublist,info)
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
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_define_subs_bc
      USE ppm_module_topo_subs2proc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_find_neigh
      USE ppm_module_tree
      USE ppm_module_alloc
      USE ppm_module_topo_box2subs
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: assig
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      INTEGER                 , INTENT(INOUT) :: nsubs,topo_id,decomp
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                 , INTENT(  OUT) :: nsublist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                           :: i,topoid,j,treetype,nbox,isub,iopt
      INTEGER, DIMENSION(:,:), POINTER  :: ineigh,subs_bc
      INTEGER, DIMENSION(1  )           :: Nmdummy,ldc
      INTEGER, DIMENSION(3,1)           :: nnodes
      INTEGER, DIMENSION(:  ), POINTER  :: nneigh,nchld
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(1,1)          :: xpdummy
      LOGICAL , DIMENSION(ppm_dim)      :: fixed
      REAL(MK), DIMENSION(3,2)          :: weights
      REAL(MK), DIMENSION(:,:), POINTER :: min_box,max_box
      CHARACTER(LEN=ppm_char)           :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_mkgeom',t0,info)
#if    __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif  __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF(ppm_debug.GT.0) THEN
         IF (.NOT. ppm_initialized) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mkgeom',  &
     &           'Please call ppm_init first!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF(ghostsize .LT. 0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &          'ghostsize must be >= 0.0',__LINE__, info)
            GOTO 9999
         ENDIF
         DO i=1,ppm_dim
            IF(max_phys(i).LE.min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'max_phys must be > min_phys',__LINE__, info)
               GOTO 9999
            ENDIF
         ENDDO
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (decomp .NE. ppm_param_decomp_user_defined) THEN
               info = ppm_error_warning
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'decomp type is set to user_defined for this assignment',&
     &             __LINE__, info)
               decomp = ppm_param_decomp_user_defined
            ENDIF
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 9999
            ENDIF
            DO i=1,nsubs
                IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 9999
                ENDIF
            ENDDO
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
            IF ((.NOT.ASSOCIATED(min_sub)).OR.(.NOT.ASSOCIATED(max_sub))) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &              'min_sub/max_sub must be allocated for user def. decomp',&
     &              __LINE__, info)
                GOTO 9999
            ENDIF
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'no subs defined in user_defined decomposition',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            !-------------------------------------------------------------------
            !  Check that the user-defined subs add up to the whole
            !  computational domain.
            !  ONE COULD DO MORE TESTS HERE.
            !-------------------------------------------------------------------
            parea = (max_phys(1)-min_phys(1))*(max_phys(2)-min_phys(2))
            IF(ppm_dim.EQ.3) THEN
               parea = parea*(max_phys(3)-min_phys(3))
            ENDIF
            sarea = 0.0_MK
            DO i=1,nsubs
               larea = (max_sub(1,i)-min_sub(1,i))*(max_sub(2,i)-min_sub(2,i))
               IF(ppm_dim.EQ.3) THEN
                  larea = larea * (max_sub(3,i)-min_sub(3,i))
               END IF
               sarea = sarea + larea
            ENDDO
            IF(ABS(sarea-parea)/parea.GE.lmyeps) THEN
               !----------------------------------------------------------------
               !  Mismatch!
               !----------------------------------------------------------------
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
     &              'faulty subdomains defined',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF   
      ENDIF

      IF (topo_id .EQ. 0) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_argument, 'ppm_topo_mkgeom',    &
     &     'topo_id was reset for non-null decomposition',__LINE__, info)
         topo_id = -1
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Dummy arguments for non-existing particles and meshes
      !-------------------------------------------------------------------------
      xpdummy(1,1)  = 0.0_MK
      Nmdummy(1)    = 0
      nnodes(1:3,1) = 0

      !-------------------------------------------------------------------------
      !  Recursive bisection
      !-------------------------------------------------------------------------
      IF (decomp.EQ.ppm_param_decomp_bisection) THEN
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK 
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
     &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
     &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &           'Bisection decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Pencils
      !-------------------------------------------------------------------------
      ELSEIF ((decomp .EQ. ppm_param_decomp_xpencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_ypencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_zpencil)) THEN
         IF (decomp.EQ.ppm_param_decomp_zpencil.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
     &           'Cannot make z pencils in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  pencil quadrisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a quad tree, binary in 2d
         treetype         = ppm_param_tree_quad
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK 
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! fix the proper direction
         fixed(1:ppm_dim) = .FALSE.
         IF (decomp .EQ. ppm_param_decomp_xpencil) fixed(1) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_ypencil) fixed(2) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_zpencil) fixed(3) = .TRUE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
     &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
     &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &           'Pencil decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Slabs
      !-------------------------------------------------------------------------
      ELSEIF ((decomp .EQ. ppm_param_decomp_xy_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_xz_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_yz_slab)) THEN
         IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
     &           'Cannot make x-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
     &           'Cannot make y-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK 
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! fix the proper directions
         fixed(1:ppm_dim) = .FALSE.
         IF (decomp .EQ. ppm_param_decomp_xy_slab) THEN
             fixed(1) = .TRUE.
             fixed(2) = .TRUE.
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_xz_slab) THEN
             fixed(1) = .TRUE.
             fixed(3) = .TRUE.
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_yz_slab) THEN
             fixed(2) = .TRUE.
             fixed(3) = .TRUE.
         ENDIF
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
     &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
     &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',    &
     &           'Slab decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Cuboids
      !-------------------------------------------------------------------------
      ELSEIF (decomp .EQ. ppm_param_decomp_cuboid) THEN
         !-------------------------------------------------------------------
         !  cuboid octasection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build an oct tree in 3d
         treetype         = ppm_param_tree_oct
         ! and a quad tree in 2d
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_quad
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK 
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
     &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
     &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &           'Cuboid decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  User provides decomposition: Do nothing
      !-------------------------------------------------------------------------
      ELSEIF (decomp .EQ. ppm_param_decomp_user_defined) THEN
          ! NOP

      !-------------------------------------------------------------------------
      !  Unknown decomposition type
      !-------------------------------------------------------------------------
      ELSE
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown decomposition type: ',decomp
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
     &       mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
     &                    min_sub,max_sub,nsubs,nneigh,ineigh,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &        'Finding neighbors failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          CALL ppm_topo_cost(xpdummy,0,min_sub,max_sub,nsubs,nnodes,  &
     &        cost,info)
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &           'Computing costs failed',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign the subdomains to processors
      !-------------------------------------------------------------------------
      IF     (assig .EQ. ppm_param_assign_internal) THEN
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, & 
     &       isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &         'Assigning subs to processors failed',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSEIF (assig .EQ. ppm_param_assign_nodal_cut .OR.    &
     &        assig .EQ. ppm_param_assign_nodal_comm .OR.   &
     &        assig .EQ. ppm_param_assign_dual_cut .OR.     &
     &        assig .EQ. ppm_param_assign_dual_comm) THEN
         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_sub,max_sub,nneigh,ineigh,cost,nsubs,&
     &       assig,sub2proc,isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &         'Assigning subs to processors using METIS failed',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSEIF (assig .EQ. ppm_param_assign_user_defined) THEN
         !-------------------------------------------------------------------
         !  user defined assignment
         !-------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mkgeom',   &
     &           'list of local subs ISUBLIST',__LINE__,info)
             GOTO 9999
         ENDIF
         isublist = ppm_param_undefined
         nsublist = 0
         DO isub=1,nsubs
             IF (sub2proc(isub) .EQ. ppm_rank) THEN
                 nsublist = nsublist + 1
                 isublist(nsublist) = isub
             ENDIF
         ENDDO
      ELSE
         !-------------------------------------------------------------------
         !  unknown assignment scheme
         !-------------------------------------------------------------------
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown assignment scheme: ',assig
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
     &       mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local 
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef,min_sub,max_sub, &
     &                        nsubs,subs_bc,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &       'finding and defining the BC of the subs failed ',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the topology internally
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                    sub2proc,nsubs,bcdef, &
     &                    isublist,nsublist,&
     &                    nneigh,ineigh,topo_id,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &        'Storing topology failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get internal topoid. If no topology has ever been defined, set 
      !  this as the current topology. 
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF (ppm_topoid .EQ. -1) ppm_topoid = topoid

      !-------------------------------------------------------------------------
      !  Dump out disgnostic files
      !-------------------------------------------------------------------------
      !IF (ppm_debug .GT. 0) THEN
      !    WRITE(mesg,'(A,I4.4)') 'part',ppm_rank
      !    OPEN(10,FILE=mesg)
      !
      !    DO j=1,nsublist
      !        i = isublist(j)
      !  
      !        ! x-y plan
      !        WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
      !        WRITE(10,'(2e12.4)') max_sub(1,i),min_sub(2,i)
      !        WRITE(10,'(2e12.4)') max_sub(1,i),max_sub(2,i)
      !        WRITE(10,'(2e12.4)') min_sub(1,i),max_sub(2,i)
      !        WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
      !        WRITE(10,'(   a  )')
      !
      !        ! y-z plan
      !        IF (ppm_dim .GT. 2) THEN
      !            WRITE(10,'(2e12.4)') min_sub(2,i),min_sub(3,i)
      !            WRITE(10,'(2e12.4)') max_sub(2,i),min_sub(3,i)
      !            WRITE(10,'(2e12.4)') max_sub(2,i),max_sub(3,i)
      !            WRITE(10,'(2e12.4)') min_sub(2,i),max_sub(3,i)
      !            WRITE(10,'(2e12.4)') min_sub(2,i),min_sub(3,i)
      !            WRITE(10,'(   a  )')
      !        ENDIF
      !    ENDDO
      !
      !    CLOSE(10)
      !ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_mkgeom',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkgeom_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkgeom_d
#endif
