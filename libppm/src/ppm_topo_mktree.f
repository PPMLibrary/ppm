      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mktree
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine creates a topology based on ppm_tree.
      !                 The user can directly specify the arguments for
      !                 ppm_tree, rather than choosing from a set of
      !                 predefined ppm_param_decomp_*. Please see header of
      !                 ppm_tree for description of the arguments. Beware:
      !                 this is for advanced users only.
      !
      !  Input        : xp(:,:)      (F) the particle positions
      !                 Npart        (I) the number of particles. Set to
      !                                  something .LE.0 if there are no
      !                                  particles.
      !                 Nm(:)        (I) global number of mesh POINTS (not
      !                                  cells) in all directions. Set to
      !                                  (0,0,0) if there is no mesh. If a
      !                                  mesh is present, the sub
      !                                  boundaries will be aligned with
      !                                  mesh planes.
      !                 storemesh    (L) .TRUE. to store the mesh
      !                                  definition internally. If .FALSE.,
      !                                  the mesh serves only as a guide
      !                                  for the decomposition and is not
      !                                  stored as a ppm compute mesh.
      !                 min_phys(:)  (F) the minimum coordinate of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinate of the 
      !                                  physical/computational domain 
      !                 treetype     (I) One of:
      !                                      ppm_param_tree_bin
      !                                      ppm_param_tree_quad
      !                                      ppm_param_tree_oct (3D only)
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
      !                 bcdef(:)     (I) the definition of the BC 
      !                 minboxsize   (F) the min size of the boxes
      !                 pruneboxes   (L) keep empty boxes or not?
      !                 weights(3,2) (F) weights for particles, mesh points
      !                                  and geometry. weights(:,1)
      !                                  defines the box costs (where to
      !                                  refine), weights(:,2) defines
      !                                  where to cut a box.
      !                 fixed(:)     (L) flags to prevent cuts in certain
      !                                  directions. (see ppm_tree)
      !                 maxvariance  (F) Maximum allowed variance of costs
      !                                  between subs
      !                 maxboxcost   (F) Maximum cost of a box.
      !                                  Subdivision will continue (if
      !                                  possible) if cost of any box is
      !                                  larger than this. Set to .LE. 0
      !                                  to disable this criterion.
      !                 pcost(:)     (F) OPTIONAL. Computational cost of
      !                                  each particle (if present). Of
      !                                  length Npart.
      !
      !  Input/output : topo_id      (I) topology identifier (user
      !                                  numbering) for which to create a
      !                                  mesh. Can be an existing one (in
      !                                  which case it is REPLACED), a new
      !                                  one or <0 (in which case ppm
      !                                  internally generates the next
      !                                  available topology ID).
      !                                  On return: actually used or
      !                                  created topology ID. Topology 0
      !                                  is reserved for the ring topology
      !                                  (null decomposition) and may only
      !                                  be used in conjunction with
      !                                  ppm_param_decomp_null and with
      !                                  particles.
      !                 mesh_id      (I) mesh identifier (user numbering)
      !                                  for default mesh (as defined by
      !                                  Nm) on the topology. If <0, ppm
      !                                  will create a number internally
      !                                  and return it here on exit. This
      !                                  parameter is only meaningful if
      !                                  storemesh is .TRUE.
      !                 min_sub(:,:) (F) the min. extent of the subdomains.
      !                                  (either user-specified on input or
      !                                  decomposition result on output)
      !                 max_sub(:,:) (F) the max. extent of the subdomains
      !                 sub2proc(:)  (I) processor affiliation of each
      !                                  subdomain. User-defined on
      !                                  input, library result on
      !                                  output.
      !                 nsubs        (I) the total number of subdomains (on
      !                                  all processors in total). Either
      !                                  user-specified on input or
      !                                  decomposoiton result on output.
      !                 cost(:)      (F) estimated cost associated with
      !                                  subdomains. Either user-specified
      !                                  on input or decomp. result. 
      !
      !  Output       : isublist(:)  (I) list of subdomains handled by the
      !                                  local processor
      !                 nsublist     (I) number of subdomains handled by
      !                                  the local processor
      !                 istart(:,:)  (I) start indices (i,j[,k]) (first
      !                                  index) of sub mesh isub (second
      !                                  index) in global mesh. This is
      !                                  only returned if there is a mesh!
      !                 ndata(:,:)   (I) number of grid POINTS in x,y[,z]
      !                                  (first index) of sub mesh isub
      !                                  (second index). Includes the
      !                                  points ON the sub boundaries. This
      !                                  is only returned if there is a
      !                                  mesh!
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_mktree.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.14  2006/09/04 18:34:56  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.13  2006/03/28 19:33:52  ivos
      !  Changed xp from INTENT(IN) to POINTER to be consistent with mkpart
      !  and mkfield.
      !
      !  Revision 1.12  2006/02/03 09:34:02  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.11  2005/02/01 16:39:27  ivos
      !  bugfix: fixed wrong symbolic constant in argument check.
      !
      !  Revision 1.10  2005/02/01 15:29:44  ivos
      !  bugfix: info test in assign_user_defined was wrong.
      !
      !  Revision 1.9  2005/01/13 13:20:05  ivos
      !  Added possiblity of user-defined sub2proc assignment.
      !
      !  Revision 1.8  2005/01/05 16:28:57  ivos
      !  bugfix: mesh can now only be stored if there is one at all.
      !
      !  Revision 1.7  2004/12/03 17:18:50  ivos
      !  Changed to use 2d array weights since ppm_tree has changed.
      !
      !  Revision 1.6  2004/11/30 15:05:21  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.5  2004/11/05 14:15:18  ivos
      !  Commented file output of subs.
      !
      !  Revision 1.4  2004/11/03 16:25:35  ivos
      !  Added pruneboxes and minboxsize to argument list.
      !
      !  Revision 1.3  2004/10/05 09:04:13  ivos
      !  security fix: topoid=0 is reserved for the ring. ppm_topo_store relies
      !  on this when determining the internal topoid. Added check that topoid 0
      !  is not used for non-ring topologies. If it is, it will be reset to -1
      !  and a new topology will be created. Updated the comment header
      !  accordingly.
      !
      !  Revision 1.2  2004/09/30 10:49:42  ivos
      !  bugfix: SIZE of Nm is now checked before accessing it when checking
      !  if there is a mesh at all. This allows the user to pass dummy pointers.
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
      SUBROUTINE ppm_topo_mktree_s(xp,Npart,Nm,storemesh,min_phys,max_phys,  &
     &    treetype,assig,bcdef,minboxsize,pruneboxes,weights,fixed,          &
     &    maxvariance,maxboxcost,topo_id,mesh_id,min_sub,max_sub,nsubs,      &
     &    cost,sub2proc,isublist,nsublist,istart,ndata,info,pcost)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mktree_d(xp,Npart,Nm,storemesh,min_phys,max_phys,  &
     &    treetype,assig,bcdef,minboxsize,pruneboxes,weights,fixed,          &
     &    maxvariance,maxboxcost,topo_id,mesh_id,min_sub,max_sub,nsubs,      &
     &    cost,sub2proc,isublist,nsublist,istart,ndata,info,pcost)
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mesh_alloc
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_mesh_store
      USE ppm_module_mesh_on_subs
      USE ppm_module_define_subs_bc
      USE ppm_module_topo_subs2proc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_find_neigh
      USE ppm_module_tree
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
      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      INTEGER                 , INTENT(IN   ) :: Npart,assig,treetype
      LOGICAL                 , INTENT(IN   ) :: storemesh,pruneboxes
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef,Nm
      REAL(MK)                , INTENT(IN   ) :: minboxsize,maxvariance
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      REAL(MK), DIMENSION(3,2), INTENT(IN   ) :: weights
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: pcost
      INTEGER                 , INTENT(INOUT) :: topo_id,mesh_id
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      INTEGER , DIMENSION(:,:), POINTER       :: istart,ndata
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                 , INTENT(  OUT) :: nsublist,nsubs
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                           :: i,topoid,j,nbox,iopt,isub
      INTEGER, DIMENSION(2  )           :: ldc
      INTEGER, DIMENSION(:,:), POINTER  :: ineigh,subs_bc
      INTEGER, DIMENSION(:  ), POINTER  :: nneigh,nchld
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(:,:), POINTER :: min_box,max_box
      CHARACTER(LEN=ppm_char)           :: mesg
      LOGICAL                           :: have_particles,have_mesh
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_mktree',t0,info)
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
             CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mktree',  &
     &           'Please call ppm_init first!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF(minboxsize .LT. 0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_topo_mktree', &
     &          'minboxsize must be >= 0.0',__LINE__, info)
            GOTO 9999
         ENDIF
         DO i=1,ppm_dim
            IF(max_phys(i).LE.min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mktree', &
     &             'max_phys must be > min_phys',__LINE__, info)
               GOTO 9999
            ENDIF
         ENDDO
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mktree', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mktree', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 9999
            ENDIF
            DO i=1,nsubs
                IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mktree', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 9999
                ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF (topo_id .EQ. 0) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_argument, 'ppm_topo_mktree',    &
     &     'topo_id was reset for non-null decomposition',__LINE__, info)
         topo_id = -1
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Check if we have particles and mesh
      !-------------------------------------------------------------------------
      have_particles = .FALSE.
      have_mesh      = .FALSE.
      IF (Npart .GT. 0) have_particles = .TRUE.
      IF (SIZE(Nm,1) .GE. ppm_dim) THEN
          IF (ppm_dim .GT. 2) THEN
              IF ((Nm(1).GT.1).AND.(Nm(2).GT.1).AND.(Nm(3).GT.1))     &
     &            have_mesh = .TRUE.
          ELSE
              IF ((Nm(1).GT.1).AND.(Nm(2).GT.1)) have_mesh = .TRUE.
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Call the tree routine
      !-------------------------------------------------------------------------
      gsvec(1:ppm_dim) = minboxsize
      ! build tree
      IF (PRESENT(pcost)) THEN
          CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,  &
     &        ppm_nproc,pruneboxes,gsvec,maxvariance,maxboxcost,fixed,  &
     &        weights,min_box,max_box,nbox,nchld,info,pcost)
      ELSE
          CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,  &
     &        ppm_nproc,pruneboxes,gsvec,maxvariance,maxboxcost,fixed,  &
     &        weights,min_box,max_box,nbox,nchld,info)
      ENDIF
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &        'General tree decomposition failed',__LINE__,info)
          GOTO 9999
      ENDIF
      ! convert tree to subs
      CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &    max_sub,nsubs,info)
      IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
     &                    min_sub,max_sub,nsubs,nneigh,ineigh,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &        'Finding neighbors failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate mesh data to avoid errors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = 0
      CALL ppm_alloc(istart,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_topo_mktree',  &
     &       'mesh start indices ISTART',__LINE__,info)
         GOTO 9999
      ENDIF 
      CALL ppm_alloc(ndata,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_topo_mktree',  &
     &       'mesh sizes NDATA',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Define meshes on the subs
      !-------------------------------------------------------------------------
      IF (have_mesh) THEN
          CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,min_sub,max_sub,nsubs, &
     &                          istart,ndata,info)
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &           'Defining meshes failed',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (PRESENT(pcost)) THEN
          CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,  &
     &        cost,info,pcost)
      ELSE
          CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,  &
     &        cost,info)
      ENDIF
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &       'Computing costs failed',__LINE__,info)
         GOTO 9999
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
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
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
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &         'Assigning subs to processors using METIS failed',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSEIF (assig .EQ. ppm_param_assign_user_defined) THEN
         !----------------------------------------------------------------------
         !  user defined assignment
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mktree',   &
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
         CALL ppm_error(ppm_err_argument,'ppm_topo_mktree',   &
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
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &       'finding and defining the BC of the subs failed ',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the topology internally
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                    sub2proc,nsubs,bcdef,isublist,nsublist,&
     &                    nneigh,ineigh,topo_id,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &        'Storing topology failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get internal topoid. If no topology has ever been defined, set 
      !  this as the current topology. 
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF (ppm_topoid .EQ. -1) ppm_topoid = topoid

      IF (storemesh .AND. have_mesh) THEN
          IF (ppm_field_topoid .EQ. -1) ppm_field_topoid = topoid

          !---------------------------------------------------------------------
          !  (Re)allocate mesh counter if number of topologies increased
          !---------------------------------------------------------------------
          ! just call grow_preserve, it will internally check if realloc is
          ! actually needed or the size is still the same...
          iopt   = ppm_param_alloc_grow_preserve
          ldc(1) = ppm_max_topoid
          CALL ppm_alloc(ppm_max_meshid,ldc,iopt,info)
          IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mktree',  &
     &           'maximum mesh id MAXMESHID',__LINE__,info)
             GOTO 9999
          ENDIF 

          !---------------------------------------------------------------------
          !  Grow ppm_meshid if needed. The following is basically a
          !  ppm_alloc_fit_preserve for ppm_meshid(ppm_max_topoid).
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit_preserve
          ldc(1) = ppm_max_topoid
          CALL ppm_mesh_alloc(ppm_meshid,ldc,iopt,info)
          IF (info .NE. 0) THEN
              CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',   &
     &            'Growing ppm_meshid list failed',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Flush mesh definitions of this topology
          !---------------------------------------------------------------------
          ppm_max_meshid(topoid) = 0
          iopt = ppm_param_dealloc
          CALL ppm_alloc(ppm_meshid(topoid)%user,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_topo_mktree',  &
     &            'User mesh number list PPM_MESHID%USER',__LINE__,info)
          ENDIF
          CALL ppm_alloc(ppm_meshid(topoid)%internal,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_topo_mktree',  &
     &            'Internal mesh number list PPM_MESHID%INTERNAL',__LINE__,info)
          ENDIF

          !---------------------------------------------------------------------
          !  Store new mesh internally
          !---------------------------------------------------------------------
          CALL ppm_mesh_store(mesh_id,topoid,nsubs,ndata,istart,Nm,info)
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mktree',  &
     &            'Storing mesh definition failed',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF           ! storemesh

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
      CALL substop('ppm_topo_mktree',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mktree_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mktree_d
#endif
