      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mkfield
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is the topology routine for meshes.
      !                 Valid decomposition types are:
      !                     ppm_param_decomp_bisection
      !                     ppm_param_decomp_xpencil 
      !                     ppm_param_decomp_ypencil 
      !                     ppm_param_decomp_zpencil 
      !                     ppm_param_decomp_xy_slab
      !                     ppm_param_decomp_xz_slab
      !                     ppm_param_decomp_yz_slab
      !                     ppm_param_decomp_cartesian
      !                     ppm_param_decomp_cuboid
      !                     ppm_param_decomp_user_defined
      !
      !                 In the user_defined case, the user must submit 
      !                 existing subdomains and all of min_sub, max_sub, 
      !                 cost, nsubs must be provided. The sub boundaries
      !                 must align with mesh planes.
      !                 For all other decompositions, particles can be used
      !                 to guide them. For this, Npart must be .GT. 0. For
      !                 Npart .LE. 0, the decomposition is purely
      !                 mesh-based.
      !                 Subs are then mapped onto processors and
      !                 both the topology and the mesh definition is
      !                 stored.
      !
      !  Input        : xp(:,:)      (F) particle positions. If present,
      !                                  domain decomposition will be
      !                                  guided by particle locations
      !                 Npart        (I) number of particles. <=0 if no xp
      !                                  is given to guide decomp
      !                 Nm(:)        (I) number of mesh POINTS (not cells)
      !                                  in each direction of the global 
      !                                  computational domain (including 
      !                                  points ON its boundaries) 
      !                 decomp       (I) type of decomposition (see above)
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
      !                 ghostsize(:) (I) the size (width) of the GL in
      !                                  units of grid spacings in all 3
      !                                  directions.
      !                 ndom         (I) number of subdomains requested.
      !                                  OPTIONAL. If not given, the number
      !                                  of subs will be equal to the
      !                                  number of processors. This is only
      !                                  relevent for decomp_cartesian and
      !                                  pencils without particles.
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
      !                                  available topology ID). Topology 0
      !                                  is reserved for the ring topology
      !                                  (null decomposition) and may only
      !                                  be used with particles.
      !                                  On return: actually used or
      !                                  created topology ID.
      !                 mesh_id      (I) mesh identifier (user numbering)
      !                                  for default mesh (as defined by
      !                                  Nm) on the topology. If <0, ppm
      !                                  will create a number internally
      !                                  and return it here on exit. 
      !                 min_sub(:,:) (F) the min. extent of the subdomains.
      !                                  (either user-specified on input or
      !                                  decomposition result on output)
      !                 max_sub(:,:) (F) the max. extent of the subdomains
      !                 sub2proc(:)  (I) processor affiliation of each
      !                                  subdomain. user-defined on
      !                                  input or library result on
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
      !                                  index) in global mesh.
      !                 ndata(:,:)   (I) number of grid POINTS in x,y[,z]
      !                                  (first index) of sub mesh isub
      !                                  (second index). Includes the
      !                                  points ON the sub boundaries.
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      : In the case of user-defined subs, a lot is
      !                 currently trusted. We may want to include further
      !                 checks in that case, e.g. if proper costs are provided,
      !                 idata and ndata actually make sense, the extents of
      !                 all subs are integer multiples of the mesh
      !                 spacing, etc...
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_mkfield.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.50  2006/06/21 17:29:43  michaebe
      !  Memory leaks 0000024 and 0000023 fixed.
      !
      !  Revision 1.49  2006/03/28 19:34:13  ivos
      !  Removed debug headers.
      !
      !  Revision 1.48  2006/03/28 12:20:24  ivos
      !  Fixed broken argument declaration.
      !
      !  Revision 1.47  2006/03/24 18:29:36  ivos
      !  fixed bug 0000017. Argument check failed if xp was not ASSOCIATED.
      !  Added check for the association status.
      !
      !  Revision 1.46  2006/02/03 09:34:01  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.45  2005/05/20 17:53:35  ivos
      !  fix: the tree decompositions now respect ndom, if present.
      !
      !  Revision 1.44  2005/02/01 16:39:26  ivos
      !  bugfix: fixed wrong symbolic constant in argument check.
      !
      !  Revision 1.43  2005/02/01 15:29:44  ivos
      !  bugfix: info test in assign_user_defined was wrong.
      !
      !  Revision 1.42  2005/01/27 09:24:25  ivos
      !  minboxsize .EQ. 0 is now allowed.
      !
      !  Revision 1.41  2005/01/13 13:20:05  ivos
      !  Added possiblity of user-defined sub2proc assignment.
      !
      !  Revision 1.40  2004/12/03 17:18:50  ivos
      !  Changed to use 2d array weights since ppm_tree has changed.
      !
      !  Revision 1.39  2004/11/30 15:05:20  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.38  2004/11/05 14:15:17  ivos
      !  Commented file output of subs.
      !
      !  Revision 1.37  2004/10/05 09:04:13  ivos
      !  security fix: topoid=0 is reserved for the ring. ppm_topo_store relies
      !  on this when determining the internal topoid. Added check that topoid 0
      !  is not used for non-ring topologies. If it is, it will be reset to -1
      !  and a new topology will be created. Updated the comment header
      !  accordingly.
      !
      !  Revision 1.36  2004/09/29 11:06:23  ivos
      !  bugfix: empty line in gnuplot output added.
      !
      !  Revision 1.35  2004/09/24 15:02:35  ivos
      !  Added particle-guided field decompositions using ppm_tree.
      !
      !  Revision 1.34  2004/09/17 12:04:39  ivos
      !  fixed argument check for Nm. Must be .GE.2 and not only .GE.0.
      !
      !  Revision 1.33  2004/07/26 14:13:06  ivos
      !  Changed name of ppm_topo_subs2proc and ppm_topo_metis_s2p after
      !  renaming them.
      !
      !  Revision 1.32  2004/07/26 11:49:56  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.31  2004/07/26 08:54:24  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.30  2004/07/26 07:42:35  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.29  2004/07/16 14:47:22  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.28  2004/07/07 11:36:00  ivos
      !  Mesh is now created by ppm (in order to guarantee compatibility) even
      !  for user-defined decompositions.
      !
      !  Revision 1.27  2004/07/01 10:46:06  ivos
      !  bugfix: ppm_topo_cost and ppm_mesh_on_subs (in mkfield) were also
      !  called for decomp.EQ.ppm_param_decomp_user_defined. They should not.
      !
      !  Revision 1.26  2004/04/20 10:57:22  walther
      !  Added min_phys and max_phys to ppm_topo_store().
      !
      !  Revision 1.25  2004/04/15 10:20:19  ivos
      !  bugfix: ppm_topoid was not set if this is the first topology the user
      !  defines.
      !
      !  Revision 1.24  2004/04/14 15:03:48  ivos
      !  Added checks for decomp to disallow ring topology for meshes.
      !
      !  Revision 1.23  2004/04/05 12:23:52  walther
      !  Added the subs_bc and call to ppm_decomp_define_subs_bc.
      !
      !  Revision 1.22  2004/04/02 10:38:49  walther
      !  ppm_topo_store() now requires bcdef being passed.
      !
      !  Revision 1.21  2004/04/01 14:54:59  walther
      !  Removed the nsubsplus again.
      !
      !  Revision 1.20  2004/03/03 09:11:11  ivos
      !  Changed CALL to ppm_decomp_cartesian to new argument list. Added CALLs to
      !  ppm_topo_cost and ppm_mesh_on_subs.
      !
      !  Revision 1.19  2004/02/25 17:26:11  ivos
      !  Added cost to argument list of ppm_decomp_metis_s2p.
      !
      !  Revision 1.18  2004/02/25 14:45:37  ivos
      !  Added new METIS assignment types ppm_param_assign_nodal_cut, 
      !  _nodal_comm, _dual_cut and _dual_comm.
      !
      !  Revision 1.17  2004/02/24 14:15:21  ivos
      !  Changed CALLs to ppm_mesh_alloc to new interface. Now has the same 
      !  argument list as ppm_alloc for regular arrays.
      !
      !
      !  Revision 1.16  2004/02/20 15:46:53  ivos
      !  bugfix: ppm_topoid and ppm_field_topoid are now correctly set the first
      !  time we create any topology.
      !
      !  Revision 1.15  2004/02/19 18:01:46  ivos
      !  Now uses pointers to avoid explicit copying of data from arrays.
      !
      !  Revision 1.14  2004/02/19 14:20:27  ivos
      !  Added routine ppm_mesh_alloc for (re)allocation of mesh number list
      !  and mesh definition user-type arrays. The corresponding code has been
      !  removed from ppm_topo_mkfield and ppm_mesh_store and the Makefile
      !  updated. The new routines are in the ppm_module_mesh.
      !
      !  Revision 1.13  2004/02/12 17:07:11  ivos
      !  Added possibility for different subs2proc assignment types and added
      !  METIS-based assignment.
      !
      !  Revision 1.12  2004/02/11 14:33:26  ivos
      !  Changed CALL to ppm_topo_store to include Nm in the argument list.
      !
      !  Revision 1.11  2004/02/06 16:40:17  ivos
      !  Bugfix: member array pointers of ppm_meshid are now nullified after
      !  mother array is allocated.
      !
      !  Revision 1.10  2004/02/04 17:22:48  ivos
      !  Added cartesian decomposition, reset and storage of mesh definitions.
      !
      !  Revision 1.9  2004/01/29 14:44:28  ivos
      !  Nuked param_decomp_cartesian.
      !
      !  Revision 1.8  2004/01/28 13:42:11  michaebe
      !  added some input parameter testing for gt 0 debug level
      !
      !  Revision 1.7  2004/01/23 17:24:18  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.6  2004/01/23 15:07:23  michaebe
      !  Corrected header infor to describe what is actually implemented.
      !
      !  Revision 1.5  2004/01/23 11:31:24  ivos
      !  Cleanup: (1) updated headers, (2) inserted ppm_error and ppm_write,
      !  (3) added argument checking, (4) added checks after every alloc.
      !
      !  Revision 1.4  2003/12/19 12:42:29  michaebe
      !  corrected log entries.
      !
      !  Revision 1.3  2003/12/19 12:32:17  michaebe
      !  wrote the part where user-defined subdomains used. added a marginal 
      !  check of the latter (volume consistency).
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_mkfield_s(xp,Npart,Nm,decomp,assig,min_phys,    &
     &               max_phys,bcdef,ghostsize,topo_id,mesh_id,min_sub,    &
     &               max_sub,cost,sub2proc,nsubs,isublist,nsublist,istart,&
     &               ndata,info,ndom,pcost)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkfield_d(xp,Npart,Nm,decomp,assig,min_phys,    &
     &               max_phys,bcdef,ghostsize,topo_id,mesh_id,min_sub,    &
     &               max_sub,cost,sub2proc,nsubs,isublist,nsublist,istart,&
     &               ndata,info,ndom,pcost)
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
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_define_subs_bc
      USE ppm_module_mesh_on_subs
      USE ppm_module_mesh_alloc
      USE ppm_module_mesh_store
      USE ppm_module_topo_subs2proc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_find_neigh
      USE ppm_module_decomp_cartesian
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
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(IN   ) :: assig
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN   ) :: pcost
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: ghostsize
      INTEGER , OPTIONAL      , INTENT(IN   ) :: ndom
      INTEGER                 , INTENT(INOUT) :: nsubs,topo_id,mesh_id,decomp
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart,ndata
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                 , INTENT(  OUT) :: nsublist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER                           :: i,topoid,ul,iopt,j,treetype,nbox
      INTEGER                           :: isub,minbox
      INTEGER, DIMENSION(1)             :: ldc
      INTEGER, DIMENSION(:,:), POINTER  :: ineigh,subs_bc
      INTEGER, DIMENSION(:  ), POINTER  :: nneigh,nchld
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps,maxvar
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec,meshdx
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
      CALL substart('ppm_topo_mkfield',t0,info)
#if    __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif  __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      ! do not care about the variance for meshes. Stop as soon as you have
      ! enough subdomains.
      maxvar = HUGE(maxvar)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF(ppm_debug.GT.0) THEN
         IF (.NOT. ppm_initialized) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mkfield',  &
     &           'Please call ppm_init first!',__LINE__,info)
             GOTO 9999
         ENDIF
         DO i=1,ppm_dim
            IF(max_phys(i).LE.min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'max_phys must be > min_phys',__LINE__, info)
               GOTO 9999
            ENDIF
            IF(Nm(i) .LT. 2) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'Nm must be > 1 in all dimensions',__LINE__, info)
               GOTO 9999
            ENDIF
            IF(ghostsize(i) .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'ghostsize must be >= 0',__LINE__, info)
               GOTO 9999
            ENDIF
         ENDDO
         IF (Npart .GT. 0 .AND. ASSOCIATED(xp)) THEN
            IF(SIZE(xp,2) .LT. Npart) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'not enough particles contained in xp',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF(SIZE(xp,1) .LT. ppm_dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'leading dimension of xp too small',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF (PRESENT(pcost)) THEN
                IF(SIZE(pcost,1) .LT. Npart) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &                 'pcost does not contain costs for all particles',&
     &                 __LINE__, info)
                   GOTO 9999
                ENDIF
            ENDIF
         ENDIF
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (decomp .NE. ppm_param_decomp_user_defined) THEN
               info = ppm_error_warning
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'decomp type is set to user_defined for this assignment',&
     &             __LINE__, info)
               decomp = ppm_param_decomp_user_defined
            ENDIF
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 9999
            ENDIF
            DO i=1,nsubs
                IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 9999
                ENDIF
            ENDDO
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'no subs defined in user_defined decomposition',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF ((.NOT.ASSOCIATED(min_sub)).OR.(.NOT.ASSOCIATED(max_sub))) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &              'min_sub/max_sub must be allocated for user def. decomp',&
     &              __LINE__, info)
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
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',   &
     &              'faulty subdomains defined',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF      
         
         !----------------------------------------------------------------------
         ! Check bcdef
         !----------------------------------------------------------------------
         !  [TODO]!!!
         !  must check for compatibility of boundary conditions
         !  and boundary condtions array must have one more 
         !  dimension, as its possible that you dont want to 
         !  have the same boundary conditions on every side
         !  of the subdomain
         !  Will be accounted for as soon as the structure of
         !  the field type is better defined
         !----------------------------------------------------------------------
      ENDIF

      IF (topo_id .EQ. 0) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_argument, 'ppm_topo_mkfield',    &
     &     'topo_id was reset for non-null decomposition',__LINE__, info)
         topo_id = -1
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Compute grid spacing
      !-------------------------------------------------------------------------
      meshdx(1) = (max_phys(1)-min_phys(1))/REAL(Nm(1)-1,MK)
      meshdx(2) = (max_phys(2)-min_phys(2))/REAL(Nm(2)-1,MK)
      IF (ppm_dim .GT. 2) THEN
          meshdx(3) = (max_phys(3)-min_phys(3))/REAL(Nm(3)-1,MK)
      ENDIF

      !-------------------------------------------------------------------------
      !  Cartesian (mesh-only) domain decomposition
      !-------------------------------------------------------------------------
      IF (decomp .EQ. ppm_param_decomp_cartesian) THEN
          IF (PRESENT(ndom)) THEN
              CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &            ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info,ndom)
          ELSE
              CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &            ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info)
          ENDIF
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &            'Cartesian decomposition failed',__LINE__,info)
              GOTO 9999
          ENDIF

      !-------------------------------------------------------------------------
      !  Recursive bisection. Can be particle-guided
      !-------------------------------------------------------------------------
      ELSEIF (decomp.EQ.ppm_param_decomp_bisection) THEN
         !-------------------------------------------------------------------
         !  recursive bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype           = ppm_param_tree_bin
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1)         = REAL(ghostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ghostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)     = REAL(ghostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
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
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make z pencils in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (Npart .LT. 1) THEN
             !------------------------------------------------------------------
             !  No particles: cartesian pencil decomposition
             !------------------------------------------------------------------
             IF (PRESENT(ndom)) THEN
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &               decomp,min_sub,max_sub,nsubs,info,ndom)
             ELSE
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &               decomp,min_sub,max_sub,nsubs,info)
             ENDIF
             IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &               'Cartesian pencil decomposition failed',__LINE__,info)
                 GOTO 9999
             ENDIF
         ELSE
             !------------------------------------------------------------------
             !  With particles: pencil quadrisection using ppm_tree
             !------------------------------------------------------------------
             ! build a quad tree, binary in 2d
             treetype           = ppm_param_tree_quad
             IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
             IF (Npart .GT. 0) THEN
                 weights(1,1:2) = 0.5_MK    ! particles have 50% weight
                 weights(2,1:2) = 0.5_MK
             ELSE
                 weights(1,1:2) = 0.0_MK    ! only mesh has weight
                 weights(2,1:2) = 1.0_MK
             ENDIF
             ! geometry has no weight
             weights(3,1:2)     = 0.0_MK
             ! fix the proper direction
             fixed(1:ppm_dim) = .FALSE.
             IF (decomp .EQ. ppm_param_decomp_xpencil) fixed(1) = .TRUE.
             IF (decomp .EQ. ppm_param_decomp_ypencil) fixed(2) = .TRUE.
             IF (decomp .EQ. ppm_param_decomp_zpencil) fixed(3) = .TRUE.
             gsvec(1)         = REAL(ghostsize(1),MK)*meshdx(1)
             gsvec(2)         = REAL(ghostsize(2),MK)*meshdx(2)
             IF (ppm_dim .GT. 2) THEN
                 gsvec(3)     = REAL(ghostsize(3),MK)*meshdx(3)
             ENDIF
             minbox = ppm_nproc
             IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
             ! build tree
             IF (PRESENT(pcost)) THEN
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,      &
     &               minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
     &               min_box,max_box,nbox,nchld,info,pcost)
             ELSE
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,      &
     &               minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
     &               min_box,max_box,nbox,nchld,info)
             ENDIF
             IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &               'Pencil decomposition failed',__LINE__,info)
                 GOTO 9999
             ENDIF
             ! convert tree to subs
             CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &           max_sub,nsubs,info)
             IF (info .NE. 0) GOTO 9999
         ENDIF

      !-------------------------------------------------------------------------
      !  Slabs
      !-------------------------------------------------------------------------
      ELSEIF ((decomp .EQ. ppm_param_decomp_xy_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_xz_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_yz_slab)) THEN
         IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make x-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make y-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
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
         gsvec(1)         = REAL(ghostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ghostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)     = REAL(ghostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',    &
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
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1)         = REAL(ghostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ghostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)         = REAL(ghostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
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
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',   &
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
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &        'Finding neighbors failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Define meshes on the subs
      !-------------------------------------------------------------------------
      CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,min_sub,max_sub,nsubs, &
     &                      istart,ndata,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &       'Defining meshes failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          IF (PRESENT(pcost)) THEN
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,  &
     &                            info,pcost)
          ELSE
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,info)
          ENDIF
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
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
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
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
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &         'Assigning subs to processors using METIS failed',__LINE__,info)
            GOTO 9999
         ENDIF
     ELSEIF (assig .EQ. ppm_param_assign_user_defined) THEN
         !--------------------------------------------------------------------
         !  user defined assignment
         !--------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mkfield',   &
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
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',   &
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
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
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
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &        'Storing topology failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get internal topoid. If no field topology has ever been defined, set 
      !  this as the current field topology. If no particle topology has
      !  been defined so far, this will also be the current particle
      !  topology.
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF (ppm_topoid .EQ. -1) ppm_topoid = topoid
      IF (ppm_field_topoid .EQ. -1) ppm_field_topoid = topoid

      !-------------------------------------------------------------------------
      !  (Re)allocate mesh counter if number of topologies increased
      !-------------------------------------------------------------------------
      ! just call grow_preserve, it will internally check if realloc is
      ! actually needed or the size is still the same...
      iopt   = ppm_param_alloc_grow_preserve
      ldc(1) = ppm_max_topoid
      CALL ppm_alloc(ppm_max_meshid,ldc,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_topo_mkfield',  &
     &       'maximum mesh id MAXMESHID',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Grow ppm_meshid if needed. The following is basically a
      !  ppm_alloc_fit_preserve for ppm_meshid(ppm_max_topoid).
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit_preserve
      ldc(1) = ppm_max_topoid
      CALL ppm_mesh_alloc(ppm_meshid,ldc,iopt,info)
      IF (info .NE. 0) THEN
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',   &
     &        'Growing ppm_meshid list failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Flush mesh definitions of this topology
      !-------------------------------------------------------------------------
      ppm_max_meshid(topoid) = 0
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ppm_meshid(topoid)%user,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield',  &
     &        'User mesh number list PPM_MESHID%USER',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ppm_meshid(topoid)%internal,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield',  &
     &        'Internal mesh number list PPM_MESHID%INTERNAL',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh internally
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(mesh_id,topoid,nsubs,ndata,istart,Nm,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &        'Storing mesh definition failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Dump out disgnostic files
      !-------------------------------------------------------------------------
!      IF (ppm_debug .GT. 0) THEN
!          WRITE(mesg,'(A,I4.4)') 'part',ppm_rank
!          OPEN(10,FILE=mesg)
!          DO ul=1,nsublist
!             i = isublist(ul)
!
!    ! x-y plan
!            WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
!            WRITE(10,'(2e12.4)') max_sub(1,i),min_sub(2,i)
!            WRITE(10,'(2e12.4)') max_sub(1,i),max_sub(2,i)
!            WRITE(10,'(2e12.4)') min_sub(1,i),max_sub(2,i)
!            WRITE(10,'(2e12.4)') min_sub(1,i),min_sub(2,i)
!            WRITE(10,'(   a  )')
!
!    ! y-z plan
!            IF (ppm_dim .GT. 2) THEN
!              WRITE(10,'(2e12.4)') min_sub(2,i),min_sub(3,i)
!              WRITE(10,'(2e12.4)') max_sub(2,i),min_sub(3,i)
!              WRITE(10,'(2e12.4)') max_sub(2,i),max_sub(3,i)
!              WRITE(10,'(2e12.4)') min_sub(2,i),max_sub(3,i)
!              WRITE(10,'(2e12.4)') min_sub(2,i),min_sub(3,i)
!              WRITE(10,'(   a  )')
!            ENDIF
!          ENDDO
!          CLOSE(10)
!      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ineigh,ldc,iopt,info)
      CALL ppm_alloc(nneigh,ldc,iopt,info)
      IF(decomp.NE.ppm_param_decomp_cartesian.AND.&
     & decomp.NE.ppm_param_decomp_user_defined) THEN
         CALL ppm_alloc(min_box,ldc,iopt,info)
         CALL ppm_alloc(max_box,ldc,iopt,info)
         CALL ppm_alloc(nchld,ldc,iopt,info)
      END IF
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield',     &
     &        'deallocation failed',__LINE__,info)
      ENDIF

9999  CONTINUE
      CALL substop('ppm_topo_mkfield',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_d
#endif
