      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_topo_mkpart
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine is the main topology routine. It performs
      !                 the decomposition of the physical space based on the
      !                 the position of the particles. The subdomains are mapped
      !                 onto the processors and a neighbour list is established.
      !
      !                 The decomposition is based on the particle positions
      !                 If nsubs is greater than zero, the subdomains 
      !                 described by (min_sub,max_sub) given by the user is 
      !                 used. If nsubs is less than one the subdomains are 
      !                 found using the decomposition defined by the option 
      !                 decomp. 
      !                 
      !                 If nsubs on input is greater than zero, the assigment 
      !                 of subdomain to the processors is used (as described 
      !                 in sub2proc) or else 
      !                 the assigment of subdomain to the processors will be
      !                 performed to assign as best as possible an equal number
      !                 of particles/grid points and at the same time
      !                 minimizing communication (see Farhat 1988, JCP 28(5),
      !                 pages: 579 -- 602).
      !
      !  Input        : xp(:,:)      (F) the position of the particles
      !                 Npart        (I) the number of particles (>0)
      !                 decomp       (I) the type of decomposition. One of
      !                                     ppm_param_decomp_pruned_cell
      !                                     ppm_param_decomp_tree
      !                                     ppm_param_decomp_bisection
      !                                     ppm_param_decomp_xpencil
      !                                     ppm_param_decomp_ypencil
      !                                     ppm_param_decomp_zpencil
      !                                     ppm_param_decomp_xy_slab
      !                                     ppm_param_decomp_xz_slab
      !                                     ppm_param_decomp_yz_slab
      !                                     ppm_param_decomp_cuboid
      !                                     ppm_param_decomp_user_defined
      !                                     ppm_param_decomp_null
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
      !                 min_phys(:)  (F) the minimum coordinates of the 
      !                                  physical/computational domain 
      !                 max_phys(:)  (F) the maximum coordinates of the 
      !                                  physical/computational domain 
      !                 bcdef(:)     (I) the definition of the BC on all
      !                                  sides
      !                 ghostsize    (F) the size (width) of the ghost
      !                                  layer.
      !                 pcost(:)     (F) OPTIONAL array of length Npart
      !                                  which specifies the computational
      !                                  cost attributed to each particle.
      !                                  If this is absent, a unity cost is
      !                                  assumed for each particle. 
      !
      !  Input/output : topo_id      (I) topology identifier (user
      !                                  numbering) for which to create
      !                                  particles. Can be an already
      !                                  existing one (in which case the
      !                                  old topology is REPLACED), a new
      !                                  one or <0 (in which case ppm
      !                                  internally generates the next
      !                                  available topology ID. Topology 0
      !                                  is reserved for the ring topology
      !                                  (null decomposition) and may only
      !                                  be used in conjunction with
      !                                  ppm_param_decomp_null.
      !                                  On return: Actually used topology
      !                                  ID (or ID of the new topology if
      !                                  <0 on input).
      !                 min_sub(:,:) (F) the min. extent of the subdomains
      !                                  (user-specified on input or 
      !                                  decomposition-result on output)
      !                 max_sub(:,:) (F) the max. extent of the subdomains
      !                 sub2proc(:)  (I) the list proc affiliations of
      !                                  all subs. user-defined on input
      !                                  or library result on output.
      !                 nsubs        (I) the total number of subdomains (on
      !                                  all processors in total). Either
      !                                  user-specified on input or a
      !                                  result of the decomposition on
      !                                  output.
      !                 cost(:)      (F) estimated cost associated with
      !                                  subdomains. Either user-defined on
      !                                  input or decomposition result on
      !                                  output.
      !
      !  Output       : isublist(:)  (I) list of subdomains handled by 
      !                                  the local processor
      !                 nsublist     (I) the number of subs assigned to the
      !                                  current processor
      !                 info         (I) return status. 0 upon success
      !
      !  Remarks      : In the case of user-defined subs, a lot is
      !                 currently trusted. We may want to include further
      !                 checks in that case, e.g. if proper costs are provided.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_mkpart.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.47  2006/09/04 18:34:56  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.46  2006/03/28 19:34:48  ivos
      !  Added argument check for non-associated xp.
      !
      !  Revision 1.45  2006/03/28 12:19:33  ivos
      !  Fixed broken argument declaration.
      !
      !  Revision 1.44  2006/03/24 18:29:36  ivos
      !  fixed bug 0000017. Argument check failed if xp was not ASSOCIATED.
      !  Added check for the association status.
      !
      !  Revision 1.43  2006/02/03 09:34:01  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.42  2005/02/01 16:39:26  ivos
      !  bugfix: fixed wrong symbolic constant in argument check.
      !
      !  Revision 1.41  2005/02/01 15:29:44  ivos
      !  bugfix: info test in assign_user_defined was wrong.
      !
      !  Revision 1.40  2005/01/13 13:20:04  ivos
      !  Added possiblity of user-defined sub2proc assignment.
      !
      !  Revision 1.39  2004/12/03 17:18:49  ivos
      !  Changed to use 2d array weights since ppm_tree has changed.
      !
      !  Revision 1.38  2004/11/30 15:05:21  ivos
      !  Added maxboxcost and part2box to ppm_tree argument list.
      !
      !  Revision 1.37  2004/11/05 11:33:34  walther
      !  Commented the ascii output of the subs at the end of the routine.
      !
      !  Revision 1.36  2004/11/03 16:24:46  ivos
      !  fixed typo in header comment.
      !
      !  Revision 1.35  2004/10/05 09:04:13  ivos
      !  security fix: topoid=0 is reserved for the ring. ppm_topo_store relies
      !  on this when determining the internal topoid. Added check that topoid 0
      !  is not used for non-ring topologies. If it is, it will be reset to -1
      !  and a new topology will be created. Updated the comment header
      !  accordingly.
      !
      !  Revision 1.34  2004/09/24 15:03:16  ivos
      !  Added slabs, cuboids, pencils and recursive bisection using ppm_tree.
      !
      !  Revision 1.33  2004/09/23 09:50:36  ivos
      !  bugfix: test of info was wrong after ppm_topo_box2subs
      !
      !  Revision 1.32  2004/09/22 17:22:42  ivos
      !  Added recursive bisection decomposition.
      !
      !  Revision 1.31  2004/08/05 16:20:44  ivos
      !  bugfix: If domain was less than cutoff in extent, cartesian
      !  decomp. failed. This is now fixed.
      !
      !  Revision 1.30  2004/07/26 14:13:06  ivos
      !  Changed name of ppm_topo_subs2proc and ppm_topo_metis_s2p after
      !  renaming them.
      !
      !  Revision 1.29  2004/07/26 08:54:25  ivos
      !  Changed names of renamed subroutines.
      !
      !  Revision 1.28  2004/07/26 07:42:36  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.27  2004/07/16 14:47:23  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.26  2004/07/01 10:46:06  ivos
      !  bugfix: ppm_topo_cost and ppm_mesh_on_subs (in mkfield) were also
      !  called for decomp.EQ.ppm_param_decomp_user_defined. They should not.
      !
      !  Revision 1.25  2004/06/15 14:48:35  ivos
      !  Relaxed argument check of ghostsize for the decomp_null case.
      !
      !  Revision 1.24  2004/04/20 10:57:22  walther
      !  Added min_phys and max_phys to ppm_topo_store().
      !
      !  Revision 1.23  2004/04/15 16:09:43  ivos
      !  Removed subs2proc assignment from decomp_null and moved it into
      !  its case in topo_mkpart. Mainly because decomp routines are not
      !  expected to do assignments (program logic).
      !
      !  Revision 1.22  2004/04/14 15:03:18  ivos
      !  topo_id is now automatically set to 0 when called with
      !  ppm_param_decomp_null (instead of producing a fatal error).
      !
      !  Revision 1.21  2004/04/13 15:22:44  oingo
      !  Added case for ppm_param_decomp_null. topo_id must be 0 to use this case
      !
      !  Revision 1.20  2004/04/05 12:23:46  walther
      !  Added the subs_bc and call to ppm_decomp_define_subs_bc.
      !
      !  Revision 1.19  2004/04/02 10:38:49  walther
      !  ppm_topo_store() now requires bcdef being passed.
      !
      !  Revision 1.18  2004/04/01 14:54:58  walther
      !  Removed the nsubsplus again.
      !
      !  Revision 1.17  2004/03/04 14:30:20  ivos
      !  bugfix: inserted missing __ in a cpp if.
      !
      !  Revision 1.16  2004/03/03 09:10:31  ivos
      !  Changed all calls to decomp_ routines to new argument list. Added CALL
      !  to ppm_topo_cost to compute the sub costs.
      !
      !  Revision 1.15  2004/02/25 17:26:10  ivos
      !  Added cost to argument list of ppm_decomp_metis_s2p.
      !
      !  Revision 1.14  2004/02/25 14:45:36  ivos
      !  Added new METIS assignment types ppm_param_assign_nodal_cut, _nodal_comm,
      !  _dual_cut and _dual_comm.
      !
      !  Revision 1.13  2004/02/20 15:46:52  ivos
      !  bugfix: ppm_topoid and ppm_field_topoid are now correctly set the first
      !  time we create any topology.
      !
      !  Revision 1.12  2004/02/18 13:55:18  ivos
      !  Corrected header comment for cost.
      !
      !  Revision 1.11  2004/02/12 17:06:35  ivos
      !  Added possibility for different subs2proc assignment types and added
      !  METIS-based assignment.
      !
      !  Revision 1.10  2004/02/04 17:23:14  ivos
      !  Added cartesian decomposition.
      !
      !  Revision 1.9  2004/01/30 10:01:35  ivos
      !  Added call to ppm_decomp_cartesian and fixed bug: output in x-z plane
      !  is now only written if ppm_dim is .GT. 2.
      !
      !  Revision 1.8  2004/01/29 14:41:47  ivos
      !  replaced decomp_cartesian with decomp_cuboid.
      !
      !  Revision 1.7  2004/01/23 17:24:18  ivos
      !  Now includes ppm_define.h for the cpp defines and no longer sets them
      !  in the Makefile.
      !
      !  Revision 1.6  2004/01/23 15:41:58  ivos
      !  Cleanup: (1) updated header, (2) added checks after all allocations,
      !  (3) inserted ppm_error and ppm_write, (4) added argument checking.
      !
      !  Revision 1.5  2003/12/12 16:12:02  ivos
      !  Removed topoid from the argument list of all decomp* subroutines as 
      !  they do not need it. Changed topoid (internal) to topo_id 
      !  (user-provided). Translation will be done by ppm_topo_store.
      !
      !  Revision 1.4  2003/12/05 11:19:32  ivos
      !  Bigfix. It now compiles.
      !
      !  Revision 1.3  2003/12/04 17:22:58  ivos
      !  Added comments about the meanings of several neighbor lists
      !
      !  Revision 1.2  2003/12/04 16:22:11  ivos
      !  Removed nneighlist and ineighlist (neighbor processor lists) from the
      !  argument list since they were never used in the whole routine anyway.
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
      SUBROUTINE ppm_topo_mkpart_s(xp,Npart,decomp,assig,min_phys,max_phys, &
     &               bcdef,ghostsize,topo_id,min_sub,max_sub,cost,sub2proc, &
     &               nsubs,isublist,nsublist,info,pcost)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkpart_d(xp,Npart,decomp,assig,min_phys,max_phys, &
     &               bcdef,ghostsize,topo_id,min_sub,max_sub,cost,sub2proc, &
     &               nsubs,isublist,nsublist,info,pcost)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_define_subs_bc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_topo_subs2proc
      USE ppm_module_find_neigh
      USE ppm_module_decomp_cartesian
      USE ppm_module_decomp_null
      USE ppm_module_decomp_tree
      USE ppm_module_decomp_pruned_cell
      USE ppm_module_tree
      USE ppm_module_topo_box2subs
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
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
      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                 , INTENT(IN   ) :: Npart,assig
      INTEGER                 , INTENT(INOUT) :: nsubs,topo_id,decomp
      INTEGER                 , INTENT(  OUT) :: nsublist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(:,:), POINTER :: ineigh,subs_bc
      INTEGER , DIMENSION(  :), POINTER :: nneigh,nchld
      INTEGER , DIMENSION(3,1)          :: nnodes
      INTEGER , DIMENSION(3)            :: ldc
      REAL(MK), DIMENSION(3,2)          :: weights
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(:,:), POINTER :: min_box,max_box
      LOGICAL , DIMENSION(ppm_dim)      :: fixed
      INTEGER                           :: i,j,k,Ntot,iopt,treetype
      INTEGER                           :: istat,topoid,nbox,isub
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps
      INTEGER, DIMENSION(ppm_dim)       :: Nm
      CHARACTER(ppm_char)               :: mesg
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_mkpart',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mkpart',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (Npart .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &            'Npart must not be negative',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (.NOT.ASSOCIATED(xp)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &            'xp is not allocated',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ASSOCIATED(xp)) THEN
             IF(SIZE(xp,2) .LT. Npart) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &            'not enough particles contained in xp',&
     &            __LINE__, info)
                GOTO 9999
             ENDIF
             IF(SIZE(xp,1) .LT. ppm_dim) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &            'leading dimension of xp too small',&
     &            __LINE__, info)
                GOTO 9999
             ENDIF
          ENDIF
          IF (PRESENT(pcost)) THEN
              IF (SIZE(pcost,1) .LT. Npart) THEN
              info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &                'pcost must be of at least length Npart',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDIF
          IF ((ghostsize.LT.0.0_MK).AND.(decomp.NE.ppm_param_decomp_null)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &            'ghostsize must not be negative',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,ppm_dim
              IF (max_phys(i) .LE. min_phys(i)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &                'max_phys must be > min_phys',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
          IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (decomp .NE. ppm_param_decomp_user_defined) THEN
               info = ppm_error_warning
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &             'decomp type is set to user_defined for this assignment',&
     &             __LINE__, info)
               decomp = ppm_param_decomp_user_defined
            ENDIF
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 9999
            ENDIF
            DO i=1,nsubs
                IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 9999
                ENDIF
            ENDDO
          ENDIF
          IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
            IF(nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
     &             'no subs defined in user_defined decomposition',&
     &             __LINE__, info)
               GOTO 9999
            ENDIF
            IF ((.NOT.ASSOCIATED(min_sub)).OR.(.NOT.ASSOCIATED(max_sub))) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart', &
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
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',   &
     &              'faulty subdomains defined',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
         IF ((topo_id .EQ. 0) .AND. (decomp .NE. ppm_param_decomp_null)) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_argument, 'ppm_topo_mkpart',  &
     &           'topo_id 0 is reserved for Null decomposition',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      IF ((decomp .EQ. ppm_param_decomp_null) .AND. (topo_id .NE. 0)) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_argument, 'ppm_topo_mkpart',    &
     &     'topo_id was set to zero for null decomposition',__LINE__, info)
         topo_id = 0
      ENDIF
     
      IF ((decomp .NE. ppm_param_decomp_null) .AND. (topo_id .EQ. 0)) THEN
         info = ppm_error_warning
         CALL ppm_error(ppm_err_argument, 'ppm_topo_mkpart',    &
     &     'topo_id was reset for non-null decomposition',__LINE__, info)
         topo_id = -1
      ENDIF
     
      !-------------------------------------------------------------------------
      !  Check that we have particles
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_AllReduce(Npart,Ntot,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_mpi_fail,'ppm_topo_mkpart',   &
     &        'MPI_AllReduce failed',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (Ntot.LT.1) THEN
         info = ppm_error_notice
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',   &
     &       'No particles in domain',__LINE__,info)
         GOTO 9999
      ENDIF 
#else
      IF (Npart.LT.1) THEN
         info = ppm_error_notice
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',   &
     &       'No particle on this processor',__LINE__,info)
         GOTO 9999
      ENDIF 
#endif

      !----------------------------------------------------------------------
      !  Dummy argument for non-existing mesh
      !----------------------------------------------------------------------
      nnodes(1:3,1) = 0

      !----------------------------------------------------------------------
      !  Perform the decomposition using various techniques
      !----------------------------------------------------------------------
      IF     (decomp.EQ.ppm_param_decomp_pruned_cell) THEN
         !-------------------------------------------------------------------
         !  a pruned cell index list
         !-------------------------------------------------------------------
         IF (PRESENT(pcost)) THEN
             CALL ppm_decomp_pruned_cell(xp,Npart,min_phys,max_phys, &
     &            ghostsize,min_sub,max_sub,nsubs,info,pcost)
         ELSE
             CALL ppm_decomp_pruned_cell(xp,Npart,min_phys,max_phys, &
     &            ghostsize,min_sub,max_sub,nsubs,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Pruned cell decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
      ELSEIF (decomp.EQ.ppm_param_decomp_tree) THEN
         !-------------------------------------------------------------------
         !  a tree data structure; use a default maxvariance of 10 pct 
         !-------------------------------------------------------------------
         IF (PRESENT(pcost)) THEN
             CALL ppm_decomp_tree(xp,Npart,min_phys,max_phys,ghostsize, &
     &                        0.1_MK,min_sub,max_sub,nsubs,info,pcost)
         ELSE
             CALL ppm_decomp_tree(xp,Npart,min_phys,max_phys,ghostsize, &
     &                        0.1_MK,min_sub,max_sub,nsubs,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Tree decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
      ELSEIF (decomp.EQ.ppm_param_decomp_bisection) THEN
         !-------------------------------------------------------------------
         !  recursive bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! particles have unit weight
         weights(1,1)   = 1.0_MK
         weights(1,2)   = 1.0_MK
         ! mesh and geometry are not considered
         weights(2,1)   = 0.0_MK
         weights(2,2)   = 0.0_MK
         weights(3,1)   = 0.0_MK
         weights(3,2)   = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! no mesh
         Nm(1:ppm_dim)    = 0
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Bisection decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      ELSEIF ((decomp .EQ. ppm_param_decomp_xpencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_ypencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_zpencil)) THEN
         IF (decomp.EQ.ppm_param_decomp_zpencil.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &           'Cannot make z pencils in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  pencil quadrisection using the general ppm_tree
         !-------------------------------------------------------------------
         !  *** HACK FOR ZPENCILS WITH OLD TREE ***
         CALL ppm_decomp_tree(xp,Npart,min_phys,max_phys,2*ghostsize,0.1_MK, &
     &       min_sub,max_sub,nsubs,info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &           'Old tree hack fucked up!',__LINE__,info)
             GOTO 9999
         ENDIF
         GOTO 333
         !  ***************************************
         ! build a quad tree, binary in 2d
         treetype         = ppm_param_tree_quad
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
         ! particles have unit weight
         weights(1,1:2)   = 1.0_MK
         ! mesh and geometry are not considered
         weights(2,1:2)   = 0.0_MK
         weights(3,1:2)   = 0.0_MK
         ! fix the proper direction
         fixed(1:ppm_dim) = .FALSE.
         IF (decomp .EQ. ppm_param_decomp_xpencil) fixed(1) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_ypencil) fixed(2) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_zpencil) fixed(3) = .TRUE.
         gsvec(1:ppm_dim) = ghostsize
         ! no mesh
         Nm(1:ppm_dim)    = 0
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Pencil decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999
333   CONTINUE

      ELSEIF ((decomp .EQ. ppm_param_decomp_xy_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_xz_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_yz_slab)) THEN
         IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &           'Cannot make x-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',  &
     &           'Cannot make y-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! particles have unit weight
         weights(1,1:2)   = 1.0_MK
         ! mesh and geometry are not considered
         weights(2,1:2)   = 0.0_MK
         weights(3,1:2)   = 0.0_MK
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
         ! no mesh
         Nm(1:ppm_dim)    = 0
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Slab decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      ELSEIF (decomp .EQ. ppm_param_decomp_cuboid) THEN
         !-------------------------------------------------------------------
         !  cuboid octasection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build an oct tree in 3d
         treetype         = ppm_param_tree_oct
         ! and a quad tree in 2d
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_quad
         ! particles have unit weight
         weights(1,1:2)   = 1.0_MK
         ! mesh and geometry are not considered
         weights(2,1:2)   = 0.0_MK
         weights(3,1:2)   = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! no mesh
         Nm(1:ppm_dim)    = 0
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Cuboid decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      ELSEIF (decomp .EQ. ppm_param_decomp_user_defined) THEN
         ! Do nothing. Just take the stuff from the user and trust the guy.
      ELSEIF (decomp .EQ. ppm_param_decomp_null) THEN
         !-------------------------------------------------------------------
         !  null decomposition
         !-------------------------------------------------------------------
         CALL ppm_decomp_null(min_phys,max_phys,min_sub,max_sub,nsubs,info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed, 'ppm_topo_mkpart',  &
     &          'Null decomposition failed', __LINE__, info)
             GOTO 9999
         ENDIF
      ELSE
         !-------------------------------------------------------------------
         !  unknown decomposition
         !-------------------------------------------------------------------
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown decomposition type: ',decomp
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',   &
     &       mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the neighbours to the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, & 
     &                    min_sub,max_sub,nsubs,nneigh,ineigh,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &       'Finding neighbors failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          IF (PRESENT(pcost)) THEN
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,nnodes,cost, &
                                  info,pcost)
          ELSE
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,nnodes,cost, &
     &            info)
          ENDIF
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &           'Computing costs failed',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assignments for non-trivial decompositions
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_null) THEN
         !----------------------------------------------------------------------
         !  Define the topology (assign the subdomains to processors)
         !----------------------------------------------------------------------
         IF     (assig .EQ. ppm_param_assign_internal) THEN
            !-------------------------------------------------------------------
            !  internal assignment routine
            !-------------------------------------------------------------------
            CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, &
     &          isublist,nsublist,info)
            IF (info.NE.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &            'Assigning subs to processors failed',__LINE__,info)
               GOTO 9999
            ENDIF
         ELSEIF (assig .EQ. ppm_param_assign_nodal_cut .OR.      &
     &           assig .EQ. ppm_param_assign_nodal_comm .OR.     &
     &           assig .EQ. ppm_param_assign_dual_cut .OR.       &
     &           assig .EQ. ppm_param_assign_dual_comm) THEN
            !-------------------------------------------------------------------
            !  use METIS library to do assignment
            !-------------------------------------------------------------------
            CALL ppm_topo_metis_s2p(min_sub,max_sub,nneigh,ineigh,cost,nsubs,&
     &          assig,sub2proc,isublist,nsublist,info)
            IF (info.NE.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &            'Assigning subs to processors using METIS failed',__LINE__,&
     &            info)
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
                CALL ppm_error(ppm_err_alloc,'ppm_topo_mkpart',   &
     &              'list of local subs ISUBLIST',__LINE__,info)
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
            CALL ppm_error(ppm_err_argument,'ppm_topo_mkpart',   &
     &          mesg,__LINE__,info)
            GOTO 9999
         ENDIF

      !-------------------------------------------------------------------------
      !  Assignment for trivial (null) decompositions. This could also be
      !  put in a subroutine...
      !-------------------------------------------------------------------------
      ELSE
          nsublist = nsubs
          !---------------------------------------------------------------------
          !  Now the assignment of this one subdomain to this processor 
          !  (i.e. every processor has the same subdomain)
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldc(1) = nsubs
          CALL ppm_alloc(sub2proc,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_mkpart',  &
     &            'alloc of sub2proc failed',__LINE__,info)
              GOTO 9999
          ENDIF
          sub2proc(1) = ppm_rank
          !---------------------------------------------------------------------
          !  Fill the list of the subs that are on this processor
          !---------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          ldc(1) = nsublist
          CALL ppm_alloc(isublist,ldc,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_topo_mkpart',  &
     &            'alloc of isublist failed',__LINE__,info)
              GOTO 9999
          ENDIF
          isublist(1) = 1
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

      ! isublist:  my subs (ppm_rank)
      ! sub2proc:  subs -> processor map
      ! ineigh:    neighbors of all subs
      !-------------------------------------------------------------------------
      !  Store the current topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                    sub2proc,nsubs,bcdef,isublist,nsublist,    &
     &                    nneigh,ineigh,topo_id,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
     &       'Storing topology failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Get internal topoid. If no topology has ever been defined, set 
      !  this as the current topology.
      !-------------------------------------------------------------------------
      topoid = ppm_internal_topoid(topo_id)
      IF (ppm_topoid .EQ. -1) ppm_topoid = topoid

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_mkpart',t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkpart_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkpart_d
#endif
