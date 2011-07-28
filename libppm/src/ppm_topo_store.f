      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_store
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine stores all relevant information about
      !                 the new topology (created by ppm_decomp routines) in 
      !                 the global arrays and creates an internal topology ID 
      !                 for it (and updates the maping lists between user 
      !                 topoid and internal ppm topoid).
      !
      !  Input        : min_phys(:)  (F) the min. extent of the comput. domain
      !                 max_phys(:)  (F) the max. extent of the comput. domain
      !                 min_sub(:,:) (F) the min. extent of the subdomains
      !                 max_sub(:,:) (F) the max. extent of the subdomains
      !                 subs_bc(:,:) (I) BC for subs on ALL processors
      !                 nneigh(:)    (I) #neighbors of each sub
      !                 ineigh(:,:)  (I) neighbors of all subs
      !                 nsubs        (I) total number of subs on all procs
      !                 bcdef(:)     (I) BC on the computational box
      !                 nsublist     (I) number of subs on current processor
      !                 sub2proc(:)  (I) assignment of subs to procs
      !                 isublist(:)  (I) list of subs handled by this
      !                                  processor
      !
      !  Input/output : topo_id      (I) user (not ppm internal!) topology
      !                                  ID. If .LT. 0 on input, the
      !                                  routine will create an automatic
      !                                  one and return it here.
      !
      !  Output       : info         (I) return status. 0 upon success.
      !
      !  Remarks      : This routine dellocates the array: subs_bc passed to it.
      !                 The choice of placing it here is somewhat arbitrary, but
      !                 since this routine is call BOTH in ppm_topo_mkpart and
      !                 ppm_topo_mkfield - why not do it in one place - here.
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_store.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.34  2006/02/03 09:34:01  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.33  2004/08/31 10:24:34  ivos
      !  Several allocations in this routine failed when a processor had
      !  nsublist.EQ.0. This is fixed now.
      !
      !  Revision 1.32  2004/07/26 07:42:36  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.31  2004/04/22 15:15:39  walther
      !  Changed to capital letter for LBOUND and bug fixed the error message for
      !  ppm_min/max_phys allocation.
      !
      !  Revision 1.30  2004/04/22 10:47:30  oingo
      !  The lbound of the part of the arrays that have the topoid as index
      !  is now only zero if the zero-topology is or has been used.
      !
      !  Revision 1.29  2004/04/20 14:03:03  walther
      !  Added an error check after the call to the invert list routines and
      !  bug fixed the storage of the ppm_max_phys.
      !
      !  Revision 1.28  2004/04/20 11:08:07  walther
      !  Added min_phys and max_phys to the argument list and now stores these
      !  values for each the topoid.
      !
      !  Revision 1.27  2004/04/14 15:05:07  ivos
      !  corrected comment header and fixed alloc of ppm_ineighsubs (lower
      !  bounds are not needed and in fact lead to problems if all subs have
      !  no neighbors. They were removed.)
      !
      !  Revision 1.26  2004/04/13 15:11:59  oingo
      !  Changed the lower bounds to zero of the arrays that hold the topoid
      !
      !  Revision 1.25  2004/04/07 09:22:26  ivos
      !  Removed debug PRINT statements.
      !
      !  Revision 1.24  2004/04/05 14:36:33  walther
      !  Bug fix: had removed an iopt too much before allocating the 
      !  ppm_isublist.
      !
      !  Revision 1.23  2004/04/05 13:11:26  walther
      !  Bug fix: forgot to allocate ppm_subs_bc.
      !
      !  Revision 1.22  2004/04/05 12:39:55  walther
      !  Added the subs_bc.
      !
      !  Revision 1.21  2004/04/02 10:38:10  walther
      !  Now passing and storing the bcdef for each topology.
      !
      !  Revision 1.20  2004/04/01 15:58:02  walther
      !  How storing the neighbouring subs of the subs belonging to the local
      !  processor (in ppm_nneighsubs() and ppm_ineighsubs()).
      !
      !  Revision 1.19  2004/04/01 07:49:42  walther
      !  Added the nsubplus.
      !
      !  Revision 1.18  2004/02/06 16:39:23  ivos
      !  Removed unused mesg variable.
      !
      !  Revision 1.17  2004/02/06 13:52:02  ivos
      !  Updated argument checking for subs.
      !
      !  Revision 1.16  2004/02/04 17:24:01  ivos
      !  Header update and cosmetics.
      !
      !  Revision 1.15  2004/01/27 12:53:48  ivos
      !  Relaxed argument check from nsublist.LE.0 to nsublist.LT.0.
      !
      !  Revision 1.14  2004/01/23 17:22:12  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.13  2004/01/08 14:21:00  ivos
      !  Updated header comments.
      !
      !  Revision 1.12  2003/12/16 13:39:33  ivos
      !  Changed ALLOCATED into ASSOCIATED to make it compile with ifc AND pgf90.
      !
      !  Revision 1.11  2003/12/16 12:26:18  ivos
      !  Made ppm_nsublist and ppm_isublist depend on topoid.
      !
      !  Revision 1.10  2003/12/12 18:03:00  ivos
      !  Fixed 2 syntax errors.
      !
      !  Revision 1.9  2003/12/12 17:45:16  ivos
      !  Removed ALLOCs for two eliminated global topology arrays.
      !
      !  Revision 1.8  2003/12/12 16:13:53  ivos
      !  Now translates user topology IDs to internal ones and (re)allocates all
      !  topology lists (and of course updates them). Can now auto-generate 
      !  topology IDs and return them to the user.
      !
      !  Revision 1.7  2003/12/11 13:56:53  ivos
      !  Bugfix: index bug in loop that determines neighboring processors removed
      !  (replaced nneigh(ppm_isublist(i)) with nneigh(i) in the 2nd DO loop (j).
      !
      !  Revision 1.6  2003/12/05 14:41:04  ivos
      !  ppm_isoptimized is now allocated and initialized here.
      !
      !  Revision 1.5  2003/12/05 11:19:05  ivos
      !  Bugfix. It now compiles.
      !
      !  Revision 1.4  2003/12/05 10:06:55  ivos
      !  Added a few comments...
      !
      !  Revision 1.3  2003/12/04 17:48:13  ivos
      !  Bugfix: ppm_subs2proc(:,topoid). topoid inserted
      !
      !  Revision 1.2  2003/12/04 17:45:24  ivos
      !  Added memory allocation of neighbor processor list ppm_ineighlist. Added
      !  loops to find my neighbors and fill in these global lists.
      !
      !  Revision 1.1.1.1  2003/11/17 15:13:45  walther
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_store_s(min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                            sub2proc,nsubs,bcdef,&
     &                            isublist,nsublist,nneigh,ineigh,topo_id,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_store_d(min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                            sub2proc,nsubs,bcdef,&
     &                            isublist,nsublist,nneigh,ineigh,topo_id,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_util_invert_list
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys,max_phys
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub,max_sub
      INTEGER , DIMENSION(:,:), POINTER       :: subs_bc
      INTEGER , DIMENSION(  :), INTENT(IN   ) :: sub2proc,isublist,bcdef
      INTEGER                 , INTENT(IN   ) :: nsubs,nsublist
      INTEGER                 , INTENT(INOUT) :: topo_id
      INTEGER                 , INTENT(  OUT) :: info
      ! number of neighbors of each sub
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: nneigh
      ! neighbors of each sub. 1st index: 1...nneigh (neighbor index of
      ! sub), 2nd index: 1..nsubs (sub of which one wants to get the
      ! neighbors)
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ineigh
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(3) :: ldc, ldl
      INTEGER                :: i,j,k,kk,iopt,isize,iproc,isin,topoid
      INTEGER                :: maxneigh,minbound,nsubmax,nsublistmax
      REAL(MK)               :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_store',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (nsubs .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'nsubs must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsublist .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'nsublist must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (nsubs .LT. nsublist) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &            'total number of subs is smaller than local one',__LINE__,info)
              GOTO 9999
          ENDIF
          DO i=1,nsubs
              DO j=1,ppm_dim
                  IF(max_sub(j,i) .LE. min_sub(j,i)) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &                    'min_sub must be < max_sub',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ENDDO
          ENDDO
          DO i=1,nsubs
              IF(sub2proc(i) .LT. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_topo_store',  &
     &                'found sub not assigned to any processor',__LINE__,info)
                  GOTO 9999
              ENDIF
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Create a new internal topology number if the user specified none
      !-------------------------------------------------------------------------
      ! if the user specified topo_id .LT. 0, auto-create one
      IF (topo_id .LT. 0) THEN
          IF (ASSOCIATED(ppm_user_topoid)) THEN
              ! find the largest user ID so far...
              kk = 0
              DO i=1,ppm_max_topoid
                  IF (ppm_user_topoid(i) .GT. kk) kk = ppm_user_topoid(i)
              ENDDO
              ! ...and increment it by 1
              topo_id = kk+1
          ELSE
              ! otherwise start at the first one
              topo_id = 1
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we are creating a new topology or altering an existing one
      !-------------------------------------------------------------------------
      ! By default assume that it is a new one
      IF (topo_id .NE. 0) THEN
          topoid = ppm_max_topoid+1
          DO i=1,ppm_max_topoid
             IF (ppm_user_topoid(i) .EQ. topo_id) topoid = i
          END DO
      ELSE
          topoid = 0
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Increase the internal topology counter if needed
      !-------------------------------------------------------------------------
      IF (topoid .GT. ppm_max_topoid) ppm_max_topoid = topoid

      !-------------------------------------------------------------------------
      !  Set the lower bound of the arrays that have the topoid as index
      !-------------------------------------------------------------------------
      IF (ASSOCIATED(ppm_user_topoid)) THEN
          minbound = LBOUND(ppm_user_topoid,1)
      ELSE
          minbound = 1
      ENDIF
      
      IF (topoid .EQ. 0) THEN
          minbound = 0
      ENDIF

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the internal topology list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = minbound
      ldc(1) = ppm_max_topoid
      CALL ppm_alloc(ppm_user_topoid,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'user topoid list PPM_USER_TOPOID',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the user-provided ID number of this topology
      !-------------------------------------------------------------------------
      ppm_user_topoid(topoid) = topo_id

      !-------------------------------------------------------------------------
      !  Update the inverse list
      !-------------------------------------------------------------------------
      CALL ppm_util_invert_list(ppm_user_topoid,ppm_internal_topoid,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_store', &
     &                  'ppm_util_invert_list failed',__LINE__,info)
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the internally stored min_phys and max_phys
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = ppm_dim
      ldc(2) = ppm_max_topoid
#if   __KIND == __SINGLE_PRECISION
      !-------------------------------------------------------------------------
      !  Single precision
      !-------------------------------------------------------------------------
      CALL ppm_alloc(ppm_min_physs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'min extent of phys PPM_MIN_PHYSS',__LINE__,info)
          GOTO 9999
      ENDIF

      CALL ppm_alloc(ppm_max_physs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'max extent of phys PPM_MAX_PHYSS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  And store the values
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         ppm_min_physs(k,topoid) = min_phys(k)
         ppm_max_physs(k,topoid) = max_phys(k)
      ENDDO
#else
      !-------------------------------------------------------------------------
      !  Double precision
      !-------------------------------------------------------------------------
      CALL ppm_alloc(ppm_min_physd,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'min extent of phys PPM_MIN_PHYSD',__LINE__,info)
          GOTO 9999
      ENDIF

      CALL ppm_alloc(ppm_max_physd,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'max extent of phys PPM_MAX_PHYSD',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  And store the values
      !-------------------------------------------------------------------------
      DO k=1,ppm_dim
         ppm_min_physd(k,topoid) = min_phys(k)
         ppm_max_physd(k,topoid) = max_phys(k)
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  (Re)allocate memory for the internally stored subdomains 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = 1
      ldl(3) = minbound
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      ldc(3) = ppm_max_topoid
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(ppm_min_subs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'min extent of subs PPM_MIN_SUBS',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_max_subs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'max extent of subs PPM_MAX_SUBS',__LINE__,info)
          GOTO 9999
      ENDIF
#else
      CALL ppm_alloc(ppm_min_subd,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'min extent of subs PPM_MIN_SUBD',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_max_subd,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'max extent of subs PPM_MAX_SUBD',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Allocate memory for internally stored total number of subdomains 
      !-------------------------------------------------------------------------
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = nsubs
      ldc(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_subs2proc,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'global subs to proc map PPM_SUBS2PROC',__LINE__,info)
          GOTO 9999
      ENDIF
      DO i=1,nsubs
         ppm_subs2proc(i,topoid) = sub2proc(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for internally stored total number of subdomains 
      !-------------------------------------------------------------------------
      ldl(1) = minbound
      ldc(1) = ppm_max_topoid
      ! global
      CALL ppm_alloc(ppm_nsubs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'total number of subs PPM_NSUBS',__LINE__,info)
          GOTO 9999
      ENDIF
      ! local on this processor
      CALL ppm_alloc(ppm_nsublist,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'local number of subs on processor PPM_NSUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate and initialize optmization status flags
      !-------------------------------------------------------------------------
      CALL ppm_alloc(ppm_isoptimized,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'comm.opt. flags PPM_ISOPTIMIZED',__LINE__,info)
          GOTO 9999
      ENDIF
      ppm_isoptimized(topoid) = .FALSE.

      !-------------------------------------------------------------------------
      !  Store the number of subdomains in the current topology
      !-------------------------------------------------------------------------
      ppm_nsubs(topoid) = nsubs

      !-------------------------------------------------------------------------
      !  Store the extent of the subdomains
      !-------------------------------------------------------------------------
      IF (ppm_dim.EQ.2) THEN
         !----------------------------------------------------------------------
         !  In two dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
#if   __KIND == __SINGLE_PRECISION
            ppm_min_subs(1,i,topoid) = min_sub(1,i)
            ppm_min_subs(2,i,topoid) = min_sub(2,i)
            ppm_max_subs(1,i,topoid) = max_sub(1,i)
            ppm_max_subs(2,i,topoid) = max_sub(2,i)
#else
            ppm_min_subd(1,i,topoid) = min_sub(1,i)
            ppm_min_subd(2,i,topoid) = min_sub(2,i)
            ppm_max_subd(1,i,topoid) = max_sub(1,i)
            ppm_max_subd(2,i,topoid) = max_sub(2,i)
#endif
         ENDDO
      ELSE
         !----------------------------------------------------------------------
         !  In three dimensions
         !----------------------------------------------------------------------
         DO i=1,nsubs
#if   __KIND == __SINGLE_PRECISION
            ppm_min_subs(1,i,topoid) = min_sub(1,i)
            ppm_min_subs(2,i,topoid) = min_sub(2,i)
            ppm_min_subs(3,i,topoid) = min_sub(3,i)
            ppm_max_subs(1,i,topoid) = max_sub(1,i)
            ppm_max_subs(2,i,topoid) = max_sub(2,i)
            ppm_max_subs(3,i,topoid) = max_sub(3,i)
#else
            ppm_min_subd(1,i,topoid) = min_sub(1,i)
            ppm_min_subd(2,i,topoid) = min_sub(2,i)
            ppm_min_subd(3,i,topoid) = min_sub(3,i)
            ppm_max_subd(1,i,topoid) = max_sub(1,i)
            ppm_max_subd(2,i,topoid) = max_sub(2,i)
            ppm_max_subd(3,i,topoid) = max_sub(3,i)
#endif
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  Store the external boundary conditions for this topology
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = 2*ppm_dim
      ldc(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_bcdef,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'number of neighbors for each sub PPM_BCDEF',__LINE__,info)
          GOTO 9999
      ENDIF
      DO i=1,2*ppm_dim
         ppm_bcdef(i,topoid) = bcdef(i) 
      ENDDO
      
      !-------------------------------------------------------------------------
      !  The MAX of nsublist and 1 is needed to avoid allocation failures
      !  if a processor has 0 subs. Same for nsubs if there are no subs at
      !  all.
      !-------------------------------------------------------------------------
      nsubmax = MAX(nsubs,1)
      nsublistmax = MAX(nsublist,1)
      
      !-------------------------------------------------------------------------
      !  Allocate memory for the number of neighbours of the subs handled by
      !  the current processor i.e, a subset of the nneigh(1:nsubs) list.
      !  The required size is nsubmax and the current max. topoid
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = nsublistmax
      ldc(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_nneighsubs,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'number of neighbors for each sub PPM_NNEIGHSUBS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  ... and store the data by looping over the subs handled by the current
      !  processor: isublist(1:nsublist)
      !-------------------------------------------------------------------------
      maxneigh = 0
      IF (nsublist .LT. 1) THEN
          ! dummy entry if we have no subs at all
          ppm_nneighsubs(1,topoid) = ppm_param_undefined
      ENDIF
      DO i=1,nsublist
         j = nneigh(isublist(i))
         ppm_nneighsubs(i,topoid) = j
         IF (j.GT.maxneigh) maxneigh = j
      ENDDO

      !-------------------------------------------------------------------------
      !  Next allocate memory for the ID of these neighbouring subs. 
      !  The required size is maxneigh, nsublist, and the maximum topoid
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldc(1) = maxneigh
      ldc(2) = nsublist 
      ldc(3) = ppm_max_topoid
      CALL ppm_alloc(ppm_ineighsubs,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'number of neighbors for each sub PPM_INEIGHSUBS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  ... and store the data by looping over the subs handled by the current
      !  processor: isublist(1:nsublist) - the neighbours are the accessed 
      !  through the ineigh(1:nsublist,)
      !-------------------------------------------------------------------------
      DO i=1,nsublist
         DO j=1,ppm_nneighsubs(i,topoid)
            ppm_ineighsubs(j,i,topoid) = ineigh(j,isublist(i))
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for BC of all the subs
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = 1
      ldl(3) = minbound
      ldc(1) = 2*ppm_dim
      ldc(2) = nsubmax
      ldc(3) = ppm_max_topoid
      CALL ppm_alloc(ppm_subs_bc,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'BCs for subs on local processor failed',__LINE__,info)
          GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  And store it ...
      !-------------------------------------------------------------------------
      IF (nsubs .LT. 1) THEN
         ! dummy entries if we have no subs at all
         ppm_subs_bc(1,1,topoid) = ppm_param_undefined
         ppm_subs_bc(2,1,topoid) = ppm_param_undefined
         ppm_subs_bc(3,1,topoid) = ppm_param_undefined
         ppm_subs_bc(4,1,topoid) = ppm_param_undefined
         IF (ppm_dim.EQ.3) THEN
            ppm_subs_bc(5,1,topoid) = ppm_param_undefined
            ppm_subs_bc(6,1,topoid) = ppm_param_undefined
         ENDIF
      ENDIF
      DO i=1,nsubs
         ppm_subs_bc(1,i,topoid) = subs_bc(1,i)
         ppm_subs_bc(2,i,topoid) = subs_bc(2,i)
         ppm_subs_bc(3,i,topoid) = subs_bc(3,i)
         ppm_subs_bc(4,i,topoid) = subs_bc(4,i)
         IF (ppm_dim.EQ.3) THEN
            ppm_subs_bc(5,i,topoid) = subs_bc(5,i)
            ppm_subs_bc(6,i,topoid) = subs_bc(6,i)
         ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  And now deallocate the local array: subs_bc
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_store',     &
     &        'BCs for subs on local processor failed',__LINE__,info)
          GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Allocate memory for internally stored number of subdomains handled by
      !  the current processor
      !-------------------------------------------------------------------------
      ppm_nsublist(topoid) = nsublist
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = nsublistmax
      ldc(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_isublist,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'list of local subs PPM_ISUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (nsublist .LT. 1) THEN
         ppm_isublist(1,topoid) = ppm_param_undefined
      ENDIF
      DO i=1,nsublist
         ppm_isublist(i,topoid) = isublist(i)
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for internally stored number of neighbor
      !  processors and their list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldl(1) = minbound
      ldc(1) = ppm_max_topoid
      CALL ppm_alloc(ppm_nneighlist,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'number of neighbors for each sub PPM_NNEIGHLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      ppm_nneighlist(topoid) = 0
      ldl(1) = 1
      ldl(2) = minbound
      ldc(1) = 26   ! just guessing :-)
      ldc(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_ineighlist,ldl,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_store',     &
     &        'sub neighbor list PPM_INEIGHLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      isize = size(ppm_ineighlist,1)

      !-------------------------------------------------------------------------
      !  Determine and store the neighbors of this processor
      !-------------------------------------------------------------------------
      DO k=1,ppm_nsublist(topoid)
          i = ppm_isublist(k,topoid)
          DO j=1,nneigh(i)
              iproc = ppm_subs2proc(ineigh(j,i),topoid)
              ! if it is not myself
              IF (iproc .NE. ppm_rank) THEN
                  isin = 0
                  DO kk=1,ppm_nneighlist(topoid)
                      IF (ppm_ineighlist(kk,topoid) .EQ. iproc) isin = 1
                  END DO
                  ! if not already in list
                  IF (isin .EQ. 0) THEN
                      ppm_nneighlist(topoid) = ppm_nneighlist(topoid) + 1  
                      IF (ppm_nneighlist(topoid) .GT. isize) THEN
                          ! kindly ask for more memory
                          isize = isize + 2
                          ldc(1) = isize
                          ldc(2) = ppm_max_topoid
                          CALL ppm_alloc(ppm_ineighlist,ldc,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_topo_store',&
     &                            'sub neighbor list PPM_INEIGHLIST',  &
     &                            __LINE__,info)
                              GOTO 9999
                          ENDIF
                      END IF
                      ! add iproc to the list of neighbors
                      ppm_ineighlist(ppm_nneighlist(topoid),topoid) = iproc
                  END IF
              END IF
          END DO
      END DO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_store',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_store_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_store_d
#endif
