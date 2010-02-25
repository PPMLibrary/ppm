      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_topo_metis_s2p
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine assigns the subdomains to the processors
      !                 using the the METIS grid partitioning library.
      !
      !  Input        : min_sub(:,:)(F) min. extent of all subdomains
      !                 max_sub(:,:)(F) max. extent of all subdomains
      !                 nneigh(:)   (I) the number of neighbours of a
      !                                 subdomain
      !                 ineigh(:,:) (I) pointers to these neighbours
      !                 cost(:)     (F) costs of all subdomains (i.e. sum
      !                                 of particle/mesh node costs in
      !                                 that subdomain. Used for
      !                                 weighting the assignment.
      !                 nsubs       (I) total number of subdomains
      !                 assig       (I) METIS assignment type. One of:
      !                                   ppm_param_assign_nodal_cut
      !                                   ppm_param_assign_nodal_comm
      !                                   ppm_param_assign_dual_cut
      !                                   ppm_param_assign_dual_comm
      !
      !  Input/output :                                            
      !
      !  Output       : sub2proc(:) (I) full list of processor 
      !                                 affiliation of subdomains
      !                 isublist(:) (I) list of subdomains assigned to
      !                                 the local processors
      !                 nsublist    (I) the number of subdomains assigned
      !                 info        (I) return status 
      !
      !  Routines     : METIS library functions
      !
      !  Remarks      : Uses ppm_myeps for real comparisons in order to
      !                 decide if two points are the same. Checks these
      !                 decisions using the ineigh list.
      !
      !                 The corner numbering is currently done using a full
      !                 search. This is inefficient and should maybe be
      !                 replaced with some more clever algorithm.
      !
      !                 Currently, we make no use of link (graph edges)
      !                 weight as well as communication weights (number of
      !                 ghost particles/mesh points that need to be sent
      !                 for each edgecut). There is definitely room for
      !                 improvement here to get even better load balancing.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_topo_metis_s2p.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2006/09/05 08:01:27  pchatela
      !  Proper scaling for REAL comparisons
      !  Added module_alloc to ppm_decomp_boxsplit
      !
      !  Revision 1.3  2006/09/04 18:34:55  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.2  2004/10/01 16:09:12  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.1  2004/07/26 14:12:22  ivos
      !  Renamed from the versions without topo to resolve name conflict
      !  with global ppm_subs2proc data array.
      !
      !  Revision 1.1  2004/07/26 08:55:30  ivos
      !  Renamed.
      !
      !  Revision 1.9  2004/07/26 07:42:37  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.8  2004/04/14 15:06:53  ivos
      !  Removed debug prints.
      !
      !  Revision 1.7  2004/03/04 14:43:12  ivos
      !  bigfix: loop for node numbering in 3D case started from sub 1 instead
      !  of 2. fixed. 3D case tested now.
      !
      !  Revision 1.6  2004/02/27 14:53:10  ivos
      !  bugfix: nsublist was not initialized without METIS library support. ifc
      !  would complain that it is INTENT(OUT) and never set...
      !
      !  Revision 1.5  2004/02/27 14:50:54  ivos
      !  bugfix: corrected wrong error code.
      !
      !  Revision 1.4  2004/02/26 14:47:31  ivos
      !  bugfix: number of on-processor subs was not counted correctly.
      !  Routine has been tested for all assign types on 2 and 3 processors with
      !  both symmetric and asymmetric processor speeds.
      !
      !  Revision 1.3  2004/02/25 17:27:37  ivos
      !  Major update: does not use METIS mesh decomposition directly any more
      !  but converts the mesh to a graph first. Added new assignment schemes
      !  nodal_cut, nodal_comm, dual_cut and dual_comm. Makes much more use of
      !  METIS flexibility and different objective functions.
      !
      !  Revision 1.2  2004/02/25 10:23:21  ivos
      !  bugfix: corner numbering was wrong. fixed and tested to work OK.
      !
      !  Revision 1.1  2004/02/12 17:08:02  ivos
      !  Initial implementation. Not tested yet.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_metis_s2p_s(min_sub,max_sub,nneigh,ineigh,    &
     &               cost,nsubs,assig,sub2proc,isublist,nsublist,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_metis_s2p_d(min_sub,max_sub,nneigh,ineigh,    &
     &               cost,nsubs,assig,sub2proc,isublist,nsublist,info)
#endif
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: cost
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: nneigh
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: ineigh
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc,isublist
      INTEGER                 , INTENT(IN   ) :: nsubs,assig
      INTEGER                 , INTENT(  OUT) :: nsublist
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(:  ), POINTER   :: elmnts,npart,nxadj,nadjncy
      INTEGER , DIMENSION(:  ), POINTER   :: vwgt,adjwgt,vsize,vote
      INTEGER , DIMENSION(:,:), POINTER   :: cornerno
      REAL(MK), DIMENSION(:,:,:), POINTER :: corner
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: tpwgt
      INTEGER , DIMENSION(3)              :: ldc
      INTEGER , DIMENSION(1)              :: maxvote
      INTEGER , DIMENSION(5)              :: metis_options
      INTEGER                             :: iopt,isub,points,etype
      INTEGER                             :: nnodes,numtype,nparts,edgecut
      INTEGER                             :: j,k,i
      INTEGER                             :: jsub,jj,ident,nvert,wgtflag,volume
      CHARACTER(ppm_char)                 :: mesg
      REAL(MK)                            :: t0,lmyeps
      REAL(ppm_kind_double)               :: minsp,maxsp,meansp
      LOGICAL                             :: lasymm
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_metis_s2p',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Error if METIS library support is not available
      !-------------------------------------------------------------------------
#ifndef __METIS
      info = ppm_error_error
      CALL ppm_error(ppm_err_nometis,'ppm_topo_metis_s2p',  &
     &    'PPM was compiled without Metis support',__LINE__,info)
      nsublist = 0
      GOTO 9999      
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (nsubs .LE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &           'nsubs must be > 0',__LINE__,info)
            GOTO 9999
         ENDIF
         IF (assig .NE. ppm_param_assign_nodal_cut .AND.  &
     &       assig .NE. ppm_param_assign_nodal_comm .AND. &
     &       assig .NE. ppm_param_assign_dual_cut   .AND. &
     &       assig .NE. ppm_param_assign_dual_comm) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &           'Invalid assignment type for METIS passed',__LINE__,info)
            GOTO 9999
         ENDIF
         DO isub=1,nsubs
             IF (nneigh(isub) .LT. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &               'nneigh must be >= 0 for all subs',__LINE__,info)
                GOTO 9999
             ENDIF
             IF (cost(isub) .LT. 0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &               'costs must be >= 0 for all subs',__LINE__,info)
                GOTO 9999
             ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the sub2proc array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs
      CALL ppm_alloc(sub2proc,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'allocating ',nsubs,' sub2procs failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we only have one processor, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .EQ. 1) THEN
         !----------------------------------------------------------------------
         !  Allocate the isublist array
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             WRITE (mesg,'(A,I10,A)') 'allocating ',nsubs,' isublist failed'
             CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &           mesg,__LINE__,info)
             GOTO 9999
         ENDIF
         !----------------------------------------------------------------------
         !  Assign all subdomains to the processor and return
         !----------------------------------------------------------------------
         nsublist = nsubs
         DO isub=1,nsubs
            isublist(isub) = isub
            sub2proc(isub) = ppm_rank
         ENDDO
         GOTO 9999
      ENDIF 

      !-------------------------------------------------------------------------
      !  Number of points per sub and element type
      !-------------------------------------------------------------------------
      points = 8
      etype  = 3   ! METIS type hexahedra
      IF (ppm_dim .EQ. 2) THEN
          points = 4
          etype  = 4   ! METIS type quadrilaterals
      ENDIF

      !-------------------------------------------------------------------------
      !  Corner subdomain arrays
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = points
      ldc(3) = nsubs
      CALL ppm_alloc(corner,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'subdomain corner list CORNER',__LINE__,info)
          GOTO 9999
      ENDIF
      ldc(1) = points
      ldc(2) = nsubs
      CALL ppm_alloc(cornerno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'subdomain corner numbers CORNERNO',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize corner numbers
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
          DO j=1,points
              cornerno(j,isub) = 0
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Store sub corners in METIS numbering order
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          DO isub=1,nsubs
              corner(1,1,isub) = min_sub(1,isub)
              corner(2,1,isub) = min_sub(2,isub)
              
              corner(1,2,isub) = min_sub(1,isub)
              corner(2,2,isub) = max_sub(2,isub)
             
              corner(1,3,isub) = max_sub(1,isub)
              corner(2,3,isub) = max_sub(2,isub)
            
              corner(1,4,isub) = max_sub(1,isub)
              corner(2,4,isub) = min_sub(2,isub)
          ENDDO
      ELSE
          DO isub=1,nsubs
              corner(1,1,isub) = min_sub(1,isub)
              corner(2,1,isub) = min_sub(2,isub)
              corner(3,1,isub) = max_sub(3,isub)
          
              corner(1,2,isub) = min_sub(1,isub)
              corner(2,2,isub) = max_sub(2,isub)
              corner(3,2,isub) = max_sub(3,isub)
         
              corner(1,3,isub) = max_sub(1,isub)
              corner(2,3,isub) = max_sub(2,isub)
              corner(3,3,isub) = max_sub(3,isub)
        
              corner(1,4,isub) = max_sub(1,isub)
              corner(2,4,isub) = min_sub(2,isub)
              corner(3,4,isub) = max_sub(3,isub)
       
              corner(1,5,isub) = min_sub(1,isub)
              corner(2,5,isub) = min_sub(2,isub)
              corner(3,5,isub) = min_sub(3,isub)
      
              corner(1,6,isub) = min_sub(1,isub)
              corner(2,6,isub) = max_sub(2,isub)
              corner(3,6,isub) = min_sub(3,isub)
     
              corner(1,7,isub) = max_sub(1,isub)
              corner(2,7,isub) = max_sub(2,isub)
              corner(3,7,isub) = min_sub(3,isub)
    
              corner(1,8,isub) = max_sub(1,isub)
              corner(2,8,isub) = min_sub(2,isub)
              corner(3,8,isub) = min_sub(3,isub)
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Number corners with same point in space having same number
      !-------------------------------------------------------------------------
      cornerno(1,1) = 1
      cornerno(2,1) = 2
      cornerno(3,1) = 3
      cornerno(4,1) = 4
      nnodes = 4
      IF (ppm_dim .EQ. 2) THEN
         DO isub=2,nsubs
            DO j=1,points
               ident = 0 
               !----------------------------------------------------------------
               !  Check all points of all neighbors
               !----------------------------------------------------------------
               DO k=1,nneigh(isub)
                  jsub = ineigh(k,isub)
                  DO jj=1,points
                     !----------------------------------------------------------
                     !  If points coincide, store point number
                     !----------------------------------------------------------
                     IF((ABS(corner(1,j,isub)-corner(1,jj,jsub)).LT.   &
     &                  lmyeps*(max_sub(1,isub)-min_sub(1,isub))).AND. &
     &                  (ABS(corner(2,j,isub)-corner(2,jj,jsub)).LT. &
     &                  lmyeps*(max_sub(2,isub)-min_sub(2,isub)))) THEN
                        IF (cornerno(jj,jsub) .GT. 0) ident = cornerno(jj,jsub)
                     ENDIF
                  ENDDO
               ENDDO
               !----------------------------------------------------------------
               !  If identical point was found: inherit number, else
               !  increase counter of distinct nodes and get new number
               !----------------------------------------------------------------
               IF (ident .GT. 0) THEN
                   cornerno(j,isub) = ident
               ELSE
                   nnodes = nnodes + 1
                   cornerno(j,isub) = nnodes
               ENDIF
            ENDDO
         ENDDO
      ELSE
         cornerno(5,1) = 5
         cornerno(6,1) = 6
         cornerno(7,1) = 7
         cornerno(8,1) = 8
         nnodes = 8
         DO isub=2,nsubs
            DO j=1,points
               ident = 0
               !----------------------------------------------------------------
               !  Check all points of all neighbors
               !----------------------------------------------------------------
               DO k=1,nneigh(isub)
                  jsub = ineigh(k,isub)
                  DO jj=1,points
                     !----------------------------------------------------------
                     !  If points coincide, store point number
                     !----------------------------------------------------------
                     IF((ABS(corner(1,j,isub)-corner(1,jj,jsub)).LT.   &
     &                  lmyeps*(max_sub(1,isub)-min_sub(1,isub))).AND. &
     &                  (ABS(corner(2,j,isub)-corner(2,jj,jsub)).LT.   &
     &                  lmyeps*(max_sub(2,isub)-min_sub(2,isub))).AND. &
     &                  (ABS(corner(3,j,isub)-corner(3,jj,jsub)).LT.   &
     &                  lmyeps*(max_sub(3,isub)-min_sub(3,isub)))) THEN
                        IF (cornerno(jj,jsub) .GT. 0) ident = cornerno(jj,jsub)
                     ENDIF
                  ENDDO
               ENDDO
               !----------------------------------------------------------------
               !  If identical point was found: inherit number, else
               !  increase counter of distinct nodes and get new number
               !----------------------------------------------------------------
               IF (ident .GT. 0) THEN
                   cornerno(j,isub) = ident
               ELSE
                   nnodes = nnodes + 1
                   cornerno(j,isub) = nnodes
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Sanity check
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (MAXVAL(cornerno) .NE. nnodes .OR. MINVAL(cornerno) .NE. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_node_number,'ppm_topo_metis_s2p', &
     &             'Nodes were missed or counted twice',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate corner list
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(corner,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',          &
     &        'subdomain corner list CORNER',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we need asymmetric partitions or not
      !-------------------------------------------------------------------------
      minsp = MINVAL(ppm_proc_speed(0:ppm_nproc-1))
      maxsp = MAXVAL(ppm_proc_speed(0:ppm_nproc-1))
      meansp = 1.0_ppm_kind_double/REAL(ppm_nproc,ppm_kind_double)
      meansp = (ABS(maxsp-minsp))/meansp
      lasymm = .FALSE. 
      ! if there is more than 5 percent difference, do it 
      IF (meansp .GT. 0.05_ppm_kind_double) lasymm = .TRUE.

      IF (ppm_debug .GT. 0) THEN
          WRITE(mesg,'(A,F6.4)') 'Slowest processor: ',minsp
          CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p',mesg,info)
          WRITE(mesg,'(A,F6.4)') 'Fastest processor: ',maxsp
          CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p',mesg,info)
          WRITE(mesg,'(A,F6.4)') 'Processor speed imbalance: ',meansp
          CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p',mesg,info)
          IF (lasymm) THEN
             CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p', &
     &          'Doing asymmetric assignment',info)
          ELSE
             CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p', &
     &          'Doing symmetric assignment',info)
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Number of vertices in the graph
      !-------------------------------------------------------------------------
      IF (assig .EQ. ppm_param_assign_nodal_cut .OR.   &
     &    assig .EQ. ppm_param_assign_nodal_comm) THEN
          nvert = nnodes
      ELSE
          nvert = nsubs
      ENDIF

      !-------------------------------------------------------------------------
      !  Number of desired graph partitions
      !-------------------------------------------------------------------------
      nparts  = ppm_nproc

      !-------------------------------------------------------------------------
      !  Allocate METIS mesh data array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs*points
      CALL ppm_alloc(elmnts,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS mesh data array ELMNTS',__LINE__,info)
          GOTO 9999
      ENDIF
      ldc(1) = nvert
      CALL ppm_alloc(npart,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS node partition vector NODEPART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(vwgt,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS vertex weight vector VWGT',__LINE__,info)
          GOTO 9999
      ENDIF
      ldc(1) = 1
      CALL ppm_alloc(adjwgt,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS vertex weight vector ADJWGT',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(vsize,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS vertex weight vector VSIZE',__LINE__,info)
          GOTO 9999
      ENDIF
  !    NULLIFY(adjwgt)   ! we make no use of link weights
  !    NULLIFY(vsize)    ! we make no use of communication weights
      ldc(1) = nvert+1
      CALL ppm_alloc(nxadj,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph structure array NXADJ',__LINE__,info)
          GOTO 9999
      ENDIF
      ldc(1) = 6*nvert
      IF (ppm_dim .EQ. 2) ldc(1) = 4*nvert
      CALL ppm_alloc(nadjncy,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph adjacency array NADJNCY',__LINE__,info)
          GOTO 9999
      ENDIF
      IF (lasymm) THEN
         ldc(1) = nparts
         CALL ppm_alloc(tpwgt,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',         &
     &           'METIS partition weights TPWGT',__LINE__,info)
             GOTO 9999
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill METIS mesh data structure
      !-------------------------------------------------------------------------
#ifdef __VECTOR
      !-------------------------------------------------------------------------
      !  For long vector length: inner loop is over subs
      !-------------------------------------------------------------------------
      DO j=1,points
          DO isub=1,nsubs
              elmnts((isub-1)*points+j) = cornerno(j,isub)
          ENDDO
      ENDDO
#else
      !-------------------------------------------------------------------------
      !  For memory stride: inner loop is over points
      !-------------------------------------------------------------------------
      DO isub=1,nsubs
          DO j=1,points
              elmnts((isub-1)*points+j) = cornerno(j,isub)
          ENDDO
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  Store partition weights for METIS if needed
      !-------------------------------------------------------------------------
      IF (lasymm) THEN
         DO i=1,nparts
            tpwgt(i) = REAL(ppm_proc_speed(i-1),ppm_kind_single)
         ENDDO 
      ENDIF

      !-------------------------------------------------------------------------
      !  Convert mesh to graph using METIS
      !-------------------------------------------------------------------------
      ! use Frotran numbering (starting from 1)
      numtype = 1
      IF (assig .EQ. ppm_param_assign_nodal_cut .OR.   &
     &    assig .EQ. ppm_param_assign_nodal_comm) THEN
         CALL METIS_MeshToNodal(nsubs,nnodes,elmnts,etype,numtype,nxadj,  &
     &      nadjncy) 
         wgtflag = 2  ! weights on vertices only / comput. weights only
         DO i=1,nvert
            vwgt(i)  = 0
         ENDDO
         DO i=1,nsubs
            DO j=1,points
               ! computation weights
               vwgt(cornerno(j,i)) = vwgt(cornerno(j,i)) + INT(cost(i))
            ENDDO
         ENDDO
      ELSEIF (assig .EQ. ppm_param_assign_dual_cut .OR.   &
     &        assig .EQ. ppm_param_assign_dual_comm) THEN
         CALL METIS_MeshToDual(nsubs,nnodes,elmnts,etype,numtype,nxadj,  &
     &      nadjncy)
         wgtflag = 2  ! weights on vertices only / comput. weights only
         DO i=1,nsubs
            ! computation weights
            vwgt(i)  = INT(cost(i))
         ENDDO
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &      'Unknown METIS assignment scheme. Bailing out.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Partition graph into nparts parts using METIS and the minimal
      !  edgecut objective
      !-------------------------------------------------------------------------
      IF (assig .EQ. ppm_param_assign_nodal_cut .OR.   &
     &    assig .EQ. ppm_param_assign_dual_cut) THEN
         IF (nparts .LE. 8) THEN
            !-------------------------------------------------------------------
            !  Use recursive algorithm since it produces better partitions
            !-------------------------------------------------------------------
            metis_options(1) = 1 ! use user-set values
            metis_options(2) = 3 ! sorted heavy-edge matching
            metis_options(3) = 1 ! region growing preprocessing
            metis_options(4) = 1 ! boundary refinement postprocessing
            metis_options(5) = 0 ! no debugging
            IF (lasymm) THEN
               CALL METIS_WPartGraphRecursive(nvert,nxadj,nadjncy,vwgt,adjwgt, &
     &            wgtflag,numtype,nparts,tpwgt,metis_options,edgecut,npart)
            ELSE
               CALL METIS_PartGraphRecursive(nvert,nxadj,nadjncy,vwgt,adjwgt, &
     &            wgtflag,numtype,nparts,metis_options,edgecut,npart)
            ENDIF
         ELSE
            !-------------------------------------------------------------------
            !  Use multilevel algorithm since it is faster
            !-------------------------------------------------------------------
            metis_options(1) = 1 ! use user-set values
            metis_options(2) = 3 ! sorted heavy-edge matching
            metis_options(3) = 1 ! initialize with recursive bisection
            metis_options(4) = 3 ! minimize connectivity of subdomains
            metis_options(5) = 0 ! no debugging
            IF (lasymm) THEN
               CALL METIS_WPartGraphKway(nvert,nxadj,nadjncy,vwgt,adjwgt,   &
     &            wgtflag,numtype,nparts,tpwgt,metis_options,edgecut,npart)
            ELSE
               CALL METIS_PartGraphKway(nvert,nxadj,nadjncy,vwgt,adjwgt,   &
     &            wgtflag,numtype,nparts,metis_options,edgecut,npart)
            ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Output diagnostics
         !----------------------------------------------------------------------
         ! egdecut is the number of mesh edges that has been cut by the
         ! partition. 
         IF (ppm_debug .GT. 0) THEN
            WRITE(mesg,'(A,I6)') 'METIS returned edgecut = ',edgecut
            CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p',mesg,info)
         ENDIF

      !-------------------------------------------------------------------------
      !  Partition graph into nparts parts using METIS and the minimal
      !  communication volume objective
      !-------------------------------------------------------------------------
      ELSEIF (assig .EQ. ppm_param_assign_nodal_comm .OR.   &
     &        assig .EQ. ppm_param_assign_dual_comm) THEN
         !----------------------------------------------------------------------
         !  For this objective METIS only has the multilevel algorithm
         !----------------------------------------------------------------------
         metis_options(1) = 1 ! use user-set values
         metis_options(2) = 3 ! sorted heavy-edge matching
         metis_options(3) = 1 ! initialize with recursive bisection
         metis_options(4) = 3 ! minimize connectivity of subdomains
         metis_options(5) = 0 ! no debugging
         IF (lasymm) THEN
            CALL METIS_WPartGraphVKway(nvert,nxadj,nadjncy,vwgt,vsize,   &
     &         wgtflag,numtype,nparts,tpwgt,metis_options,volume,npart)
         ELSE
            CALL METIS_PartGraphVKway(nvert,nxadj,nadjncy,vwgt,vsize,   &
     &         wgtflag,numtype,nparts,metis_options,volume,npart)
         ENDIF

         !----------------------------------------------------------------------
         !  Output diagnostics
         !----------------------------------------------------------------------
         ! egdecut is the number of mesh edges that has been cut by the
         ! partition. 
         IF (ppm_debug .GT. 0) THEN
            WRITE(mesg,'(A,I6)') 'METIS returned communication volume = ',volume
            CALL ppm_write(ppm_rank,'ppm_topo_metis_s2p',mesg,info)
         ENDIF
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &      'Unknown METIS assignment scheme. Bailing out.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find subs2proc assignment based on graph decomposition
      !-------------------------------------------------------------------------
      IF (assig .EQ. ppm_param_assign_nodal_cut .OR.   &
     &    assig .EQ. ppm_param_assign_nodal_comm) THEN
         ! a vertex is a mesh node: take a majority vote of all corners of
         ! a sub to determine its processor affiliation
         ldc(1) = nparts
         CALL ppm_alloc(vote,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',         &
     &           'majority vote count VOTE',__LINE__,info)
             GOTO 9999
         ENDIF
         DO i=1,nsubs
            vote(1:nparts) = 0
            DO j=1,points
               ! the processor this node was assigned to 
               jj = npart(cornerno(j,i))
               vote(jj) = vote(jj) + 1
            ENDDO
            ! majority vote of all corners of this sub
            maxvote = MAXLOC(vote)
            sub2proc(i) = maxvote(1)-1  ! MPI ranks start at 0
         ENDDO
      ELSEIF (assig .EQ. ppm_param_assign_dual_cut .OR.   &
     &        assig .EQ. ppm_param_assign_dual_comm) THEN
         ! a graph vertex is a sub
         DO i=1,nsubs
            sub2proc(i) = npart(i)-1    ! MPI ranks start at 0
         ENDDO
      ELSE
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_topo_metis_s2p', &
     &      'Unknown METIS assignment scheme. Bailing out.',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate METIS work memory
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      IF (lasymm) THEN
         CALL ppm_alloc(tpwgt,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',         &
     &           'METIS graph adjacency array TPWGT',__LINE__,info)
         ENDIF
      ENDIF
      CALL ppm_alloc(nadjncy,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph adjacency array NADJNCY',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nxadj,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph structure array NXADJ',__LINE__,info)
      ENDIF
      CALL ppm_alloc(vsize,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph structure array VSIZE',__LINE__,info)
      ENDIF
      CALL ppm_alloc(adjwgt,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS graph structure array ADJWGT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(vwgt,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS vertex weight vector VWGT',__LINE__,info)
      ENDIF
      CALL ppm_alloc(npart,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS node partition vector NODEPART',__LINE__,info)
      ENDIF
      CALL ppm_alloc(elmnts,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'METIS mesh data array ELMNTS',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Deallocate corner numbers
      !-------------------------------------------------------------------------
      iopt   = ppm_param_dealloc
      CALL ppm_alloc(cornerno,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_topo_metis_s2p',            &
     &        'subdomain corner numbers CORNERNO',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Count number of subs assigned to local processor
      !-------------------------------------------------------------------------
      nsublist = 0
      DO isub=1,nsubs
          IF (sub2proc(isub) .EQ. ppm_rank) nsublist = nsublist + 1
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate the isublist array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsublist
      CALL ppm_alloc(isublist,ldc,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          WRITE (mesg,'(A,I10,A)') 'allocating ',nsublist,' isublist failed'
          CALL ppm_error(ppm_err_alloc,'ppm_topo_metis_s2p',            &
     &                   mesg,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Fill list of local subs
      !-------------------------------------------------------------------------
      nsublist = 0
      DO isub=1,nsubs
          IF (sub2proc(isub) .EQ. ppm_rank) THEN
              nsublist = nsublist + 1
              isublist(nsublist) = isub
          ENDIF
      ENDDO
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_metis_s2p',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_metis_s2p_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_metis_s2p_d
#endif
