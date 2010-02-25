      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_util_commopt
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Determine an approximately optimal communication
      !                 sequence for each processor to SendRecv data with
      !                 its neighbors only. Such that: no conflicts occur
      !                 (i.e. A wants to send to B, but B is currently busy
      !                 receiving from C) and that a minimum number of
      !                 communication rounds are needed. This is done by
      !                 using the Vizing approximation to the minimal edge
      !                 coloring problem of the processor topology graph.
      !
      !  Input        : topoid     (I) topology ID to be optimized
      !                                (internal numbering)
      !
      !  Output       : info       (I) return status. 0 upon success.
      !
      !  Routines     : vizing_coloring (libvizing)
      !
      !  Remarks      : 
      !
      !  References   : V.G. Vizing, On an estimate of the chromatic class
      !                 of a p-graph. Discret. Analiz. 3, 1964, pp.25-30.
      !                 (In Russian).
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_commopt.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:19:01  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.16  2006/07/03 12:58:36  ivos
      !  indentation cosmetics.
      !
      !  Revision 1.15  2004/10/01 16:33:39  ivos
      !  cosmetics.
      !
      !  Revision 1.14  2004/10/01 16:09:13  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.13  2004/08/31 13:29:59  ivos
      !  changed argument checks to use ppm_check_topoid and ppm_check_meshid.
      !
      !  Revision 1.12  2004/07/26 07:42:32  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.11  2004/03/19 08:24:22  hiebers
      !  minor change at including files
      !
      !  Revision 1.10  2004/01/23 17:22:13  ivos
      !  Cleanup: (1) updated header, (2) inserted ppm_write and ppm_error, (3)
      !  inserted checks after every allocate, (4) added argument checking.
      !
      !  Revision 1.9  2003/12/17 14:00:57  ivos
      !  bug fix: size check for ilinks changed from nlinks to 2*nlinks.
      !
      !  Revision 1.8  2003/12/16 13:38:48  ivos
      !  bug fix: length of ilinks set to 2*nlinks (because each link has 2 
      !  nodes). Eliminated SIGSEGV in vizing_coloring.
      !
      !  Revision 1.7  2003/12/16 08:48:16  ivos
      !  Replaced grow_preserve with fit_preserve when shrinking lists before
      !  calling C++.
      !
      !  Revision 1.6  2003/12/11 14:01:05  ivos
      !  Bugfixes: (1) MPI bug: ranks .GT. 0 tried to access data that was not
      !  allocated on them. Fixed by placing an IF(rank.EQ.0) around it. (2)
      !  ppm_icommseq now always includes the processor itself as its first entry
      !  and ppm_ncommseq is increased by 1 as map_part_send needs it that way.
      !
      !  Revision 1.5  2003/12/09 11:27:52  hiebers
      !  merged
      !
      !  Revision 1.4  2003/12/09 11:18:41  ivos
      !  Bugfix: status is now in an #ifdef block for __MPI.
      !
      !  Revision 1.3  2003/12/09 09:35:39  ivos
      !  Changed INTENT of info to INOUT.
      !
      !  Revision 1.2  2003/12/09 08:56:50  ivos
      !  First complete version of the communication optimizer utility.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_util_commopt(topoid,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_check_topoid
      IMPLICIT NONE
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
      ! The ID of the topology which is to be optimized
      INTEGER                 , INTENT(IN   ) :: topoid
      ! return status
      INTEGER                 , INTENT(INOUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(2)                 :: ldu
      INTEGER                               :: iopt
      ! number of neighborhood relations in total
      INTEGER                               :: nlinks
      ! DISTINCT neighbor pairs, i.e. 1<->2 and 2<->1 is the same and only
      ! listed once. even indices are first points, odd ones second ones of
      ! the same edges.
      INTEGER, DIMENSION(:  ) , POINTER     :: ilinks
      ! optimal edge coloring determined. sequence of triples (p1,p2,c),...
      ! with p1 and p2 being the 2 vertices of each edge and c its color.
      INTEGER, DIMENSION(:  ) , POINTER     :: optres
      ! number of neighbors of all every CPU. index: MPI rank
      INTEGER, DIMENSION(:  ) , POINTER     :: nneighprocs
      ! all neighbors of all processors. 1st index: neighbor nr., 2nd:
      ! processor rank
      INTEGER, DIMENSION(:,:) , POINTER     :: ineighprocs
      INTEGER                               :: i,j,maxneigh,isize,ii,isin
      ! processor ranks
      INTEGER                               :: p1,p2
      ! min and max of assigned colors
      INTEGER                               :: mincolor,maxcolor
      LOGICAL                               :: valid
#ifdef __MPI
      ! MPI comm status
      INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_commopt',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_check_topoid(ppm_param_id_internal,topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_commopt',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if there are more than 1 processor. If not, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .LT. 2) THEN
          !---------------------------------------------------------------------
          !  Allocate memory for communication protocols
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_grow_preserve
          ldu(1) = 1
          ldu(2) = ppm_max_topoid
          CALL ppm_alloc(ppm_icommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'communication sequence PPM_ICOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF
          iopt   = ppm_param_alloc_grow_preserve
          ldu(1) = ppm_max_topoid
          CALL ppm_alloc(ppm_ncommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'number of comm.rounds PPM_NCOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Set the trivial protocol: only myself
          !---------------------------------------------------------------------
          ppm_ncommseq(topoid) = 1
          ppm_icommseq(1,topoid) = ppm_rank
          GOTO 9999
      END IF

#ifdef __MPI

      !-------------------------------------------------------------------------
      !  Allocate memory for number of neighbors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc
      CALL ppm_alloc(nneighprocs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'number of neighbors NNEIGHPROCS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Rank 0 receives the numbers of neighbors from all other processors
      !-------------------------------------------------------------------------
      maxneigh = 0
      IF (ppm_rank .GT. 0) THEN
          CALL MPI_Send(ppm_nneighlist(topoid),1,MPI_INTEGER,0,ppm_rank,   &
     &                  ppm_comm,info)
      ELSE
          nneighprocs(1) = ppm_nneighlist(topoid)
          maxneigh = nneighprocs(1)
          DO i=1,ppm_nproc-1
              CALL MPI_Recv(nneighprocs(i+1),1,MPI_INTEGER,i,i,   &
     &                      ppm_comm,status,info)
              IF (nneighprocs(i+1) .GT. maxneigh) maxneigh=nneighprocs(i+1)
          END DO
      END IF

      !-------------------------------------------------------------------------
      !  Allocate memory for neighbor lists
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = maxneigh
      ldu(2) = ppm_nproc
      CALL ppm_alloc(ineighprocs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'neighbor list INEIGHPROCS',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Rank 0 receives all neighbor lists
      !-------------------------------------------------------------------------
      IF (ppm_rank .GT. 0) THEN
          CALL MPI_Send(ppm_ineighlist(1:ppm_nneighlist(topoid),topoid), &
     &                  ppm_nneighlist(topoid),  &
     &                  MPI_INTEGER,0,ppm_rank,ppm_comm,info)
      ELSE
          ineighprocs(1:nneighprocs(1),1) =     &
     &                       ppm_ineighlist(1:ppm_nneighlist(topoid),topoid)
          DO i=1,ppm_nproc-1
              CALL MPI_Recv(ineighprocs(1:nneighprocs(i+1),i+1),   &
     &                      nneighprocs(i+1),MPI_INTEGER,i,i,ppm_comm,  &
     &                      status,info)
          END DO
      END IF

      !-------------------------------------------------------------------------
      !  Rank 0: Build graph and call optimizer
      !-------------------------------------------------------------------------
      IF (ppm_rank .EQ. 0) THEN
          !---------------------------------------------------------------------
          !  Build graph edges as UNIQUE neighbor pairs
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          isize  = 2*ppm_nproc   ! initial guess
          ldu(1) = isize
          CALL ppm_alloc(ilinks,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'links ILINKS',__LINE__,info)
              GOTO 9999
          ENDIF
          nlinks = 0
          DO i=1,ppm_nproc
              DO j=1,nneighprocs(i)
                  ! Check is the reverse link is already in the list
                  isin = 0
                  DO ii=1,nlinks
                      IF ((ilinks(2*ii-1) .EQ. (ineighprocs(j,i)+1)) .AND.  &
     &                    (ilinks(2*ii) .EQ. i)) isin = 1
                  END DO
                  ! add link if reverse is not already in the list
                  IF (isin .EQ. 0) THEN
                      nlinks = nlinks + 1
                      ! every link will cause two entries in ilinks !!
                      IF (2*nlinks .GT. isize) THEN
                          isize  = isize + 10
                          iopt   = ppm_param_alloc_grow_preserve
                          ldu(1) = isize
                          CALL ppm_alloc(ilinks,ldu,iopt,info)
                          IF (info .NE. 0) THEN
                              info = ppm_error_fatal
                              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',&
     &                            'links ILINKS',__LINE__,info)
                              GOTO 9999
                          ENDIF
                      END IF
                      ilinks(2*nlinks-1)  = i
                      ! in MPI numbering starts at 0, but here at 1 :-)
                      ilinks(2*nlinks)    = ineighprocs(j,i)+1
                  END IF
              END DO
          END DO

          !---------------------------------------------------------------------
          !  Allocate memory for optimal coloring result
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          ldu(1) = 3*nlinks
          CALL ppm_alloc(optres,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'optimal result OPTRES',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Shrink ilinks to the correct size. C does not like it
          !  otherwise.
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit_preserve
          ldu(1) = 2*nlinks
          CALL ppm_alloc(ilinks,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'final links ILINKS',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          ! Call the C++ wrapper to the Vizing code obtained from:
          ! The Stony Brook Algorithm Repository 
          ! http://www.cs.sunysb.edu/~algorith/implement/stony/distrib/Vizing/
          !---------------------------------------------------------------------
          CALL vizing_coloring(ppm_nproc,nlinks,ilinks,optres)
          !---------------------------------------------------------------------
          !  optres now contains the result as a sequence of nlinks triples
          !  (p1,p2,c) where each triple is one edge (determined by its two
          !  vertices p1 and p2) with a color c. Numbering of colors and
          !  vertices starts at 1.
          !---------------------------------------------------------------------

          !---------------------------------------------------------------------
          !  Determine smallest and largest assigned color
          !---------------------------------------------------------------------
          mincolor = HUGE(1)
          maxcolor = -HUGE(1)
          DO i=3,3*nlinks,3
              IF (optres(i) .LT. mincolor) mincolor = optres(i)
              IF (optres(i) .GT. maxcolor) maxcolor = optres(i)
          END DO
              
          !---------------------------------------------------------------------
          !  Allocate memory for neighbor communication lists
          !---------------------------------------------------------------------
          ii = maxcolor-mincolor+1
          iopt   = ppm_param_alloc_grow
          ldu(1) = ii
          ldu(2) = ppm_nproc
          CALL ppm_alloc(ineighprocs,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'neighbor comm. list INEIGHPROCS',__LINE__,info)
              GOTO 9999
          ENDIF
          ! initialize array. -1 means: no communication in this round
          ineighprocs(1:ii,1:ppm_nproc) = -1

          !---------------------------------------------------------------------
          !  Assemble communication protocol
          !---------------------------------------------------------------------
          ii = 0
          DO j=mincolor,maxcolor
              ii = ii + 1
              DO i=3,3*nlinks,3
                  ! if this link has the correct color
                  IF (optres(i) .EQ. j) THEN
                      p1 = optres(i-1)
                      p2 = optres(i-2)
                      ! MPI ranks start at 0.
                      ineighprocs(ii,p2) = p1-1
                      ineighprocs(ii,p1) = p2-1
                  END IF
              END DO
          END DO
      ENDIF     ! ppm_rank .EQ. 0

      !-------------------------------------------------------------------------
      !  Allocate memory for communication protocols
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = ppm_max_topoid
      CALL ppm_alloc(ppm_ncommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'number of comm. sequences PPM_NCOMMSEQ',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Distribute protocol to all processors
      !-------------------------------------------------------------------------
      ! First, broadcast the number of rounds to everybody
      IF (ppm_rank .EQ. 0) ppm_ncommseq(topoid) = ii
          CALL MPI_Bcast(ppm_ncommseq(topoid),1,MPI_INTEGER,0,ppm_comm,info)

      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ! the +1 is due to the fact that the processor itself also needs 
      ! to be in the list (as element 1) in order for map_part_send to 
      ! work properly
      ppm_ncommseq(topoid) = ppm_ncommseq(topoid) + 1
      ldu(1) = ppm_ncommseq(topoid)   
      ldu(2) = ppm_max_topoid
      CALL ppm_alloc(ppm_icommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'final communication sequence PPM_ICOMMSEQ',__LINE__,info)
          GOTO 9999
      ENDIF

      ! Then distribute the individual optimized neighbor lists
      IF (ppm_rank .GT. 0) THEN
          CALL MPI_Recv(ppm_icommseq(2:ppm_ncommseq(topoid),topoid), &
     &                  ppm_ncommseq(topoid)-1,  &
     &                  MPI_INTEGER,0,ppm_rank,ppm_comm,status,info)
      ELSE
          ppm_icommseq(2:ppm_ncommseq(topoid),topoid) = ineighprocs(1:ii,1)
          DO i=1,ppm_nproc-1
              CALL MPI_Send(ineighprocs(1:ii,i+1),ii,MPI_INTEGER,i,i,  &
     &                      ppm_comm,info)
          END DO
      END IF

      !-------------------------------------------------------------------------
      !  Every processor must also communicate to itself in order for
      !  map_part_send to work properly
      !-------------------------------------------------------------------------
      ppm_icommseq(1,topoid) = ppm_rank

      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      ppm_isoptimized(topoid) = .TRUE.

      !-------------------------------------------------------------------------
      !  Deallocate temporary storage
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      IF (ppm_rank .EQ. 0) THEN
          CALL ppm_alloc(optres,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_util_commopt',     &
     &            'optimal result OPTRES',__LINE__,info)
          ENDIF
          CALL ppm_alloc(ilinks,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_util_commopt',     &
     &            'links ILINKS',__LINE__,info)
          ENDIF
      END IF
      CALL ppm_alloc(ineighprocs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_commopt',     &
     &        'neighbor list INEIGHPROCS',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nneighprocs,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_util_commopt',     &
     &        'number of neighbors NNEIGHPROCS',__LINE__,info)
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_commopt',t0,info)
      RETURN
      END SUBROUTINE ppm_util_commopt
