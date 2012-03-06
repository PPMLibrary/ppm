      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_util_commopt
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
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

      SUBROUTINE ppm_util_commopt(topoid,info)
      !!! Determine an approximately optimal communication
      !!! sequence for each processor to SendRecv data with
      !!! its neighbors only. Such that: no conflicts occur
      !!! (i.e. A wants to send to B, but B is currently busy
      !!! receiving from C) and that a minimum number of
      !!! communication rounds are needed. This is done by
      !!! using the Vizing approximation to the minimal edge
      !!! coloring problem of the processor topology graph.
      !!!
      !!! [NOTE]
      !!! The current implementation relies on external routine written in C++
      !!! and its Fortran wrapper.
      !!!
      !!! Reference: V.G. Vizing, On an estimate of the chromatic class
      !!! of a p-graph. Discret. Analiz. 3, 1964, pp.25-30.
      !!! (In Russian).
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_check_id
      USE ppm_module_typedef
      USE ppm_module_color_edge
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
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
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology to be optimized
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success


      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                 :: t0
      INTEGER, DIMENSION(2)                 :: ldu
      INTEGER, DIMENSION(1)                 :: ldu_1D
      INTEGER                               :: iopt
      ! number of neighborhood relations in total
      INTEGER                               :: nlinks
      ! DISTINCT neighbor pairs, i.e. 1<->2 and 2<->1 is the same and only
      ! listed once. even indices are first points, odd ones second ones of
      ! the same edges.
      INTEGER, DIMENSION(:  ) , POINTER     :: ilinks      => NULL()
      ! optimal edge coloring determined. sequence of triples (p1,p2,c),...
      ! with p1 and p2 being the 2 vertices of each edge and c its color.
      INTEGER, DIMENSION(:  ) , POINTER     :: optres      => NULL()
      ! number of neighbors of all every CPU. index: MPI rank
      INTEGER, DIMENSION(:  ) , POINTER     :: nneighprocs => NULL()
      ! all neighbors of all processors in 1D for MPI_Alltoall
      INTEGER, DIMENSION(:  ) , POINTER     :: ineighprocs_1D => NULL()
      ! all neighbors of all processors. 1st index: neighbor nr., 2nd:
      ! processor rank
      INTEGER, DIMENSION(:,:) , POINTER     :: ineighprocs => NULL()
      INTEGER, DIMENSION(:)   , POINTER     :: weights => NULL()
      ! Sources in the graph
      INTEGER, DIMENSION(:)   , POINTER     :: sources=>NULL()
      ! Source weights
      INTEGER, DIMENSION(:)   , POINTER     :: sw=>NULL()
      ! Destinations (targets) in the graph
      INTEGER, DIMENSION(:)   , POINTER     :: dests=>NULL()
      ! Destination weights (edge weights)
      INTEGER, DIMENSION(:)   , POINTER     :: dw=>NULL()
      INTEGER                               :: i,j,maxneigh,isize,ii,isin
      ! processor ranks
      INTEGER                               :: p1,p2
      ! min and max of assigned colors
      INTEGER                               :: mincolor,maxcolor
      ! offset for the MPI_Alltoall send/recv buffer
      INTEGER                               :: offset
      ! new communicator that will be returned after MPI_Graph comm is
      ! created
      INTEGER                               :: newcomm
      ! Total number of neighbors update by rank 0
      INTEGER                               :: total_nneigh
      ! Torsten's library needs this 'n'
      INTEGER                               :: n,reorderedcomm
      LOGICAL                               :: valid
      LOGICAL                               :: file_exists
      CHARACTER(LEN=ppm_char)               :: filename
      INTEGER                               :: inneigh,outneigh,weighted
#ifdef __MPI
      ! MPI comm status
      INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
#endif
      TYPE(ppm_t_topo),      POINTER        :: topo => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
!      EXTERNAL mpix_dist_graph_create
!      EXTERNAL mpix_dist_graph_neighbors_count
!      EXTERNAL mpix_dist_graph_neighbors
!      EXTERNAL tpm_topomap
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_commopt',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Check if there are more than 1 processor. If not, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc .LT. 2) THEN
          !---------------------------------------------------------------------
          !  Allocate memory for communication protocols
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          ldu(1) = 1
          CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &            'communication sequence PPM_ICOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Set the trivial protocol: only myself
          !---------------------------------------------------------------------
          topo%ncommseq = 1
          topo%icommseq(1) = ppm_rank
          GOTO 9999
      END IF

#ifdef __MPI

      !-------------------------------------------------------------------------
      !  Allocate memory for number of neighbors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu_1D(1) = ppm_nproc
      CALL ppm_alloc(nneighprocs,ldu_1D,iopt,info)
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
          CALL MPI_Send(topo%nneighproc,1,MPI_INTEGER,0,ppm_rank,   &
     &                  ppm_comm,info)
      ELSE
          nneighprocs(1) = topo%nneighproc
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
      !  TODO: Should we use a MPI_Alltoall here?
      !-------------------------------------------------------------------------
      IF (ppm_rank .GT. 0) THEN
          CALL MPI_Send(topo%ineighproc(1:topo%nneighproc), &
     &                  topo%nneighproc,  &
     &                  MPI_INTEGER,0,ppm_rank,ppm_comm,info)
      ELSE
          ineighprocs(1:nneighprocs(1),1) =     &
     &                       topo%ineighproc(1:topo%nneighproc)
          DO i=1,ppm_nproc-1
              CALL MPI_Recv(ineighprocs(1:nneighprocs(i+1),i+1),   &
     &                      nneighprocs(i+1),MPI_INTEGER,i,i,ppm_comm,  &
     &                      status,info)
          END DO
      END IF

      !-------------------------------------------------------------------------
      !  Rank 0 broadcasts the total number of neighbors of all the
      !  processors to other MPI processes
      !-------------------------------------------------------------------------
      total_nneigh = SUM(nneighprocs)
      CALL MPI_Bcast(total_nneigh,1,MPI_INTEGER,0,ppm_comm,info)
      CALL MPI_Barrier(ppm_comm,info)

      !-----------------------------------------------------------------------
      !  Allocate memory for number of neighbors in 1D so that MPI_Bcast
      !  would work fine
      !-----------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu_1D(1) = total_nneigh   ! at this point every process has the same
                                 ! total_nneigh
      CALL ppm_alloc(ineighprocs_1D,ldu_1D,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &      'number of neighbors NNEIGHPROCS_1D',__LINE__,info)
        GOTO 9999
      ENDIF
      ineighprocs_1D=0
      IF (ppm_rank .EQ. 0) THEN
!        DO i=0,ppm_nproc-1
!            print*,'neigh of rank=',i,'->',ineighprocs(1:nneighprocs(i+1),i+1)
!        ENDDO
        ineighprocs_1D(1:nneighprocs(1))= topo%ineighproc(1:topo%nneighproc)
        DO i=1,ppm_nproc-1
            offset = i*nneighprocs(i+1)+1
            ineighprocs_1D(offset:offset+nneighprocs(i+1)-1)= &
     &                               ineighprocs(1:nneighprocs(i+1),i+1)
        ENDDO
!        print*,'rank=0->',ineighprocs_1D(1:nneighprocs(1))
!        DO i=1,ppm_nproc-1
!            offset = i*nneighprocs(i+1)+1
!            print*,'NEW neigh of rank=',i,'->',ineighprocs_1D(offset:offset+ &
!     &                      nneighprocs(i+1)-1)
!        ENDDO
      ENDIF
      !-----------------------------------------------------------------------
      !  Rank 0 broadcasts all neighborlists to everyone so that every rank
      !  can create the same MPI_Graph
      !-----------------------------------------------------------------------
      ! Wait for rank 0 to finish up his work
      CALL MPI_Barrier(ppm_comm,info)
      CALL MPI_Bcast(ineighprocs_1D,total_nneigh,MPI_INTEGER,0,ppm_comm,info)
      CALL MPI_Bcast(nneighprocs,ppm_nproc,MPI_INTEGER,0,ppm_comm,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &      'links ILINKS',__LINE__,info)
        GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  The topology file should be provided so that LIBTOPOMAP can be used
      !  TODO: Update the path of the file
      !-------------------------------------------------------------------------
      filename="ppm_processor_topology.txt"
      INQUIRE(FILE=filename, EXIST=file_exists)
      IF (file_exists) THEN
      !-------------------------------------------------------------------------
      !  Using Torsten Hoefler's LIBTOPOMAP to find a better mapping of MPI
      !  processes onto physical processors
      !-------------------------------------------------------------------------

          iopt   = ppm_param_alloc_fit
          ldu(1) = total_nneigh
          ldu(2) = 0
          CALL ppm_alloc(weights,ldu,iopt,info)
          IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
         &      'number of neighbors NNEIGHPROCS_1D',__LINE__,info)
            GOTO 9999
          ENDIF
      ! TODO:
      ! if the graph has weighted edges this 'weights' array needs to be
      ! modified
          weights=1
          n=1
          CALL mpix_dist_graph_create(ppm_comm,n,ppm_rank,nneighprocs,&
     &     ineighprocs_1D,weights,MPI_INFO_NULL,0,topo%mpi_graph_id,info)
          IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &                   'mpix_dist_graph_create failed',__LINE__,info)
            GOTO 9999
          ENDIF
          ppm_comm = topo%mpi_graph_id

          CALL tpm_topomap(ppm_comm,filename,50,ppm_rank,info)

          !  Find neighbors in the graph
          weighted = 1
          CALL mpix_dist_graph_neighbors_count(ppm_comm,inneigh,outneigh,  &
     &                                     weighted,info)
          IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &      'mpix_dist_graph_neighbors_count failed',__LINE__,info)
            GOTO 9999
          ENDIF
          !---------------------------------------------------------------------
          !  Allocate memory for graph sources and destinations as well as
          !  weights. Use always ppm_param_alloc_fit
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          ldu(1) = inneigh
          CALL ppm_alloc(sources,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'sources',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(sw,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'sources',__LINE__,info)
              GOTO 9999
          ENDIF

          ldu(1) = outneigh
          CALL ppm_alloc(dests,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'sources',__LINE__,info)
              GOTO 9999
          ENDIF
          CALL ppm_alloc(dw,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'sources',__LINE__,info)
              GOTO 9999
          ENDIF

          CALL mpix_dist_graph_neighbors(ppm_comm,inneigh,sources,sw,outneigh, &
     &                               dests,dw,info)
          IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',       &
     &      'mpix_dist_graph_neighbors failed',__LINE__,info)
            GOTO 9999
          ENDIF

          CALL MPI_Comm_split(ppm_comm, 0, ppm_rank, reorderedcomm,info);
          topo%mpi_graph_id = reorderedcomm
          ppm_comm = reorderedcomm
      ELSE ! There is no provided processor topology file
          !---------------------------------------------------------------------
          !  Else do not use Libtopomap to create a weighted graph
          !---------------------------------------------------------------------
          CALL MPI_Graph_create(ppm_comm,ppm_nproc,nneighprocs,&
     &                      ineighprocs_1D,0,topo%mpi_graph_id,info)

      ENDIF
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
          ! Call ppm_color_edge to color edges of the given graph such that
          ! every edge of a node has distinct color, using minimum number of
          ! of colors possible (which is degree + 1)
          !---------------------------------------------------------------------
          CALL ppm_color_edge(ppm_nproc,ilinks,optres,info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_util_commopt',     &
     &            'edge coloring failed',__LINE__,info)
              GOTO 9999
          ENDIF
          !CALL vizing_coloring(ppm_nproc,nlinks,ilinks,optres)
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
      !  Distribute protocol to all processors
      !-------------------------------------------------------------------------
      ! First, broadcast the number of rounds to everybody
      IF (ppm_rank .EQ. 0) topo%ncommseq = ii
      CALL MPI_Bcast(topo%ncommseq,1,MPI_INTEGER,0,ppm_comm,info)

      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ! the +1 is due to the fact that the processor itself also needs
      ! to be in the list (as element 1) in order for map_part_send to
      ! work properly
      topo%ncommseq = topo%ncommseq + 1
      ldu(1) = topo%ncommseq
      CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_commopt',     &
     &        'final communication sequence PPM_ICOMMSEQ',__LINE__,info)
          GOTO 9999
      ENDIF

      ! Then distribute the individual optimized neighbor lists
      IF (ppm_rank .GT. 0) THEN
          CALL MPI_Recv(topo%icommseq(2:topo%ncommseq), &
     &                  topo%ncommseq-1,  &
     &                  MPI_INTEGER,0,ppm_rank,ppm_comm,status,info)
      ELSE
          topo%icommseq(2:topo%ncommseq) = ineighprocs(1:ii,1)
          DO i=1,ppm_nproc-1
              CALL MPI_Send(ineighprocs(1:ii,i+1),ii,MPI_INTEGER,i,i,  &
     &                      ppm_comm,info)
          END DO
      END IF

      !-------------------------------------------------------------------------
      !  Every processor must also communicate to itself in order for
      !  map_part_send to work properly
      !-------------------------------------------------------------------------
      topo%icommseq(1) = ppm_rank

      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      topo%isoptimized = .TRUE.

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
      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_commopt',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_util_commopt
