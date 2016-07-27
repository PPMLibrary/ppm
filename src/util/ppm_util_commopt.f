      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_util_commopt
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
      USE ppm_module_mpi
      USE ppm_module_topo_typedef
      USE ppm_module_color
      IMPLICIT NONE

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: topoid
      !!! The topology to be optimized
      INTEGER, INTENT(INOUT) :: info
      !!! return status, 0 on success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(1)                :: ldl
      INTEGER, DIMENSION(1)                :: ldu
      INTEGER                              :: iopt
      ! number of neighborhood relations in total
      INTEGER                              :: nlinks
      ! DISTINCT neighbor pairs, i.e. 1<->2 and 2<->1 is the same and only
      ! listed once. even indices are first points, odd ones second ones of
      ! the same edges.
      INTEGER, DIMENSION(:),   ALLOCATABLE :: ilinks,ilinkst
      ! optimal edge coloring determined. sequence of triples (p1,p2,c),...
      ! with p1 and p2 being the 2 vertices of each edge and c its color.
      INTEGER, DIMENSION(:),   ALLOCATABLE :: optres
      ! number of neighbors of all every CPU. index: MPI rank
      INTEGER, DIMENSION(:),   ALLOCATABLE :: nneighprocs
      ! all neighbors of all processors. 1st index: neighbor nr., 2nd:
      ! processor rank
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ineighprocs


      INTEGER, DIMENSION(:),   ALLOCATABLE :: coloring
      ! processor colors
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: ineighcolors

      INTEGER, DIMENSION(:),   ALLOCATABLE :: counts
      INTEGER, DIMENSION(:),   ALLOCATABLE :: displ
      INTEGER                              :: i
      INTEGER                              :: j
      INTEGER                              :: k
      INTEGER                              :: ii
      INTEGER                              :: maxneigh,isize,isin
      ! processor ranks
      INTEGER                              :: p1,p2
      ! min and max of assigned colors
      INTEGER                              :: mincolor,maxcolor

      CHARACTER(LEN=*), PARAMETER :: caller="ppm_util_commopt"

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Check if there are more than 1 processor. If not, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc.LT.2) THEN
         !---------------------------------------------------------------------
         !  Allocate memory for communication protocols
         !---------------------------------------------------------------------
         iopt=ppm_param_alloc_fit
         ldu  = 1
         CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
         or_fail_alloc("communication sequence PPM_ICOMMSEQ",ppm_error=ppm_error_fatal)

         !---------------------------------------------------------------------
         !  Set the trivial protocol: only myself
         !---------------------------------------------------------------------
         topo%ncommseq = 1
         topo%icommseq = ppm_rank

         !---------------------------------------------------------------------
         !  One processor has only one color
         !---------------------------------------------------------------------
         topo%ineighcolor = 1

         GOTO 9999
      ENDIF

#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Allocate memory for number of neighbors
      !-------------------------------------------------------------------------
      ALLOCATE(nneighprocs(ppm_nproc),SOURCE=0,STAT=info)
      or_fail_alloc("Failed to allocate number of neighbors NNEIGHPROCS",ppm_error=ppm_error_fatal)

      ALLOCATE(displ(ppm_nproc),counts(ppm_nproc),STAT=info)
      or_fail_alloc("displ, counts",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Rank 0 receives the numbers of neighbors from all other processors
      !-------------------------------------------------------------------------
      CALL MPI_Gather(topo%nneighproc,1,MPI_INTEGER,nneighprocs,1,MPI_INTEGER,0,ppm_comm,info)
      or_fail_MPI("MPI_Gather")

      maxneigh=MAXVAL(nneighprocs)

      counts=maxneigh

      DO i=1,ppm_nproc
         displ(i)=(i-1)*maxneigh
      ENDDO

      !-------------------------------------------------------------------------
      !  Allocate memory for neighbor lists
      !-------------------------------------------------------------------------
      ALLOCATE(ineighprocs(maxneigh,ppm_nproc),STAT=info)
      or_fail_alloc("neighbor list INEIGHPROCS",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for neighbor colors
      !-------------------------------------------------------------------------
      ALLOCATE(ineighcolors(0:maxneigh+1,ppm_nproc),STAT=info)
      or_fail_alloc("neighbor colors INEIGHCOLORS",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Rank 0 receives all neighbor lists
      !-------------------------------------------------------------------------
      CALL MPI_Gatherv(topo%ineighproc,topo%nneighproc,MPI_INTEGER, &
      &    ineighprocs,counts,displ,MPI_INTEGER,0,ppm_comm,info)
      or_fail_MPI("MPI_Gather")

      !-------------------------------------------------------------------------
      !  Rank 0: Build graph and call optimizer
      !-------------------------------------------------------------------------
      IF (ppm_rank.EQ.0) THEN
         !---------------------------------------------------------------------
         !  Build graph edges as UNIQUE neighbor pairs
         !---------------------------------------------------------------------
         isize=2*ppm_nproc   ! initial guess
         ALLOCATE(ilinks(isize),STAT=info)
         or_fail_alloc("links ILINKS",ppm_error=ppm_error_fatal)

         nlinks = 0
         DO i=1,ppm_nproc
            DO j=1,nneighprocs(i)
               ! Check is the reverse link is already in the list
               isin = 0
               DO ii=1,nlinks
                  IF ((ilinks(2*ii-1).EQ.(ineighprocs(j,i)+1)).AND. &
                  &   (ilinks(2*ii).EQ.i)) isin = 1
               ENDDO
               ! add link if reverse is not already in the list
               IF (isin.EQ.0) THEN
                  nlinks = nlinks + 1
                  ! every link will cause two entries in ilinks !!
                  IF (2*nlinks.GT.isize) THEN
                     ALLOCATE(ilinkst(isize+10),STAT=info)
                     or_fail_alloc("links ILINKS",ppm_error=ppm_error_fatal)

                     FORALL (k=1:isize) ilinkst(k)=ilinks(k)

                     CALL MOVE_ALLOC(ilinkst,ilinks)

                     isize = isize + 10
                   ENDIF

                   ilinks(2*nlinks-1) = i
                   ! in MPI numbering starts at 0, but here at 1 :-)
                   ilinks(2*nlinks)   = ineighprocs(j,i)+1
                ENDIF !(isin.EQ.0)
            ENDDO !j=1,nneighprocs(i)
         ENDDO !i=1,ppm_nproc

         !---------------------------------------------------------------------
         !  Allocate memory for optimal coloring result
         !---------------------------------------------------------------------
         ALLOCATE(optres(3*nlinks),STAT=info)
         or_fail_alloc("optimal result OPTRES",ppm_error=ppm_error_fatal)

         !---------------------------------------------------------------------
         !  Shrink ilinks to the correct size. C does not like it
         !  otherwise.
         !---------------------------------------------------------------------
         k=2*nlinks
         ALLOCATE(ilinkst(k),SOURCE=ilinks(1:k),STAT=info)
         or_fail_alloc("final links ILINKS",ppm_error=ppm_error_fatal)

         CALL MOVE_ALLOC(ilinkst,ilinks)

         !---------------------------------------------------------------------
         ! Call ppm_color_edge to color edges of the given graph such that
         ! every edge of a node has distinct color, using minimum number of
         ! of colors possible (which is degree + 1)
         !---------------------------------------------------------------------
         CALL ppm_color_edge(ppm_nproc,ilinks,optres,info)
         or_fail('edge coloring failed')
         !---------------------------------------------------------------------
         !  optres now contains the result as a sequence of nlinks triples
         !  (p1,p2,c) where each triple is one edge (determined by its two
         !  vertices p1 and p2) with a color c. Numbering of colors and
         !  vertices starts at 1.
         !---------------------------------------------------------------------

         !-------------------------------------------------------------------------
         !  Deallocate temporary storage
         !-------------------------------------------------------------------------
         DEALLOCATE(ilinks,STAT=info)
         or_fail_dealloc("links ILINKS")

         !-------------------------------------------------------------------------
         !  Vertex coloring
         !-------------------------------------------------------------------------
         ALLOCATE(coloring(ppm_nproc),STAT=info)
         or_fail_alloc("coloring",ppm_error=ppm_error_fatal)

         ineighprocs=ineighprocs+1

         CALL ppm_color_vertex(ppm_nproc,nneighprocs,ineighprocs,coloring,info)
         or_fail("ppm_color_vertex")

         ii=MAXVAL(coloring)
         DO i=1,ppm_nproc
            ineighcolors(0,i)=coloring(i)
            DO j=1,nneighprocs(i)
               k=ineighprocs(j,i)
               ineighcolors(j,i)=coloring(k)
            ENDDO
            ineighcolors(j,i)=ii
         ENDDO

         !---------------------------------------------------------------------
         !  Allocate memory for neighbor communication lists
         !---------------------------------------------------------------------
         DEALLOCATE(ineighprocs,coloring,STAT=info)
         or_fail_dealloc("ineighprocs & coloring",ppm_error=ppm_error_fatal)

         !---------------------------------------------------------------------
         !  Determine smallest and largest assigned color
         !---------------------------------------------------------------------
         mincolor =  ppm_big_i
         maxcolor = -ppm_big_i
         DO i=3,3*nlinks,3
            IF (optres(i).LT.mincolor) mincolor = optres(i)
            IF (optres(i).GT.maxcolor) maxcolor = optres(i)
         ENDDO

         !---------------------------------------------------------------------
         !  Allocate memory for neighbor communication lists
         !---------------------------------------------------------------------
         ii = maxcolor-mincolor+1
         ALLOCATE(ineighprocs(ii,ppm_nproc),STAT=info)
         or_fail_alloc('neighbor comm. list INEIGHPROCS',ppm_error=ppm_error_fatal)

         ! initialize array. -1 means: no communication in this round
         ineighprocs=-1

         !This is the displ data for MPI_Scatterv
         DO i=1,ppm_nproc
            displ(i)=(i-1)*ii
         ENDDO

         !---------------------------------------------------------------------
         !  Assemble communication protocol
         !---------------------------------------------------------------------
         ii = 0
         DO j=mincolor,maxcolor
            ii = ii + 1
            DO i=3,3*nlinks,3
               ! if this link has the correct color
               IF (optres(i).EQ.j) THEN
                  p1 = optres(i-1)
                  p2 = optres(i-2)
                  ! MPI ranks start at 0.
                  ineighprocs(ii,p2) = p1-1
                  ineighprocs(ii,p1) = p2-1
               ENDIF
            ENDDO
         ENDDO

         topo%ncommseq = ii

         !-------------------------------------------------------------------------
         !  Deallocate temporary storage
         !-------------------------------------------------------------------------
         DEALLOCATE(optres,STAT=info)
         or_fail_dealloc("optimal result OPTRES")
      ENDIF ! ppm_rank.EQ.0

      !-------------------------------------------------------------------------
      !  Distribute protocol to all processors
      !-------------------------------------------------------------------------
      ! First, broadcast the number of rounds to everybody
      CALL MPI_Bcast(topo%ncommseq,1,MPI_INTEGER,0,ppm_comm,info)
      or_fail_MPI("MPI_Bcast")

      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_grow_preserve
      ! the +1 is due to the fact that the processor itself also needs
      ! to be in the list (as element 1) in order for map_part_send to
      ! work properly
      ldu(1) = topo%ncommseq+1
      CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
      or_fail_alloc('final communication sequence PPM_ICOMMSEQ',ppm_error=ppm_error_fatal)

      ! Then distribute the individual optimized neighbor lists
      counts=topo%ncommseq
      CALL MPI_Scatterv(ineighprocs,counts,displ,MPI_INTEGER, &
      &    topo%icommseq(2),topo%ncommseq,MPI_INTEGER,0,ppm_comm,info)
      or_fail_MPI("MPI_Scatterv")

      topo%ncommseq=topo%ncommseq+1

      !-------------------------------------------------------------------------
      !  Every processor must also communicate to itself in order for
      !  map_part_send to work properly
      !-------------------------------------------------------------------------
      topo%icommseq(1)=ppm_rank

      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      topo%isoptimized=.TRUE.


      maxneigh=maxneigh+2
      DO i=1,ppm_nproc
         displ(i)=(i-1)*maxneigh
      ENDDO

      nneighprocs=nneighprocs+2

      CALL MPI_Scatterv(ineighcolors,nneighprocs,displ,MPI_INTEGER, &
      &    topo%ineighcolor(0),topo%nneighproc+2,MPI_INTEGER,0,ppm_comm,info)
      or_fail_MPI("MPI_Scatterv")

      !-------------------------------------------------------------------------
      !  Deallocate temporary storage
      !-------------------------------------------------------------------------
      DEALLOCATE(ineighprocs,ineighcolors,nneighprocs,STAT=info)
      or_fail_dealloc("neighbor list INEIGHPROCS, ineighcolors & number of neighbors NNEIGHPROCS")

      DEALLOCATE(displ,counts,STAT=info)
      or_fail_dealloc("displ & counts")
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IMPLICIT NONE
        LOGICAL :: valid
        CALL ppm_check_topoid(topoid,valid,info)
        IF (.NOT.valid) THEN
           fail("topoid out of range",exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_util_commopt
