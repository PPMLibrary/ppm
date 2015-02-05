      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_topo_metis_s2p
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_metis_s2p_s(min_phys,max_phys,min_sub,max_sub, &
      &          nneigh,ineigh,cost,nsubs,assig,sub2proc,isublist,nsublist,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_metis_s2p_d(min_phys,max_phys,min_sub,max_sub, &
      &          nneigh,ineigh,cost,nsubs,assig,sub2proc,isublist,nsublist,info)
#endif
      !!! This routine assigns the subdomains to the processors
      !!! using the METIS grid partitioning library version 5.
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:  ), POINTER       :: min_phys
      !!! Minimum of physical extend of the computational domain
      REAL(MK), DIMENSION(:  ), POINTER       :: max_phys
      !!! Maximum of physical extend of the computational domain
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub
      !!! Min. extent of all subdomains
      REAL(MK), DIMENSION(:,:), POINTER       :: max_sub
      !!! Max. extent of all subdomains
      INTEGER,  DIMENSION(:  ), POINTER       :: nneigh
      !!! The number of neighbours of a subdomain
      INTEGER,  DIMENSION(:,:), POINTER       :: ineigh
      !!! Pointers to these neighbours
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      !!! Costs of all subdomains (i.e. sum of particle/mesh node costs in
      !!! that subdomain. Used for weighting the assignment.
      INTEGER,                  INTENT(IN   ) :: nsubs
      !!! Total number of subdomains
      INTEGER,                  INTENT(IN   ) :: assig
      !!! METIS assignment type. One of:
      !!!
      !!! * ppm_param_assign_metis_cut
      !!!  Partitioning a graph of subdomains into k parts minimizing the edgecut
      !!! * ppm_param_assign_metis_comm
      !!!  Partitioning a graph of subdomains into k parts minimizing the total communication volume
      INTEGER,  DIMENSION(:  ), POINTER       :: sub2proc
      !!! Full list of processor affiliation of subdomains
      INTEGER,  DIMENSION(:  ), POINTER       :: isublist
      !!! List of subdomains assigned to the local processors

      INTEGER,                  INTENT(  OUT) :: nsublist
      !!! The number of subdomains assigned
      INTEGER,                  INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK),              DIMENSION(ppm_dim)    :: len_phys
      ! physical extent of comput. domain
      REAL(MK)                                     :: t0,lmyeps
      REAL(MK)                                     :: minsp,maxsp,meansp
      REAL(MK)                                     :: sx,sy,sz
      REAL(MK)                                     :: ex,ey,ez
      REAL(MK)                                     :: lx,ly,lz,ll,kk
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: tpwgts => NULL()
      !This is an array of size nparts*ncon that specifies the desired weight
      !for each partition and constraint.
      REAL(ppm_kind_single), DIMENSION(:), POINTER :: ubvec => NULL()
      !Array for load imbalance

      INTEGER                        :: nvtxs
      !The number of vertices in the graph
      INTEGER                        :: ncon
      !The number of balancing constraints. It should be at least 1.
      INTEGER, DIMENSION(:), POINTER :: xadj   => NULL()
      !The adjacency structure of the graph.
      INTEGER, DIMENSION(:), POINTER :: adjncy => NULL()
      !The adjacency structure of the graph.
      INTEGER, DIMENSION(:), POINTER :: vwgt   => NULL()
      !The weights of the vertices
      INTEGER, DIMENSION(:), POINTER :: vsize  => NULL()
      !The size of the vertices for computing the total communication volume
      INTEGER, DIMENSION(:), POINTER :: adjwgt => NULL()
      !The weights of the edges
      INTEGER                        :: nparts
      !The number of parts to partition the graph.
      !(which is the number of processors)
      INTEGER                        :: objval
      !This variable stores the edge-cut or the total communication
      !volume of the partitioning solution
      INTEGER, DIMENSION(:), POINTER :: part   => NULL()
      !This is a vector which tells each vertex belongs to which partition or processor
      INTEGER, DIMENSION(40)         :: options
      !The maximum length of the options array is 40
      ! (METIS_OPTION_PTYPE   ==1)=0,1     (Multilevel recursive bisectioning,
      !                                     Multilevel k-way partitioning.)
      ! (METIS_OPTION_OBJTYPE ==2)=0,1     (Edge-cut minimization,
      !                                     Total communication volume minimization.)
      ! (METIS_OPTION_CTYPE   ==3)=0,1     (Random matching,Sorted heavy-edge matching.)
      ! (METIS_OPTION_IPTYPE  ==4)=0,1,2,3 (Grows a bisection using a greedy strategy,
      !                                     Computes a bisection at random followed by a refinement,
      !                                     Derives a separator from an edge cut,
      !                                     Grow a bisection using a greedy node-based strategy)
      ! (METIS_OPTION_RTYPE   ==5)=0,1,2,3 (FM-based cut refinement,
      !                                     Greedy-based cut and volume refinement,
      !                                     Two-sided node FM refinement,
      !                                     One-sided node FM refinement.)
      ! (METIS_OPTION_DBGLVL  ==6)=0,.     (The default value is 0
      !                                    (no debugging/progress information))
      ! (METIS_OPTION_NITER,
      ! (METIS_OPTION_NCUTS,
      ! (METIS_OPTION_SEED,
      ! (METIS_OPTION_NO2HOP,
      ! (METIS_OPTION_MINCONN ==11)=0,1 (Does not explicitly minimize the
      !                                  maximum connectivity,
      !                                  Explicitly minimize the maximum connectivity.)
      ! (METIS_OPTION_CONTIG  ==12)=0,1 (Does not force contiguous partitions,
      !                                  Forces contiguous partitions.)
      ! (METIS_OPTION_COMPRESS,
      ! (METIS_OPTION_CCORDER,
      ! (METIS_OPTION_PFACTOR,
      ! (METIS_OPTION_NSEPS,
      ! (METIS_OPTION_UFACTOR,
      ! (METIS_OPTION_NUMBERING ==18)=0,1 (C-style,Fortran-style.)
      INTEGER, DIMENSION(1)            :: ldc
      INTEGER                          :: iopt,isub,jsub
      INTEGER                          :: i,j,nedgs

      CHARACTER(ppm_char) :: caller='ppm_topo_metis_s2p'

      LOGICAL :: lasymm
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Error if METIS library support is not available
      !-------------------------------------------------------------------------
#ifndef __METIS
      nsublist = 0

      fail("PPM was compiled without Metis support", &
      & ppm_err_nometis,ppm_error=ppm_error_fatal)
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate the sub2proc array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs
      CALL ppm_alloc(sub2proc,ldc,iopt,info)
      or_fail_alloc("sub2procs allocation failed",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Check if we only have one processor, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc.EQ.1) THEN
         !----------------------------------------------------------------------
         !  Allocate the isublist array
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         or_fail_alloc("isublist allocation failed",ppm_error=ppm_error_fatal)

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

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      nvtxs=nsubs
      !The number of vertices in the graph.

      ncon=1
      !The number of balancing constraints.

      iopt=ppm_param_alloc_fit
      ldc(1)=nvtxs+1
      CALL ppm_alloc(xadj,ldc,iopt,info)
      or_fail_alloc('METIS adjacency structure of a graph, array XADJ',ppm_error=ppm_error_fatal)

      nedgs=SUM(nneigh(1:nsubs))
      ldc(1)=nedgs
      !The number of edges in the graph
      CALL ppm_alloc(adjncy,ldc,iopt,info)
      or_fail_alloc('METIS adjacency structure of a graph, array ADJNCY',ppm_error=ppm_error_fatal)

      ldc(1)=nvtxs*ncon
      CALL ppm_alloc(vwgt,ldc,iopt,info)
      or_fail_alloc('METIS weights of the vertices, array VWGT',ppm_error=ppm_error_fatal)

      !NULLIFY(vsize)
      ldc(1)=nvtxs
      CALL ppm_alloc(vsize,ldc,iopt,info)
      or_fail_alloc('METIS size of the vertices for the total communication volume, array VSIZE',ppm_error=ppm_error_fatal)

      FORALL (i=1:nvtxs) vsize(i)=0

      ldc(1)=nedgs
      !The number of edges in the graph
      CALL ppm_alloc(adjwgt,ldc,iopt,info)
      or_fail_alloc('METIS weights of the edges, array ADJWGT',ppm_error=ppm_error_fatal)

      nparts=ppm_nproc
      !The number of parts to partition the graph which is the number of processors

      !-------------------------------------------------------------------------
      !  Check if we need asymmetric partitions or not
      !-------------------------------------------------------------------------
      minsp =MINVAL(ppm_proc_speed(0:nparts-1))
      maxsp =MAXVAL(ppm_proc_speed(0:nparts-1))
      meansp=ABS(maxsp-minsp)*REAL(nparts,MK)
      lasymm=meansp.GT.0.05_MK
      ! if there is more than 5 percent difference, do it

      IF (ppm_debug.GT.0) THEN
         stdout_f('(A,F7.4)',"Slowest processor: ",minsp)
         stdout_f('(A,F7.4)',"Fastest processor: ",maxsp)
         stdout_f('(A,F7.4)',"Processor speed imbalance: ",meansp)
         IF (lasymm) THEN
            stdout("doing Asymmetric assignment")
         ELSE
            stdout("doing Symmetric assignment")
         ENDIF
      ENDIF

      IF (lasymm) THEN
         ldc(1) = nparts
         CALL ppm_alloc(tpwgts,ldc,iopt,info)
         or_fail_alloc('METIS desired weight for each (partition) processor, array TPWGTS',ppm_error=ppm_error_fatal)

         !-------------------------------------------------------------------------
         !  Store partition weights for METIS if needed
         !-------------------------------------------------------------------------
         DO i=1,nparts
            tpwgts(i) = REAL(ppm_proc_speed(i-1),ppm_kind_single)
         ENDDO
      ELSE
         NULLIFY(tpwgts)
      ENDIF

      NULLIFY(ubvec)
      !Array for load imbalance

      CALL METIS_SetDefaultOptions(options)
      !Initializes the options array into its default values.

      ldc(1) = nvtxs
      CALL ppm_alloc(part,ldc,iopt,info)
      or_fail_alloc('METIS partition vector of the graph which tells each vertex belongs to which partition or processor, array PART',ppm_error=ppm_error_fatal)

      !----------------------------------------------------------------
      !  Fill the adjacency structure of the graph
      !----------------------------------------------------------------
      xadj(1)=1
      DO isub=2,nsubs+1
         xadj(isub)=xadj(isub-1)+nneigh(isub-1)
      ENDDO

      j=0
      DO isub=1,nsubs
         DO i=1,nneigh(isub)
            j=j+1
            adjncy(j)=ineigh(i,isub)
         ENDDO
      ENDDO

      FORALL (isub=1:nsubs) vwgt(isub)=INT(cost(isub)) ! computation weights

      len_phys(1:ppm_dim) = max_phys(1:ppm_dim) - min_phys(1:ppm_dim)

      sx=MINVAL(max_sub(1,1:nsubs)-min_sub(1,1:nsubs))/len_phys(1)
      sy=MINVAL(max_sub(2,1:nsubs)-min_sub(2,1:nsubs))/len_phys(2)
      !Finding the smallest edge, we can find the biggest factor of 10
      !in order to convert the REAL communication Coefficient into INTEGER weight


      !----------------------------------------------------------------
      !  Fill the weights of the edge array (adjwgt).
      !  Here we consider the length (area) of the neighbor subdomain
      !  as weight for the edge which connects two neighbour subdomains
      !  by an edge.
      !  The edge (area) is a coefficient which is related to the approximate
      !  amount of communications between two adjacent neighbor subdomains.
      !
      !  Fill the total communication volume array (vsize).
      !  Total communication on each vertex is considered as a summation
      !  of all weights of the edge array of neighbour subdomains.
      !----------------------------------------------------------------

      SELECT CASE (ppm_dim)
      CASE (2)
         ll=MIN(sx,sy)

         kk=1._MK
         DO WHILE (ll.LE.1.0_MK)
            kk=kk*10._MK
         ENDDO
         !Coefficient to convert REAL number into INTEGER value

         j=0
         DO isub=1,nsubs
            DO i=1,nneigh(isub)
               jsub=ineigh(i,isub)

               sx=MAX(min_sub(1,isub),min_sub(1,jsub))
               ex=MIN(max_sub(1,isub),max_sub(1,jsub))
               lx=ABS(ex-sx)

               sy=MAX(min_sub(2,isub),min_sub(2,jsub))
               ey=MIN(max_sub(2,isub),max_sub(2,jsub))
               ly=ABS(ey-sy)

               j=j+1
               IF (lx.LE.lmyeps) THEN
                  adjwgt(j)=INT(ly/len_phys(2)*kk)
               ELSE IF (ly.LE.lmyeps) THEN
                  adjwgt(j)=INT(lx/len_phys(1)*kk)
               ENDIF

               vsize(isub)=vsize(isub)+adjwgt(j)

            ENDDO
         ENDDO
      CASE (3)
         sz=MINVAL(max_sub(3,1:nsubs)-min_sub(3,1:nsubs))/len_phys(3)

         ll=MIN(sx*sy,sx*sz,sy*sz)

         kk=1._MK
         DO WHILE (ll.LE.1.0_MK)
            kk=kk*10._MK
         ENDDO

         j=0
         DO isub=1,nsubs
            DO i=1,nneigh(isub)
               jsub=ineigh(i,isub)

               sx=MAX(min_sub(1,isub),min_sub(1,jsub))
               ex=MIN(max_sub(1,isub),max_sub(1,jsub))
               lx=ABS(ex-sx)

               sy=MAX(min_sub(2,isub),min_sub(2,jsub))
               ey=MIN(max_sub(2,isub),max_sub(2,jsub))
               ly=ABS(ey-sy)

               sz=MAX(min_sub(3,isub),min_sub(3,jsub))
               ez=MIN(max_sub(3,isub),max_sub(3,jsub))
               lz=ABS(ez-sz)

               j=j+1
               IF (lx.LE.lmyeps) THEN
                  adjwgt(j)=INT(ly*lz/(len_phys(2)*len_phys(3))*kk)
               ELSE IF (ly.LE.lmyeps) THEN
                  adjwgt(j)=INT(lx*lz/(len_phys(1)*len_phys(3))*kk)
               ELSE IF (lz.LE.lmyeps) THEN
                  adjwgt(j)=INT(lx*ly/(len_phys(1)*len_phys(2))*kk)
               ENDIF

               vsize(isub)=vsize(isub)+adjwgt(j)

            ENDDO
         ENDDO
      END SELECT

      options( 1)=MERGE(0,1,ppm_nproc.LE.8)
      options( 3)=1 !Sorted heavy-edge matching
      options( 4)=0 !Grows a bisection using a greedy strategy.
      options( 5)=3 !One-sided node FM refinement
      options( 6)=0 !(no debugging/progress information)
      options(18)=1 !Fortran-style numbering is assumed that starts from 1.

      !-------------------------------------------------------------------------
      !  Partition graph into nparts parts using METIS and the minimal
      !  edgecut objective, egdecut is the number of mesh edges that
      !  has been cut by the partition.
      !-------------------------------------------------------------------------
      SELECT CASE (assig)
      CASE (ppm_param_assign_metis_cut)
         options(2)=0 !Edge-cut minimization
      CASE (ppm_param_assign_metis_comm)
         options(2)=1 !Total communication volume minimization.
      CASE DEFAULT
         fail('Unknown METIS assignment scheme. Bailing out.')
      END SELECT

      SELECT CASE (ppm_nproc)
      CASE (2:8)
         CALL METIS_PartGraphRecursive(nvtxs,ncon,xadj,adjncy,vwgt, &
         &    vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part)

      CASE DEFAULT
         CALL METIS_PartGraphKway(nvtxs,ncon,xadj,adjncy,vwgt, &
         &    vsize,adjwgt,nparts,tpwgts,ubvec,options,objval,part)

      END SELECT

      !----------------------------------------------------------------------
      !  Output diagnostics
      !----------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         SELECT CASE (assig)
         CASE (ppm_param_assign_metis_cut)
            stdout_f('(A,I6)',"METIS returned edgecut = ",objval)
         CASE (ppm_param_assign_metis_comm)
            stdout_f('(A,I6)',"METIS returned communication volume = ",objval)
         END SELECT
      ENDIF

      !-------------------------------------------------------------------------
      !  Find subs2proc assignment based on graph decomposition
      !-------------------------------------------------------------------------
      ! a graph vertex is a sub
      FORALL (isub=1:nsubs) sub2proc(isub) = part(isub)-1    ! MPI ranks start at 0

      !-------------------------------------------------------------------------
      !  Deallocate METIS work memory
      !-------------------------------------------------------------------------
      iopt=ppm_param_dealloc
      CALL ppm_alloc(xadj,ldc,iopt,info)
      or_fail_dealloc('METIS adjacency structure of a graph, array XADJ',exit_point=no)

      CALL ppm_alloc(adjncy,ldc,iopt,info)
      or_fail_dealloc('METIS adjacency structure of a graph, array ADJNCY',exit_point=no)

      CALL ppm_alloc(vwgt,ldc,iopt,info)
      or_fail_dealloc('METIS weights of the vertices, vector VWGT',exit_point=no)

      CALL ppm_alloc(vsize,ldc,iopt,info)
      or_fail_dealloc('METIS size of the vertices for the total communication volume, array VSIZE',exit_point=no)

      CALL ppm_alloc(adjwgt,ldc,iopt,info)
      or_fail_dealloc('METIS weights of the edges, array ADJWGT',exit_point=no)

      IF (lasymm) THEN
         CALL ppm_alloc(tpwgts,ldc,iopt,info)
         or_fail_dealloc('METIS graph adjacency array TPWGTS',exit_point=no)
      ELSE
         NULLIFY(tpwgts)
      ENDIF

      CALL ppm_alloc(part,ldc,iopt,info)
      or_fail_dealloc('METIS partition vector of the graph which tells each vertex belongs to which partition or processor, array PART',exit_point=no)

      !-------------------------------------------------------------------------
      !  Count number of subs assigned to local processor
      !-------------------------------------------------------------------------
      nsublist = COUNT(sub2proc(1:nsubs).EQ.ppm_rank)

      !-------------------------------------------------------------------------
      !  Allocate the isublist array
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsublist
      CALL ppm_alloc(isublist,ldc,iopt,info)
      or_fail_alloc('isublist allocation failed',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Fill list of local subs
      !-------------------------------------------------------------------------
      isublist(1:nsublist)=PACK(sub2proc(1:nsubs),sub2proc(1:nsubs).EQ.ppm_rank)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (nsubs .LE. 0) THEN
            fail('nsubs must be > 0',exit_point=8888)
         ENDIF
         IF (assig .NE. ppm_param_assign_metis_cut .AND.  &
         &   assig .NE. ppm_param_assign_metis_comm) THEN
            fail('Invalid assignment type for METIS passed',exit_point=8888)
         ENDIF
         DO isub=1,nsubs
            IF (nneigh(isub) .LT. 0) THEN
               fail('nneigh must be >= 0 for all subs',exit_point=8888)
            ENDIF
            IF (cost(isub) .LT. 0) THEN
               fail('costs must be >= 0 for all subs',exit_point=8888)
            ENDIF
         ENDDO
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_metis_s2p_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_metis_s2p_d
#endif
