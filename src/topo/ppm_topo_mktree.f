      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mktree
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

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_topo_mktree_s(topoid,meshid,xp,Npart,assig,min_phys,        &
      &          max_phys,bcdef,cost,minboxsize,pruneboxes,weights,fixed,        &
      &          maxvariance,maxboxcost,istart,ndata,Nm,storemesh,treetype,info, &
      &          pcost,user_minsub,user_maxsub,user_nsubs,user_sub2proc)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mktree_d(topoid,meshid,xp,Npart,assig,min_phys,        &
      &          max_phys,bcdef,cost,minboxsize,pruneboxes,weights,fixed,        &
      &          maxvariance,maxboxcost,istart,ndata,Nm,storemesh,treetype,info, &
      &          pcost,user_minsub,user_maxsub,user_nsubs,user_sub2proc)
#endif
      !!! This routine creates a topology based on ppm_tree.
      !!! The user can directly specify the arguments for ppm_tree, rather
      !!! than choosing from a set of predefined `ppm_param_decomp_*`.
      !!! Please see documentation of `ppm_tree` for description of the
      !!! arguments.
      !!!
      !!! [WARNING]
      !!! This is for advanced users only.
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data_mesh
      !USE ppm_module_mesh_alloc
      USE ppm_module_mesh_define
      USE ppm_module_mesh_on_subs
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
      INTEGER                 , INTENT(INOUT) :: topoid
      !!! ID number identifying the topology.
      !!! If topoid == 0 on input a new topology is created and the new topoid
      !!! is returned here, else the indicated toplogy is replaced.
      !!!
      !!! [CAUTION]
      !!! *SEMANTICS CHANGED:* `topoid = 0` is *not anymore* reserved for the
      !!! ring topology (null decomposition). "Ring topologies" are non
      !!! geometric and need no setup. The user can perform ppm ring shift
      !!! operations without having to first define a topology.
      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      !!! Particle positions. If present, domain decomposition will be
      !!! guided by particle locations
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles. Set to something .LE.0 if there are no
      !!! particles.
      INTEGER                 , INTENT(IN   ) :: assig
      !!! The type of subdomain-to-processor assignment. One of:
      !!!
      !!! *  ppm_param_assign_internal
      !!! *  ppm_param_assign_metis_cut
      !!! *  ppm_param_assign_metis_comm
      !!! *  ppm_param_assign_user_defined
      !!!
      !!! [NOTE]
      !!! The latter uses the external library METIS and is only
      !!! available if ppm was compiled with METIS support.
      REAL(MK), DIMENSION(:  ), POINTER       :: min_phys
      !!! Minimum of physical extend of the computational domain (double)
      !!!
      !!! first index is ppm_dim
      REAL(MK), DIMENSION(:  ), POINTER       :: max_phys
      !!! Maximum of physical extend of the computational domain (double)
      !!!
      !!! first index is ppm_dim
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      !!! Boundary conditions for the topology
      !!!
      !!! NOTE: first index is 1-6 (each of the faces)
      !!! - west  : 1
      !!! - east  : 2
      !!! - south : 3
      !!! - north : 4
      !!! - bottom: 5
      !!! - top   : 6
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      !!! Estimated cost associated with subdomains. Either user-defined on
      !!! input or decomposition result on output. The cost of a subdomain
      !!! is given by its volume.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Array of length Npart which specifies the computational
      !!! cost attributed to each particle. If this is absent, a unity cost is
      !!! assumed for each particle.
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER :: user_minsub
      !!! Mimimum of extension of subs. OPTIONAL parameter used if
      !!! decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER :: user_maxsub
      !!! maximum of extension of subs.
      !!! Used if decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      INTEGER                 , OPTIONAL          :: user_nsubs
      !!! Total number of subs on all processors.
      !!! Used when decomp is user defined
      INTEGER, DIMENSION(:  ),  OPTIONAL, POINTER :: user_sub2proc
      !!! Subdomain to processor assignment.
      !!! Used if assignment is user defined.
      !!!
      !!! index: subID (global)
      REAL(MK)                , INTENT(IN   ) :: minboxsize
      !!! The min size of the boxes
      LOGICAL                 , INTENT(IN   ) :: pruneboxes
      !!! Whether to keep empty boxes or not?
      REAL(MK), DIMENSION(3,2), INTENT(IN   ) :: weights
      !!! Weights for particles, mesh points and geometry. `weights(:,1)`
      !!! defines the box costs (where to refine), `weights(:,2)` defines
      !!! where to cut a box.
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      !!! Flags to prevent cuts in certain directions. (see ppm_tree)
      REAL(MK)                , INTENT(IN   ) :: maxvariance
      !!! Maximum allowed variance of costs between subs
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      !!! Maximum cost of a box. Subdivision will continue (if  possible)
      !!! if cost of any box is larger than this. Set to <= 0
      !!! to disable this criterion.
      INTEGER                 , INTENT(INOUT) :: meshid
      !!! Mesh identifier (user numbering) for default mesh (as defined by
      !!! Nm) on the topology. If <0, ppm will create a number internally
      !!! and return it here on exit.
      INTEGER , DIMENSION(:,:), POINTER       :: istart
      !!! Start indices (i,j[,k]) (first index) of sub mesh isub (second
      !!! index) in global mesh.
      INTEGER , DIMENSION(:,:), POINTER       :: ndata
      !!! Number of grid points in x,y[,z] (first index) of sub mesh isub
      !!! (second index). Includes the points ON the sub boundaries.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Number of mesh points (not cells) in each direction of the global
      !!! computational domain (including points ON its boundaries)
      LOGICAL                 , INTENT(IN   ) :: storemesh
      !!! If `TRUE` store the mesh definition internally. If `FALSE`,
      !!! the mesh serves only as a guide for the decomposition and is not
      !!! stored as a ppm compute mesh.
      INTEGER                 , INTENT(IN   ) :: treetype
      !!! One of:
      !!!
      !!! *ppm_param_tree_bin
      !!! *ppm_param_tree_quad
      !!! *ppm_param_tree_oct (3D only)
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(:,:), POINTER :: min_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: min_sub => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_sub => NULL()

      INTEGER                          :: i,nbox,iopt,isub
      INTEGER, DIMENSION(2  )          :: ldc
      INTEGER, DIMENSION(:,:), POINTER :: ineigh   => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: subs_bc  => NULL()
      INTEGER, DIMENSION(:  ), POINTER :: nneigh   => NULL()
      INTEGER, DIMENSION(:  ), POINTER :: nchld    => NULL()
      INTEGER                          :: nsublist, nsubs
      INTEGER, DIMENSION(  :), POINTER :: isublist => NULL()
      INTEGER, DIMENSION(:  ), POINTER :: sub2proc => NULL()

      CHARACTER(LEN=ppm_char) :: caller="ppm_topo_mktree"
      CHARACTER(LEN=ppm_char) :: mesg

      LOGICAL :: have_particles,have_mesh

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      ! If the user defined nsubs then use those
      nsubs=MERGE(user_nsubs,0,PRESENT(user_nsubs))

      IF (PRESENT(user_minsub)) THEN
         min_sub => user_minsub
      ENDIF
      IF (PRESENT(user_maxsub)) THEN
         max_sub => user_maxsub
      ENDIF
      IF (PRESENT(user_sub2proc)) THEN
         sub2proc => user_sub2proc
      ENDIF

      !-------------------------------------------------------------------------
      !  Check if we have particles and mesh
      !-------------------------------------------------------------------------
      have_particles=Npart.GT.0
      have_mesh     =.FALSE.
      IF (SIZE(Nm,1).GE.ppm_dim) THEN
         have_mesh=ALL(Nm(1:ppm_dim).GT.1)
      ENDIF

      !-------------------------------------------------------------------------
      !  Call the tree routine
      !-------------------------------------------------------------------------
      gsvec(1:ppm_dim) = minboxsize

      ! build tree
      IF (PRESENT(pcost)) THEN
         CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
         &    ppm_nproc,pruneboxes,gsvec,maxvariance,maxboxcost, &
         &    fixed,weights,min_box,max_box,nbox,nchld,info,pcost)
      ELSE
         CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
         &    ppm_nproc,pruneboxes,gsvec,maxvariance,maxboxcost, &
         &    fixed,weights,min_box,max_box,nbox,nchld,info)
      ENDIF
      or_fail('General tree decomposition failed')

      ! convert tree to subs
      CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox, &
      &    min_sub,max_sub,nsubs,info)
      IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,nneigh,ineigh,gsvec,info)
      or_fail('Finding neighbors failed')

      !-------------------------------------------------------------------------
      !  Allocate mesh data to avoid errors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = 0
      CALL ppm_alloc(istart,ldc,iopt,info)
      or_fail_alloc('mesh start indices ISTART',ppm_error=ppm_error_fatal)

      CALL ppm_alloc(ndata,ldc,iopt,info)
      or_fail_alloc('mesh sizes NDATA',ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Define meshes on the subs
      !-------------------------------------------------------------------------
      IF (have_mesh) THEN
         CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,min_sub,max_sub,nsubs, &
         &   istart,ndata,info)
         or_fail('Defining meshes failed')
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (PRESENT(pcost)) THEN
         CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,info,pcost)
      ELSE
         CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,info)
      ENDIF
      or_fail('Computing costs failed')

      !-------------------------------------------------------------------------
      !  Assign the subdomains to processors
      !-------------------------------------------------------------------------
      NULLIFY(isublist)

      SELECT CASE (assig)
      CASE (ppm_param_assign_internal)
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, &
         &    isublist,nsublist,info)
         or_fail('Assigning subs to processors failed')

      CASE (ppm_param_assign_metis_cut,ppm_param_assign_metis_comm)
         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_phys,max_phys,min_sub,max_sub, &
         &    cost,nneigh,ineigh,nsubs,sub2proc,isublist,nsublist,assig,info)
         or_fail('Assigning subs to processors using METIS failed')

      CASE (ppm_param_assign_user_defined)
         !----------------------------------------------------------------------
         !  user defined assignment
         !----------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         or_fail_alloc('list of local subs ISUBLIST',ppm_error=ppm_error_fatal)

         isublist = ppm_param_undefined
         nsublist = 0
         DO isub=1,nsubs
            IF (sub2proc(isub) .EQ. ppm_rank) THEN
               nsublist = nsublist + 1
               isublist(nsublist) = isub
            ENDIF
         ENDDO

      CASE DEFAULT
         !-------------------------------------------------------------------
         !  unknown assignment scheme
         !-------------------------------------------------------------------
         WRITE(mesg,'(A,I5)') 'Unknown assignment scheme: ',assig
         fail(mesg)

      END SELECT

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,subs_bc,info)
      or_fail('finding and defining the BC of the subs failed ')

      !-------------------------------------------------------------------------
      !  Store the topology internally
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(topoid,min_phys,max_phys,min_sub,max_sub,subs_bc, &
      &    sub2proc,nsubs,bcdef,minboxsize,isublist,nsublist,nneigh,ineigh,info,0)
      or_fail('Storing topology failed')

      IF (storemesh.AND.have_mesh) THEN
         !---------------------------------------------------------------------
         !  Create and store new mesh internally
         !---------------------------------------------------------------------
         CALL ppm_mesh_define(topoid,meshid,Nm,info)
         or_fail('Storing mesh definition failed')
      ENDIF ! storemesh

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

      iopt = ppm_param_dealloc
      CALL ppm_alloc(nchld,ldc,iopt,info)
      or_fail_dealloc("nchld")

      CALL ppm_alloc(min_box,ldc,iopt,info)
      or_fail_dealloc("min_box")

      CALL ppm_alloc(max_box,ldc,iopt,info)
      or_fail_dealloc("max_box")

      IF (PRESENT(user_minsub)) THEN
         IF (ASSOCIATED(min_sub,user_minsub)) THEN
            NULLIFY(min_sub)
         ELSE
            CALL ppm_alloc(min_sub,ldc,iopt,info)
            or_fail_dealloc("min_sub")
         ENDIF
      ELSE
         CALL ppm_alloc(min_sub,ldc,iopt,info)
         or_fail_dealloc("min_sub")
      ENDIF
      IF (PRESENT(user_maxsub)) THEN
         IF (ASSOCIATED(max_sub,user_maxsub)) THEN
            NULLIFY(max_sub)
         ELSE
            CALL ppm_alloc(max_sub,ldc,iopt,info)
            or_fail_dealloc("max_sub")
         ENDIF
      ELSE
         CALL ppm_alloc(max_sub,ldc,iopt,info)
         or_fail_dealloc("max_sub")
      ENDIF
      IF (PRESENT(user_sub2proc)) THEN
         IF (ASSOCIATED(sub2proc,user_sub2proc)) THEN
            NULLIFY(sub2proc)
         ELSE
            CALL ppm_alloc(sub2proc,ldc,iopt,info)
            or_fail_dealloc("sub2proc")
         ENDIF
      ELSE
         CALL ppm_alloc(sub2proc,ldc,iopt,info)
         or_fail_dealloc("sub2proc")
      ENDIF

      CALL ppm_alloc(ineigh,ldc,iopt,info)
      or_fail_dealloc("ineigh")

      CALL ppm_alloc(nneigh,ldc,iopt,info)
      or_fail_dealloc("nneigh")

      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      or_fail_dealloc("subs_bc")

      CALL ppm_alloc(isublist,ldc,iopt,info)
      or_fail_dealloc("isublist")

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!",ppm_err_ppm_noinit,exit_point=8888)
         ENDIF
         IF (minboxsize .LT. 0.0_MK) THEN
            fail("minboxsize must be >= 0.0",exit_point=8888)
         ENDIF
         DO i=1,ppm_dim
            IF (max_phys(i).LE.min_phys(i)) THEN
               fail("max_phys must be > min_phys",exit_point=8888)
            ENDIF
         ENDDO
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (user_nsubs .LE. 0) THEN
               fail("no subs defined in user_defined assignment",exit_point=8888)
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
               fail("sub2proc must be allocated for user defined assignment",exit_point=8888)
            ENDIF
            DO i=1,user_nsubs
               IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                  fail("invalid processor specified in sub2proc",exit_point=8888)
               ENDIF
            ENDDO
         ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mktree_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mktree_d
#endif
