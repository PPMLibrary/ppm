      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_topo_mkpart
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
      SUBROUTINE ppm_topo_mkpart_s(topoid,xp,Npart,decomp,assig, &
      &          min_phys,max_phys,bcdef,ghostsize,cost,info,    &
      &          pcost,user_minsub,user_maxsub,user_nsubs,       &
      &          user_sub2proc)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkpart_d(topoid,xp,Npart,decomp,assig, &
      &          min_phys,max_phys,bcdef,ghostsize,cost,info,    &
      &          pcost,user_minsub,user_maxsub,user_nsubs,       &
      &          user_sub2proc)
#endif
      !!! This routine is the topology creation routine for particles.
      !!! It performs the decomposition of the physical space based on the
      !!! position of the particles. The subdomains are mapped
      !!! onto the processors and a neighbour list is created.
      !!!
      !!! The decomposition is based on the particle positions
      !!! If user_nsubs is greater than zero, the subdomains
      !!! described by `(min_sub,max_sub)` given by the user is used.
      !!! If user_nsubs is less than one, the subdomains are
      !!! found using the decomposition defined by the option
      !!! decomp.
      !!!
      !!! If user_nsubs on input is greater than zero, the assigment
      !!! of subdomain to the processors is used (as described
      !!! in sub2proc) or else
      !!! the assigment of subdomain to the processors will be
      !!! performed to assign as best as possible an equal number
      !!! of particles/grid points and at the same time
      !!! minimizing communication (see Farhat 1988, JCP 28(5),
      !!! pages: 579 -- 602).
      !!!
      !!! If topoid is 0 a new topology is created. Else the topoid is used to
      !!! locate an existing topology and overwrite it
      !!!
      !!! [NOTE]
      !!! In the case of user-defined subs, a lot is
      !!! currently trusted. We may want to include further
      !!! checks in that case, e.g. if proper costs are provided.
      !!!
      !!! .References
      !!! *************************************************************
      !!! - C. Farhat, A Simple and Efficiant Automatic FEM Decomposer.
      !!!   Computers and Structures. 1988. 28(5), 579-602.
      !!! *************************************************************

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_mpi
      USE ppm_module_decomp
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,                            INTENT(INOUT) :: topoid
      !!! ID number identifying the topology.
      !!! If topoid == 0 on input a new topology is created and the new topoid
      !!! is returned here, else the indicated toplogy is replaced.
      !!!
      !!! [CAUTION]
      !!! *SEMANTICS CHANGED:* `topoid = 0` is *not anymore* reserved for the
      !!! ring topology (null decomposition). "Ring topologies" are non
      !!! geometric and need no setup. The user can perform ppm ring shift
      !!! operations without having to first define a topology.
      REAL(MK), DIMENSION(:,:),           POINTER       :: xp
      !!! Position of particles
      INTEGER,                            INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER,                            INTENT(IN   ) :: decomp
      !!! The type of decomposition. One of:
      !!!
      !!! *  ppm_param_decomp_pruned_cell
      !!! *  ppm_param_decomp_tree
      !!! *  ppm_param_decomp_bisection
      !!! *  ppm_param_decomp_xpencil
      !!! *  ppm_param_decomp_ypencil
      !!! *  ppm_param_decomp_zpencil
      !!! *  ppm_param_decomp_xy_slab
      !!! *  ppm_param_decomp_xz_slab
      !!! *  ppm_param_decomp_yz_slab
      !!! *  ppm_param_decomp_cuboid
      !!! *  ppm_param_decomp_user_defined
      !!!
      !!! [NOTE]
      !!! There is no more ppm_param_decomp_null. If the user does not
      !!! want to define a geometric decomposition there is no need for a
      !!! topology.
      INTEGER,                            INTENT(IN   ) :: assig
      !!! The type of subdomain-to-processor assignment. One of:
      !!!
      !!! *  ppm_param_assign_internal
      !!! *  ppm_param_assign_metis_cut
      !!! *  ppm_param_assign_metis_comm
      !!! *  ppm_param_assign_user_defined
      !!!
      !!! *nodal* and *dual* assignments  use the external library METIS
      !!! and are only available if ppm was compiled with METIS support.


      REAL(MK), DIMENSION(:),             POINTER       :: min_phys
      !!! Minimum of physical extend of the computational domain
      !!!
      !!! first index is ppm_dim
      REAL(MK), DIMENSION(:),             POINTER       :: max_phys
      !!! Maximum of physical extend of the computational domain
      !!!
      !!! first index is ppm_dim
      INTEGER,  DIMENSION(:),             INTENT(IN   ) :: bcdef
      !!! Boundary conditions for the topology
      !!!
      !!! NOTE: first index is 1-6 (each of the faces)
      !!! - west  : 1
      !!! - east  : 2
      !!! - south : 3
      !!! - north : 4
      !!! - bottom: 5
      !!! - top   : 6
      REAL(MK),                           INTENT(IN   ) :: ghostsize
      !!! Size (width) of the ghost layer.
      !!!
      !!! TIP: If the particles are moving and you intend to use the ppm
      !!! neighbor lists,  this argument must be equal the particle interaction
      !!! radius plus the skin parameter
      REAL(MK), DIMENSION(:),             POINTER       :: cost
      !!! Estimated cost associated with subdomains. Either user-defined on
      !!! input or decomposition result on output.
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success

      REAL(MK), DIMENSION(:),   OPTIONAL, INTENT(IN   ) :: pcost
      !!! Array of length Npart which specifies the computational
      !!! cost attributed to each particle. If this is absent, a unity cost is
      !!! assumed for each particle.

      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER       :: user_minsub
      !!! Mimimum of extension of subs.
      !!! Used when decomp is user defined
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER       :: user_maxsub
      !!! maximum of extension of subs
      !!! Used when decomp is user defined
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      INTEGER,                  OPTIONAL                :: user_nsubs
      !!! Total number of subs on all processors.
      !!! Used when decomp is user defined
      INTEGER,  DIMENSION(:),   OPTIONAL, POINTER       :: user_sub2proc
      !!! subdomain to processor assignment. index: subID (global)
      !!! Used when assign is user defined
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                          :: t0
      REAL(MK)                          :: parea,sarea,larea,lmyeps
      REAL(MK), DIMENSION(3,2)          :: weights
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(:,:), POINTER :: min_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: min_sub => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_sub => NULL()

      INTEGER                          :: i,Ntot,iopt,treetype
      INTEGER                          :: nbox,isub
      INTEGER, DIMENSION(3,1)          :: nnodes
      INTEGER, DIMENSION(3)            :: ldc
      INTEGER, DIMENSION(ppm_dim)      :: Nm
      INTEGER                          :: nsublist
      INTEGER                          :: nsubs
      INTEGER, DIMENSION(:,:), POINTER :: ineigh   => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: subs_bc  => NULL()
      INTEGER, DIMENSION(  :), POINTER :: nneigh   => NULL()
      INTEGER, DIMENSION(  :), POINTER :: nchld    => NULL()
      INTEGER, DIMENSION(:  ), POINTER :: sub2proc => NULL()
      INTEGER, DIMENSION(  :), POINTER :: isublist => NULL()

      CHARACTER(ppm_char) :: cbuf
      CHARACTER(ppm_char) :: caller='ppm_topo_mkpart'

      LOGICAL, DIMENSION(ppm_dim) :: fixed
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Check that we have particles
      !-------------------------------------------------------------------------
#ifdef __MPI
      CALL MPI_AllReduce(Npart,Ntot,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
      or_fail_MPI('MPI_AllReduce failed')

      IF (Ntot.LT.1) THEN
         fail('No particles in domain',ppm_error=ppm_error_notice)
      ENDIF
#else
      IF (Npart.LT.1) THEN
        fail('No particle on this processor',ppm_error=ppm_error_notice)
      ENDIF
#endif

      ! If the user defined nsubs then use those
      IF (PRESENT(user_nsubs)) THEN
         nsubs = user_nsubs
      ELSE
         nsubs = 0
      ENDIF
      IF (PRESENT(user_minsub)) THEN
         min_sub => user_minsub
      ENDIF
      IF (PRESENT(user_maxsub)) THEN
         max_sub => user_maxsub
      ENDIF
      IF (PRESENT(user_sub2proc)) THEN
         sub2proc => user_sub2proc
      ENDIF

      !----------------------------------------------------------------------
      !  Dummy argument for non-existing mesh
      !----------------------------------------------------------------------
      nnodes(1:3,1) = 0

      !----------------------------------------------------------------------
      !  Perform the decomposition using various techniques
      !----------------------------------------------------------------------
      SELECT CASE (decomp)
      CASE (ppm_param_decomp_pruned_cell)
         !-------------------------------------------------------------------
         !  a pruned cell index list
         !-------------------------------------------------------------------
         IF (PRESENT(pcost)) THEN
            CALL ppm_decomp_pruned_cell(xp,Npart,min_phys,max_phys, &
            &    ghostsize,min_sub,max_sub,nsubs,info,pcost)
         ELSE
            CALL ppm_decomp_pruned_cell(xp,Npart,min_phys,max_phys, &
            &    ghostsize,min_sub,max_sub,nsubs,info)
         ENDIF
         or_fail('Pruned cell decomposition failed')

      CASE (ppm_param_decomp_tree)
         !-------------------------------------------------------------------
         !  a tree data structure; use a default maxvariance of 10 pct
         !-------------------------------------------------------------------
         IF (PRESENT(pcost)) THEN
            CALL ppm_decomp_tree(xp,Npart,min_phys,max_phys,ghostsize, &
            &    0.1_MK,min_sub,max_sub,nsubs,info,pcost)
         ELSE
            CALL ppm_decomp_tree(xp,Npart,min_phys,max_phys,ghostsize, &
            &    0.1_MK,min_sub,max_sub,nsubs,info)
         ENDIF
         or_fail('Tree decomposition failed')

      CASE (ppm_param_decomp_bisection)
         !-------------------------------------------------------------------
         !  recursive bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype = ppm_param_tree_bin
         ! particles have unit weight
         weights(1,1) = 1.0_MK
         weights(1,2) = 1.0_MK
         ! mesh and geometry are not considered
         weights(2,1) = 0.0_MK
         weights(2,2) = 0.0_MK
         weights(3,1) = 0.0_MK
         weights(3,2) = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! no mesh
         Nm(1:ppm_dim) = 0
         ! build tree
         IF (PRESENT(pcost)) THEN
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Bisection decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      CASE (ppm_param_decomp_xpencil, &
      &     ppm_param_decomp_ypencil, &
      &     ppm_param_decomp_zpencil)
         IF ((decomp.EQ.ppm_param_decomp_zpencil).AND.(ppm_dim.LT.3)) THEN
            fail('Cannot make z pencils in 2D!',ppm_err_argument)
         ENDIF
         !-------------------------------------------------------------------
         !  pencil quadrisection using the general ppm_tree
         !-------------------------------------------------------------------
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
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Pencil decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      CASE (ppm_param_decomp_xy_slab, &
      &     ppm_param_decomp_xz_slab, &
      &     ppm_param_decomp_yz_slab)
         IF ((decomp.EQ.ppm_param_decomp_xz_slab).AND.(ppm_dim.LT.3)) THEN
            fail('Cannot make x-z slabs in 2D!',ppm_err_argument)
         ENDIF
         IF ((decomp.EQ.ppm_param_decomp_yz_slab).AND.(ppm_dim.LT.3)) THEN
            fail('Cannot make y-z slabs in 2D!',ppm_err_argument)
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
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,      &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights, &
            &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,      &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights, &
            &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Slab decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      CASE (ppm_param_decomp_cuboid)
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
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Cuboid decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      CASE (ppm_param_decomp_user_defined)
         ! Do nothing. Just take the stuff from the user and trust the guy.

      CASE DEFAULT
         !-------------------------------------------------------------------
         !  unknown decomposition
         !-------------------------------------------------------------------
         WRITE(cbuf,'(A,I5)') 'Unknown decomposition type: ',decomp
         fail(cbuf,ppm_err_argument)

      END SELECT

      !-------------------------------------------------------------------------
      !  Find the neighbours to the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,nneigh,ineigh,gsvec,info)
      or_fail('Finding neighbors failed')

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
         IF (PRESENT(pcost)) THEN
            CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,nnodes,cost, &
            &    info,pcost)
         ELSE
            CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,nnodes,cost,info)
         ENDIF
         or_fail('Computing costs failed')
      ENDIF

      !----------------------------------------------------------------------
      !  Define the topology (assign the subdomains to processors)
      !----------------------------------------------------------------------
      SELECT CASE (assig)
      CASE (ppm_param_assign_internal)
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, &
         &    isublist,nsublist,info)
         or_fail('Assigning subs to processors failed')

      CASE (ppm_param_assign_metis_cut, &
      &     ppm_param_assign_metis_comm)
         gsvec(1:ppm_dim) = ghostsize

         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_phys,max_phys,min_sub,max_sub, &
         &    cost,gsvec,nneigh,ineigh,nsubs,sub2proc,isublist,     &
         &    nsublist,assig,info)
         or_fail('Assigning subs to processors using METIS failed')

      CASE (ppm_param_assign_user_defined)
         !-------------------------------------------------------------------
         !  user defined assignment
         !-------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         or_fail_alloc('list of local subs ISUBLIST')

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
         WRITE(cbuf,'(A,I5)') 'Unknown assignment scheme: ',assig
         fail(cbuf,ppm_err_argument)

      END SELECT

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,subs_bc,info)
      or_fail('finding and defining the BC of the subs failed ')

      ! isublist:  my subs (ppm_rank)
      ! sub2proc:  subs -> processor map
      ! ineigh:    neighbors of all subs
      !-------------------------------------------------------------------------
      !  Store the current topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(topoid,min_phys,max_phys,min_sub,max_sub,  &
      &    subs_bc,sub2proc,nsubs,bcdef,ghostsize,isublist,nsublist, &
      &    nneigh,ineigh,info,decomp)
      or_fail('Storing topology failed')

      iopt = ppm_param_dealloc
      CALL ppm_alloc(nchld,ldc,iopt,info)
      or_fail_dealloc('nchld')

      CALL ppm_alloc(min_box,ldc,iopt,info)
      or_fail_dealloc('min_box')

      CALL ppm_alloc(max_box,ldc,iopt,info)
      or_fail_dealloc('max_box')

      CALL ppm_alloc(ineigh,ldc,iopt,info)
      or_fail_dealloc("ineigh")

      CALL ppm_alloc(nneigh,ldc,iopt,info)
      or_fail_dealloc("nneigh")

      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      or_fail_dealloc("subs_bc")

      CALL ppm_alloc(isublist,ldc,iopt,info)
      or_fail_dealloc("isublist")

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
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
             fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
          ENDIF ! test ppm_initialized
          IF (Npart .LT. 0) THEN
             fail('Npart must not be negative!',exit_point=8888)
          ENDIF ! test Npart > 0
          IF (.NOT.ASSOCIATED(xp)) THEN
             fail('xp is not allocated!',exit_point=8888)
          ENDIF ! test xp associated
          IF (SIZE(xp,1) .LT. ppm_dim) THEN
             fail('leading dimension of xp too small!',exit_point=8888)
          ENDIF ! test xp enough dims
          IF (SIZE(xp,2) .LT. Npart) THEN
             fail('not enough particles contained in xp!',exit_point=8888)
          ENDIF ! test xp >= Npart
          IF (topoid.NE.0) THEN
             IF ((topoid.GT.SIZE(ppm_topo)).OR.(topoid.LT.1)) THEN
                fail('topoid indexing outside ppm_topo!',exit_point=8888)
             ENDIF
          ENDIF
          IF (PRESENT(pcost)) THEN
             IF (SIZE(pcost,1) .LT. Npart) THEN
                fail('pcost must be of at least length Npart!',exit_point=8888)
             ENDIF ! test pcost >= Npart
          ENDIF ! test pcost present
          IF (ghostsize.LT.0.0_MK) THEN
             fail('ghostsize must not be negative!',exit_point=8888)
          ENDIF ! test ghostsize
          DO i=1,ppm_dim
             IF (max_phys(i) .LE. min_phys(i)) THEN
                fail('max_phys must be > min_phys!',exit_point=8888)
             ENDIF ! test valid dims
          ENDDO ! for each dimension
          IF (assig .EQ. ppm_param_assign_user_defined) THEN
             IF (decomp .NE. ppm_param_decomp_user_defined) THEN
                fail('decomp type set to user_defined for this assignment!',exit_point=8888)
             ENDIF ! test decomp user_defined
             IF (.NOT.PRESENT(user_nsubs)) THEN
                fail('user_nsubs must be provided if assignment user defined!',exit_point=8888)
             ENDIF ! test user_nsubs present
             IF (user_nsubs .LE. 0) THEN
                fail('user_nsubs has to be > 0!',exit_point=8888)
             ENDIF ! test user_nsubs > 0
             IF (.NOT.PRESENT(user_minsub)) THEN
                fail('min_sub must be provided if assignment user defined!',exit_point=8888)
             ENDIF ! test min_subs present
             IF (.NOT.ASSOCIATED(user_minsub)) THEN
                fail('min_sub has to be associated!',exit_point=8888)
             ENDIF ! test min_sub asoc
             IF (SIZE(user_minsub,2).LT.user_nsubs) THEN
                fail('min_sub has to be of length at least user_nsubs!',exit_point=8888)
             ENDIF ! test size of min_sub
             IF (.NOT.PRESENT(user_maxsub)) THEN
                fail('max_sub must be provided if assignment user defined!',exit_point=8888)
             ENDIF ! test max_subs present
             IF (.NOT.ASSOCIATED(user_maxsub)) THEN
                fail('max_sub has to be associated!',exit_point=8888)
             ENDIF ! test max_sub asoc
             IF (SIZE(user_maxsub,2).LT.user_nsubs) THEN
                fail('max_sub has to be of length at least user_nsubs!',exit_point=8888)
             ENDIF ! test size of max_sub
             IF (.NOT.PRESENT(user_sub2proc)) THEN
                fail('sub2proc must be provided if assignment user defined!',exit_point=8888)
             ENDIF ! test sub2proc present
             IF (.NOT.ASSOCIATED(user_sub2proc)) THEN
                fail('sub2proc has to be associated!',exit_point=8888)
             ENDIF ! test sub2proc asoc
             IF (SIZE(user_sub2proc).LT.user_nsubs) THEN
                fail('sub2proc must be at least length user_nsubs!',exit_point=8888)
             ENDIF ! test size of sub2proc
             DO i=1,user_nsubs
                IF ((user_sub2proc(i).LT.0).OR.(user_sub2proc(i).GE.ppm_nproc)) THEN
                   fail('invalid processor specified in sub2proc!',exit_point=8888)
                ENDIF ! test if processor assign valid
             ENDDO ! for each sub
          ENDIF ! do user defined assignment checks
          IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
             IF (user_nsubs .LE. 0) THEN
                fail('no subs defined in user_defined decomposition!',exit_point=8888)
             ENDIF ! test user_nsubs
             !-------------------------------------------------------------
             !  Check that the user-defined subs add up to the whole
             !  computational domain.
             !  TODO: ONE COULD DO MORE TESTS HERE.
             !-------------------------------------------------------------
             parea = (max_phys(1)-min_phys(1))*(max_phys(2)-min_phys(2))
             IF (ppm_dim.EQ.3) THEN
                parea = parea*(max_phys(3)-min_phys(3))
             ENDIF
             sarea = 0.0_MK
             DO i=1,user_nsubs
                larea = (user_maxsub(1,i)-user_minsub(1,i))*(user_maxsub(2,i)-user_minsub(2,i))
                IF (ppm_dim.EQ.3) THEN
                   larea = larea*(user_maxsub(3,i)-user_minsub(3,i))
                END IF
                sarea = sarea + larea
             ENDDO
             IF (ABS(sarea-parea)/parea.GE.lmyeps) THEN
                !---------------------------------------------------------
                !  Mismatch!
                !---------------------------------------------------------
                fail('faulty subdomains defined!',exit_point=8888)
             ENDIF
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkpart_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkpart_d
#endif
