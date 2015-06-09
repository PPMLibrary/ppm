      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mkfield
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
      SUBROUTINE ppm_topo_mkfield_s(topoid,meshid,xp,Npart,decomp,assig,     &
      &          min_phys,max_phys,bcdef,ighostsize,cost,Nm,info,ndom,pcost, &
      &          user_minsub,user_maxsub,user_nsubs,user_sub2proc,Offset)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkfield_d(topoid,meshid,xp,Npart,decomp,assig,     &
      &          min_phys,max_phys,bcdef,ighostsize,cost,Nm,info,ndom,pcost, &
      &          user_minsub,user_maxsub,user_nsubs,user_sub2proc,Offset)
#endif
      !!! This routine is the topology creation routine for meshes.
      !!! In the user_defined case, the user must pass
      !!! existing subdomains and all of `user_minsub`, `user_maxsub`,
      !!! `cost`, `user_nsubs` must be provided. The sub boundaries
      !!! must align with mesh planes.
      !!! For all other decompositions, particles can be used
      !!! to guide them. For this, `Npart` must be > 0. For
      !!! Npart <= 0, the decomposition is purely mesh-based.
      !!! Subs are then mapped onto processors and
      !!! both the topology and the mesh definition is stored.
      !!!
      !!! [NOTE]
      !!! In the case of user-defined subs, a lot is
      !!! currently trusted. We may want to include further
      !!! checks in that case, e.g. if proper costs are provided,
      !!! idata and ndata actually make sense, the extents of
      !!! all subs are integer multiples of the mesh spacing, etc...

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_mesh_on_subs
      USE ppm_module_mesh_define
      USE ppm_module_decomp
      IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
      INTEGER,  PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
      INTEGER,  PARAMETER :: MK = ppm_kind_double
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
      INTEGER,                            INTENT(INOUT) :: meshid
      !!! Mesh identifier (user numbering) for default mesh (as defined by
      !!! Nm) on the topology. If <0, ppm will create a number internally
      !!! and return it here on exit.
      REAL(MK), DIMENSION(:,:),           POINTER       :: xp
      !!! Particle positions. If present, domain decomposition will be
      !!! guided by particle locations
      INTEGER,                            INTENT(IN   ) :: Npart
      !!! Number of particles. <=0 if no xp is given to guide the decomposition
      INTEGER,                            INTENT(IN   ) :: decomp
      !!! The valid decomposition types are:
      !!!
      !!! * ppm_param_decomp_bisection
      !!! * ppm_param_decomp_xpencil
      !!! * ppm_param_decomp_ypencil
      !!! * ppm_param_decomp_zpencil
      !!! * ppm_param_decomp_xy_slab
      !!! * ppm_param_decomp_xz_slab
      !!! * ppm_param_decomp_yz_slab
      !!! * ppm_param_decomp_cartesian
      !!! * ppm_param_decomp_cuboid
      !!! * ppm_param_decomp_user_defined
      INTEGER,                            INTENT(IN   ) :: assig
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
      REAL(MK), DIMENSION(:),             POINTER       :: min_phys
      !!! Minimum of physical extend of the computational domain (double)
      !!!
      !!! First index is ppm_dim
      REAL(MK), DIMENSION(:),             POINTER       :: max_phys
      !!! Maximum of physical extend of the computational domain (double)
      !!!
      !!! First index is ppm_dim
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
      INTEGER,  DIMENSION(:),             INTENT(IN   ) :: ighostsize
      !!! Size (width) of the ghost layer.
      REAL(MK), DIMENSION(:),             POINTER       :: cost
      !!! Estimated cost associated with subdomains. Either user-defined on
      !!! input or decomposition result on output. The cost of a subdomain
      !!! is given by its volume.
      INTEGER,  DIMENSION(:),             INTENT(IN   ) :: Nm
      !!! Number of *mesh points* (not cells) in each direction of the global
      !!! computational domain (including points ON its boundaries)
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.
      INTEGER,                  OPTIONAL, INTENT(IN   ) :: ndom
      !!! Number of subdomains requested. If not given, the number
      !!! of subs will be equal to the number of processors. This is only
      !!! relevent for decomp_cartesian and pencils without particles.
      REAL(MK), DIMENSION(:),   OPTIONAL, INTENT(IN   ) :: pcost
      !!! Array of length Npart which specifies the computational
      !!! cost attributed to each particle. If this is absent, a unity cost is
      !!! assumed for each particle.
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER       :: user_minsub
      !!! Mimimum of extension of subs.
      !!! Used if decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)
      !!! 2nd: subID
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER       :: user_maxsub
      !!! Maximum of extension of subs.
      !!! Used if decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)
      !!! 2nd: subID
      INTEGER,                  OPTIONAL, INTENT(IN   ) :: user_nsubs
      !!! Total number of subs on all processors.
      !!! parameter used when decomp is user defined
      INTEGER,  DIMENSION(:),   OPTIONAL, POINTER       :: user_sub2proc
      !!! Subdomain to processor assignment. index: subID (global)
      !!!  parameter used if assignment is user defined.
      REAL(MK), DIMENSION(:),   OPTIONAL, INTENT(IN   ) :: Offset
      !!! Offset in each dimension

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                                  :: t0
      REAL(MK)                                  :: parea,sarea,larea
      REAL(MK)                                  :: lmyeps,maxvar
      REAL(MK),              DIMENSION(ppm_dim) :: gsvec
      REAL(ppm_kind_double), DIMENSION(ppm_dim) :: h,Offst,righostsize
      REAL(MK),              DIMENSION(3,2)          :: weights
      REAL(MK),              DIMENSION(:,:), POINTER :: min_box => NULL()
      REAL(MK),              DIMENSION(:,:), POINTER :: max_box => NULL()
      REAL(MK),              DIMENSION(:,:), POINTER :: min_sub  => NULL()
      REAL(MK),              DIMENSION(:,:), POINTER :: max_sub  => NULL()

      INTEGER                          :: i,k
      INTEGER                          :: iopt,treetype,nbox
      INTEGER                          :: isub,minbox
      INTEGER, DIMENSION(1)            :: ldc
      INTEGER, DIMENSION(:,:), POINTER :: ineigh   => NULL()
      INTEGER, DIMENSION(:),   POINTER :: nneigh   => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: subs_bc  => NULL()
      INTEGER, DIMENSION(:),   POINTER :: nchld    => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: istart   => NULL()
      INTEGER, DIMENSION(:,:), POINTER :: ndata    => NULL()
      INTEGER, DIMENSION(:),   POINTER :: isublist => NULL()
      INTEGER, DIMENSION(:),   POINTER :: sub2proc => NULL()
      INTEGER                          :: nsubs
      INTEGER                          :: nsublist

      CHARACTER(LEN=ppm_char) :: msg
      CHARACTER(LEN=ppm_char) :: caller = 'ppm_topo_mkfield'

      LOGICAL, DIMENSION(ppm_dim) :: fixed
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

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
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

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

      IF (PRESENT(Offset)) THEN
#if    __KIND == __SINGLE_PRECISION
         Offst(1:ppm_dim)=REAL(Offset(1:ppm_dim),ppm_kind_double)
#elif  __KIND == __DOUBLE_PRECISION
         Offst(1:ppm_dim)=Offset(1:ppm_dim)
#endif
      ELSE
         Offst=0.0_ppm_kind_double
      ENDIF

      righostsize=REAL(ighostsize,ppm_kind_double)

      !-------------------------------------------------------------------------
      !  Compute grid spacing
      !-------------------------------------------------------------------------
      h(1) = (max_phys(1)-min_phys(1))/REAL(Nm(1)-1,ppm_kind_double)
      h(2) = (max_phys(2)-min_phys(2))/REAL(Nm(2)-1,ppm_kind_double)
      IF (ppm_dim.GT.2) THEN
         h(3) = (max_phys(3)-min_phys(3))/REAL(Nm(3)-1,ppm_kind_double)
      ENDIF

      !check for round-off problems and fix them if necessary
      DO k=1,ppm_dim
         DO WHILE (min_phys(k)+(Nm(k)-1)*h(k).LT.max_phys(k))
            h(k)=h(k)+EPSILON(h(k))
         ENDDO
      ENDDO
      check_true(<#ALL(min_phys(1:ppm_dim)+(Nm(1:ppm_dim)-1)*h(1:ppm_dim).GE.max_phys(1:ppm_dim))#>,"round-off problem in mesh creation")

      !-------------------------------------------------------------------------
      !  Cartesian (mesh-only) domain decomposition
      !-------------------------------------------------------------------------
      SELECT CASE (decomp)
      CASE (ppm_param_decomp_cartesian)
         IF (PRESENT(ndom)) THEN
            CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
            &    ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info,ndom)
         ELSE
            CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
            &    ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info)
         ENDIF
         or_fail('Cartesian decomposition failed')

      !-------------------------------------------------------------------------
      !  Recursive bisection. Can be particle-guided
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_bisection)
         !-------------------------------------------------------------------
         !  recursive bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype = ppm_param_tree_bin
         IF (Npart .GT. 0) THEN
            weights(1,1:2) = 0.5_MK    ! particles have 50% weight
            weights(2,1:2) = 0.5_MK    ! mesh has 50% weight
         ELSE
            weights(1,1:2) = 0.0_MK    ! only mesh has weight
            weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = REAL(righostsize(1:ppm_dim)*h(1:ppm_dim),MK)

         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
            CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
            &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
            &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Bisection decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Pencils
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_xpencil, &
      &     ppm_param_decomp_ypencil, &
      &     ppm_param_decomp_zpencil)
         IF (Npart.LT.1) THEN
             !------------------------------------------------------------------
             !  No particles: cartesian pencil decomposition
             !------------------------------------------------------------------
             IF (PRESENT(ndom)) THEN
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
                 &    decomp,min_sub,max_sub,nsubs,info,ndom)
             ELSE
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,         &
                 &    decomp,min_sub,max_sub,nsubs,info)
             ENDIF
             or_fail('Cartesian pencil decomposition failed')

         ELSE
             !------------------------------------------------------------------
             !  With particles: pencil quadrisection using ppm_tree
             !------------------------------------------------------------------
             ! build a quad tree, binary in 2d
             treetype           = ppm_param_tree_quad
             IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
             IF (Npart .GT. 0) THEN
                 weights(1,1:2) = 0.5_MK    ! particles have 50% weight
                 weights(2,1:2) = 0.5_MK    ! mesh has 50% weight
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
             gsvec(1) = REAL(righostsize(1)*h(1),MK)
             gsvec(2) = REAL(righostsize(2)*h(2),MK)
             IF (ppm_dim .GT. 2) THEN
                gsvec(3)= REAL(righostsize(3)*h(3),MK)
             ENDIF
             minbox = ppm_nproc
             IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
             ! build tree
             IF (PRESENT(pcost)) THEN
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
                 &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
                 &    min_box,max_box,nbox,nchld,info,pcost)
             ELSE
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
                 &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
                 &    min_box,max_box,nbox,nchld,info)
             ENDIF
             or_fail('Pencil decomposition failed')

             ! convert tree to subs
             CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,  &
             &    max_sub,nsubs,info)
             IF (info .NE. 0) GOTO 9999
         ENDIF

      !-------------------------------------------------------------------------
      !  Slabs
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_xy_slab, &
      &     ppm_param_decomp_xz_slab, &
      &     ppm_param_decomp_yz_slab)
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype = ppm_param_tree_bin
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
         gsvec(1)    = REAL(righostsize(1)*h(1),MK)
         gsvec(2)    = REAL(righostsize(2)*h(2),MK)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)= REAL(righostsize(3)*h(3),MK)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
             &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
             &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
             &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
             &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Slab decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999
      !-------------------------------------------------------------------------
      !  Cuboids
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_cuboid)
         !-------------------------------------------------------------------
         !  cuboid octasection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build an oct tree in 3d
         treetype         = ppm_param_tree_oct
         ! and a quad tree in 2d
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_quad
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK    ! mesh has 50% weight
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1)     = REAL(righostsize(1)*h(1),MK)
         gsvec(2)     = REAL(righostsize(2)*h(2),MK)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3) = REAL(righostsize(3)*h(3),MK)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
             &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
             &    min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
             &    minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
             &    min_box,max_box,nbox,nchld,info)
         ENDIF
         or_fail('Cuboid decomposition failed')

         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  User provides decomposition: Do nothing
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_user_defined)
         ! NOP
      !-------------------------------------------------------------------------
      !  Unknown decomposition type
      !-------------------------------------------------------------------------
      CASE DEFAULT
         stdout_f('(A,I5)',"Unknown decomposition type: ",decomp)
         fail(cbuf,ppm_error=ppm_error_fatal)

      END SELECT

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,nneigh,ineigh,gsvec,info)
      or_fail('Finding neighbors failed')

      !-------------------------------------------------------------------------
      !  Define meshes on the subs
      !-------------------------------------------------------------------------
      CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,min_sub, &
      &    max_sub,nsubs,istart,ndata,info,Offset=Offst)
      or_fail('Defining meshes failed')

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp.NE.ppm_param_decomp_user_defined) THEN
         IF (PRESENT(pcost)) THEN
            CALL ppm_topo_cost(xp,Npart,min_sub,max_sub, &
            &    nsubs,ndata,cost,info,pcost)
         ELSE
            CALL ppm_topo_cost(xp,Npart,min_sub,max_sub, &
            &    nsubs,ndata,cost,info)
         ENDIF
         or_fail('Computing costs failed')
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign the subdomains to processors
      !-------------------------------------------------------------------------
      SELECT CASE (assig)
      CASE (ppm_param_assign_internal)
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs, &
         &    sub2proc,isublist,nsublist,info)
         or_fail('Assigning subs to processors failed')

      CASE (ppm_param_assign_metis_cut,  &
      &     ppm_param_assign_metis_comm)
         gsvec(1:ppm_dim) = REAL(righostsize(1:ppm_dim)*h(1:ppm_dim),MK)

         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_phys,max_phys,min_sub,max_sub, &
         &    cost,gsvec,nneigh,ineigh,nsubs,sub2proc,isublist,     &
         &    nsublist,assig,info)
         or_fail('Assigning subs to processors using METIS failed')

      CASE (ppm_param_assign_user_defined)
         !--------------------------------------------------------------------
         !  user defined assignment
         !--------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         or_fail_alloc('list of local subs ISUBLIST')

         CALL ppm_alloc(sub2proc,ldc,iopt,info)
         or_fail_alloc('sub to processor assignment list SUB2PROC')

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
         stdout_f('(A,I5)',"Unknown assignment scheme: ",assig)
         fail(cbuf,ppm_error=ppm_error_fatal)

      END SELECT

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,subs_bc,info)
      or_fail('finding and defining the BC of the subs failed ')

      !-------------------------------------------------------------------------
      !  Store the topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(topoid,min_phys,max_phys,min_sub,max_sub,subs_bc, &
      &    sub2proc,nsubs,bcdef,0._MK,isublist,nsublist,nneigh,ineigh,info,decomp)
      or_fail('Storing topology failed')

      !-------------------------------------------------------------------------
      !  Create and store the new mesh internally
      !-------------------------------------------------------------------------
      CALL ppm_mesh_define(topoid,meshid,Nm,info,Offset=Offst,ghostsize=ighostsize)
      or_fail('Storing mesh definition failed')

      iopt = ppm_param_dealloc
      CALL ppm_alloc(ineigh,ldc,iopt,info)
      or_fail_dealloc("ineigh")

      CALL ppm_alloc(nneigh,ldc,iopt,info)
      or_fail_dealloc("nneigh")

      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      or_fail_dealloc("subs_bc")

      CALL ppm_alloc(istart,ldc,iopt,info)
      or_fail_dealloc("istart")

      CALL ppm_alloc(ndata,ldc,iopt,info)
      or_fail_dealloc("ndata")

      CALL ppm_alloc(isublist,ldc,iopt,info)
      or_fail_dealloc("isublist")

      IF (decomp.NE.ppm_param_decomp_cartesian .AND. &
      &   decomp.NE.ppm_param_decomp_user_defined) THEN
         CALL ppm_alloc(nchld,ldc,iopt,info)
         or_fail_dealloc('nchld')

         CALL ppm_alloc(min_box,ldc,iopt,info)
         or_fail_dealloc('min_box')

         CALL ppm_alloc(max_box,ldc,iopt,info)
         or_fail_dealloc('max_box')
      END IF

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
        IF (.NOT.ppm_initialized) THEN
           fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        DO i=1,ppm_dim
           IF (max_phys(i).LE.min_phys(i)) THEN
              fail('max_phys must be > min_phys',exit_point=8888)
           ENDIF
           IF (Nm(i).LT.2) THEN
              fail('Nm must be > 1 in all dimensions',exit_point=8888)
           ENDIF
           IF (ighostsize(i).LT. 0) THEN
              fail('ighostsize must be >= 0',exit_point=8888)
           ENDIF
        ENDDO
        IF (Npart.GT.0.AND.ASSOCIATED(xp)) THEN
           IF (SIZE(xp,2).LT.Npart) THEN
              fail('not enough particles contained in xp',exit_point=8888)
           ENDIF
           IF (SIZE(xp,1).LT.ppm_dim) THEN
              fail('leading dimension of xp too small',exit_point=8888)
           ENDIF
           IF (PRESENT(pcost)) THEN
              IF (SIZE(pcost,1).LT.Npart) THEN
                 fail('pcost does not contain costs for all particles',exit_point=8888)
              ENDIF
           ENDIF
        ENDIF
        IF (assig.EQ.ppm_param_assign_user_defined) THEN
           IF (decomp.NE.ppm_param_decomp_user_defined) THEN
              fail('decomp type is set to user_defined for this assignment',exit_point=8888)
           ENDIF
           IF (user_nsubs.LE.0) THEN
              fail('no subs defined in user_defined assignment',exit_point=8888)
           ENDIF
           IF (.NOT.ASSOCIATED(user_sub2proc)) THEN
              fail('sub2proc must be allocated for user defined assignment',exit_point=8888)
           ENDIF
           DO i=1,user_nsubs
              IF ((user_sub2proc(i).LT.0).OR.(user_sub2proc(i).GE.ppm_nproc)) THEN
                 fail('invalid processor specified in sub2proc',exit_point=8888)
              ENDIF
           ENDDO
        ENDIF
        IF (decomp.EQ.ppm_param_decomp_user_defined) THEN
           IF (user_nsubs.LE.0) THEN
              fail('no subs defined in user_defined decomposition',exit_point=8888)
           ENDIF
           IF ((.NOT.ASSOCIATED(user_minsub)).OR. &
           &   (.NOT.ASSOCIATED(user_maxsub))) THEN
               fail('min_sub/max_sub must be allocated for user def. decomp',exit_point=8888)
           ENDIF
           !-------------------------------------------------------------------
           !  Check that the user-defined subs add up to the whole
           !  computational domain.
           !  ONE COULD DO MORE TESTS HERE.
           !-------------------------------------------------------------------
           parea = (max_phys(1)-min_phys(1))*(max_phys(2)-min_phys(2))
           IF (ppm_dim.EQ.3) THEN
              parea = parea*(max_phys(3)-min_phys(3))
           ENDIF
           sarea = 0.0_MK
           DO i=1,user_nsubs
              larea = (user_maxsub(1,i)-user_minsub(1,i))* &
              &       (user_maxsub(2,i)-user_minsub(2,i))
              IF (ppm_dim.EQ.3) THEN
                 larea = larea * (user_maxsub(3,i)-user_minsub(3,i))
              END IF
              sarea = sarea + larea
           ENDDO
           IF (ABS(sarea-parea)/parea.GE.lmyeps) THEN
              !----------------------------------------------------------------
              !  Mismatch!
              !----------------------------------------------------------------
              fail('faulty subdomains defined',exit_point=8888)
           ENDIF
        ENDIF

        IF (decomp.EQ.ppm_param_decomp_zpencil.AND.ppm_dim.LT.3) THEN
           fail('Cannot make z pencils in 2D!',exit_point=8888)
        ENDIF
        IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
           fail('Cannot make x-z slabs in 2D!',exit_point=8888)
        ENDIF
        IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
           fail('Cannot make y-z slabs in 2D!',exit_point=8888)
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
      8888 CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_d
#endif
