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
      SUBROUTINE ppm_topo_mkfield_s(topoid,meshid,xp,Npart,decomp,assig,      &
     &               min_phys,max_phys,bcdef,ighostsize,cost,                  &
     &               Nm,info,ndom,pcost,user_minsub,user_maxsub, &
     &               user_nsubs,user_sub2proc)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkfield_d(topoid,meshid,xp,Npart,decomp,assig,      &
     &               min_phys,max_phys,bcdef,ighostsize,cost,                  &
     &               Nm,info,ndom,pcost,user_minsub,user_maxsub, &
     &               user_nsubs,user_sub2proc)
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
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_define_subs_bc
      USE ppm_module_mesh_on_subs
      USE ppm_module_mesh_alloc
      USE ppm_module_mesh_store
      USE ppm_module_topo_subs2proc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_find_neigh
      USE ppm_module_decomp
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
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(INOUT) :: topoid
      !!! ID number identifying the topology.                                  +
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
      !!! Number of particles. <=0 if no xp is given to guide the decomposition
      INTEGER                 , INTENT(IN   ) :: decomp
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
      INTEGER                 , INTENT(IN   ) :: assig
      !!! The type of subdomain-to-processor assignment. One of:
      !!!
      !!! *  ppm_param_assign_internal
      !!! *  ppm_param_assign_nodal_cut
      !!! *  ppm_param_assign_nodal_comm
      !!! *  ppm_param_assign_dual_cut
      !!! *  ppm_param_assign_dual_comm
      !!! *  ppm_param_assign_user_defined
      !!!
      !!! [NOTE]
      !!! The latter uses the external library METIS and is only
      !!! available if ppm was compiled with METIS support.
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_phys
      !!! Minimum of physical extend of the computational domain (double)
      !!!
      !!! First index is ppm_dim
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: max_phys
      !!! Maximum of physical extend of the computational domain (double)
      !!!
      !!! First index is ppm_dim
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
      INTEGER,  DIMENSION(:  ), INTENT(IN   ) :: ighostsize
      !!! Size (width) of the ghost layer.
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      !!! Estimated cost associated with subdomains. Either user-defined on
      !!! input or decomposition result on output. The cost of a subdomain
      !!! is given by its volume.
      INTEGER                 , INTENT(INOUT) :: meshid
      !!! Mesh identifier (user numbering) for default mesh (as defined by
      !!! Nm) on the topology. If <0, ppm will create a number internally
      !!! and return it here on exit.
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      !!! Number of *mesh points* (not cells) in each direction of the global
      !!! computational domain (including points ON its boundaries)
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success.
      INTEGER , OPTIONAL      , INTENT(IN   ) :: ndom
      !!! Number of subdomains requested. If not given, the number
      !!! of subs will be equal to the number of processors. This is only
      !!! relevent for decomp_cartesian and pencils without particles.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Array of length Npart which specifies the computational
      !!! cost attributed to each particle. If this is absent, a unity cost is
      !!! assumed for each particle.
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER :: user_minsub
      !!! Mimimum of extension of subs.
      !!! Used if decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      REAL(MK), DIMENSION(:,:), OPTIONAL, POINTER :: user_maxsub
      !!! Maximum of extension of subs.
      !!! Used if decomp is user defined.
      !!!
      !!! 1st index: x,y,(z)                                                   +
      !!! 2nd: subID
      INTEGER                 , OPTIONAL          :: user_nsubs
      !!! Total number of subs on all processors.
      !!! parameter used when decomp is user defined
      INTEGER, DIMENSION(:  ),  OPTIONAL, POINTER :: user_sub2proc
      !!! Subdomain to processor assignment. index: subID (global)
      !!!  parameter used if assignment is user defined.


      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                           :: i,iopt,treetype,nbox
      INTEGER                           :: isub,minbox
      INTEGER, DIMENSION(1)             :: ldc
      INTEGER, DIMENSION(:,:), POINTER  :: ineigh  => NULL()
      INTEGER, DIMENSION(:,:), POINTER  :: subs_bc => NULL()
      INTEGER, DIMENSION(:  ), POINTER  :: nneigh  => NULL()
      INTEGER, DIMENSION(:  ), POINTER  :: nchld   => NULL()
      INTEGER, DIMENSION(:,:), POINTER  :: istart  => NULL()
      INTEGER, DIMENSION(:,:), POINTER  :: ndata   => NULL()
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps,maxvar
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec,meshdx
      LOGICAL , DIMENSION(ppm_dim)      :: fixed
      REAL(MK), DIMENSION(3,2)          :: weights
      REAL(MK), DIMENSION(:,:), POINTER :: min_box => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box => NULL()
      CHARACTER(LEN=ppm_char)           :: mesg
      INTEGER                           :: nsublist
      INTEGER , DIMENSION(  :), POINTER :: isublist => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: min_sub
      REAL(MK), DIMENSION(:,:), POINTER :: max_sub
      INTEGER                           :: nsubs
      INTEGER, DIMENSION(:  ), POINTER  :: sub2proc
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_mkfield',t0,info)
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
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      ! If the user defined nsubs then use those
      IF (PRESENT(user_nsubs)) THEN
          nsubs = user_nsubs
      ELSE
          nsubs = 0
      ENDIF
      IF (PRESENT(user_minsub)) THEN
          min_sub => user_minsub
      ELSE
          NULLIFY(min_sub)
      ENDIF
      IF (PRESENT(user_maxsub)) THEN
          max_sub => user_maxsub
      ELSE
          NULLIFY(max_sub)
      ENDIF
      IF (PRESENT(user_sub2proc)) THEN
          sub2proc => user_sub2proc
      ELSE
          NULLIFY(sub2proc)
      ENDIF
      !-------------------------------------------------------------------------
      !  Compute grid spacing
      !-------------------------------------------------------------------------
      meshdx(1) = (max_phys(1)-min_phys(1))/REAL(Nm(1)-1,MK)
      meshdx(2) = (max_phys(2)-min_phys(2))/REAL(Nm(2)-1,MK)
      IF (ppm_dim .GT. 2) THEN
          meshdx(3) = (max_phys(3)-min_phys(3))/REAL(Nm(3)-1,MK)
      ENDIF

      !-------------------------------------------------------------------------
      !  Cartesian (mesh-only) domain decomposition
      !-------------------------------------------------------------------------
      IF (decomp .EQ. ppm_param_decomp_cartesian) THEN
          IF (PRESENT(ndom)) THEN
              CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &            ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info,ndom)
          ELSE
              CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &            ppm_param_decomp_cuboid,min_sub,max_sub,nsubs,info)
          ENDIF
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &            'Cartesian decomposition failed',__LINE__,info)
              GOTO 9999
          ENDIF

      !-------------------------------------------------------------------------
      !  Recursive bisection. Can be particle-guided
      !-------------------------------------------------------------------------
      ELSEIF (decomp.EQ.ppm_param_decomp_bisection) THEN
         !-------------------------------------------------------------------
         !  recursive bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype           = ppm_param_tree_bin
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1)         = REAL(ighostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ighostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)     = REAL(ighostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,       &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &           'Bisection decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Pencils
      !-------------------------------------------------------------------------
      ELSEIF ((decomp .EQ. ppm_param_decomp_xpencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_ypencil) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_zpencil)) THEN
         IF (decomp.EQ.ppm_param_decomp_zpencil.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make z pencils in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (Npart .LT. 1) THEN
             !------------------------------------------------------------------
             !  No particles: cartesian pencil decomposition
             !------------------------------------------------------------------
             IF (PRESENT(ndom)) THEN
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,  &
     &               decomp,min_sub,max_sub,nsubs,info,ndom)
             ELSE
                 CALL ppm_decomp_cartesian(Nm,min_phys,max_phys,         &
     &               decomp,min_sub,max_sub,nsubs,info)
             ENDIF
             IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',   &
     &               'Cartesian pencil decomposition failed',__LINE__,info)
                 GOTO 9999
             ENDIF
         ELSE
             !------------------------------------------------------------------
             !  With particles: pencil quadrisection using ppm_tree
             !------------------------------------------------------------------
             ! build a quad tree, binary in 2d
             treetype           = ppm_param_tree_quad
             IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
             IF (Npart .GT. 0) THEN
                 weights(1,1:2) = 0.5_MK    ! particles have 50% weight
                 weights(2,1:2) = 0.5_MK
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
             gsvec(1)         = REAL(ighostsize(1),MK)*meshdx(1)
             gsvec(2)         = REAL(ighostsize(2),MK)*meshdx(2)
             IF (ppm_dim .GT. 2) THEN
                 gsvec(3)     = REAL(ighostsize(3),MK)*meshdx(3)
             ENDIF
             minbox = ppm_nproc
             IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
             ! build tree
             IF (PRESENT(pcost)) THEN
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
     &               minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
     &               min_box,max_box,nbox,nchld,info,pcost)
             ELSE
                 CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,   &
     &               minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,  &
     &               min_box,max_box,nbox,nchld,info)
             ENDIF
             IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',   &
     &               'Pencil decomposition failed',__LINE__,info)
                 GOTO 9999
             ENDIF
             ! convert tree to subs
             CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,  &
     &           max_sub,nsubs,info)
             IF (info .NE. 0) GOTO 9999
         ENDIF

      !-------------------------------------------------------------------------
      !  Slabs
      !-------------------------------------------------------------------------
      ELSEIF ((decomp .EQ. ppm_param_decomp_xy_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_xz_slab) .OR.   &
     &        (decomp .EQ. ppm_param_decomp_yz_slab)) THEN
         IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make x-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',  &
     &           'Cannot make y-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
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
         gsvec(1)         = REAL(ighostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ighostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)     = REAL(ighostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',    &
     &           'Slab decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999
      !-------------------------------------------------------------------------
      !  Cuboids
      !-------------------------------------------------------------------------
      ELSEIF (decomp .EQ. ppm_param_decomp_cuboid) THEN
         !-------------------------------------------------------------------
         !  cuboid octasection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build an oct tree in 3d
         treetype         = ppm_param_tree_oct
         ! and a quad tree in 2d
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_quad
         IF (Npart .GT. 0) THEN
             weights(1,1:2) = 0.5_MK    ! particles have 50% weight
             weights(2,1:2) = 0.5_MK
         ELSE
             weights(1,1:2) = 0.0_MK    ! only mesh has weight
             weights(2,1:2) = 1.0_MK
         ENDIF
         ! geometry has no weight
         weights(3,1:2)     = 0.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1)         = REAL(ighostsize(1),MK)*meshdx(1)
         gsvec(2)         = REAL(ighostsize(2),MK)*meshdx(2)
         IF (ppm_dim .GT. 2) THEN
             gsvec(3)         = REAL(ighostsize(3),MK)*meshdx(3)
         ENDIF
         minbox = ppm_nproc
         IF (PRESENT(ndom)) minbox = MAX(ppm_nproc,ndom)
         ! build tree
         IF (PRESENT(pcost)) THEN
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info,pcost)
         ELSE
             CALL ppm_tree(xp,Npart,Nm,min_phys,max_phys,treetype,    &
     &           minbox,.FALSE.,gsvec,maxvar,-1.0_MK,fixed,weights,   &
     &           min_box,max_box,nbox,nchld,info)
         ENDIF
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',    &
     &           'Cuboid decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
     &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  User provides decomposition: Do nothing
      !-------------------------------------------------------------------------
      ELSEIF (decomp .EQ. ppm_param_decomp_user_defined) THEN
          ! NOP

      !-------------------------------------------------------------------------
      !  Unknown decomposition type
      !-------------------------------------------------------------------------
      ELSE
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown decomposition type: ',decomp
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',         &
     &       mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
     &                    min_sub,max_sub,nsubs,nneigh,ineigh,gsvec,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',      &
     &        'Finding neighbors failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Define meshes on the subs
      !-------------------------------------------------------------------------
      CALL ppm_mesh_on_subs(Nm,min_phys,max_phys,min_sub,max_sub,nsubs,&
     &                      istart,ndata,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &       'Defining meshes failed',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          IF (PRESENT(pcost)) THEN
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,  &
     &                            info,pcost)
          ELSE
              CALL ppm_topo_cost(xp,Npart,min_sub,max_sub,nsubs,ndata,cost,info)
          ENDIF
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &           'Computing costs failed',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign the subdomains to processors
      !-------------------------------------------------------------------------
      IF     (assig .EQ. ppm_param_assign_internal) THEN
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, &
     &       isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',   &
     &         'Assigning subs to processors failed',__LINE__,info)
            GOTO 9999
         ENDIF
      ELSEIF (assig .EQ. ppm_param_assign_nodal_cut .OR.    &
     &        assig .EQ. ppm_param_assign_nodal_comm .OR.   &
     &        assig .EQ. ppm_param_assign_dual_cut .OR.     &
     &        assig .EQ. ppm_param_assign_dual_comm) THEN
         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_sub,max_sub,nneigh,ineigh,cost,nsubs,&
     &       assig,sub2proc,isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &         'Assigning subs to processors using METIS failed',__LINE__,&
     &          info)
            GOTO 9999
         ENDIF
     ELSEIF (assig .EQ. ppm_param_assign_user_defined) THEN
         !--------------------------------------------------------------------
         !  user defined assignment
         !--------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mkfield',   &
     &           'list of local subs ISUBLIST',__LINE__,info)
                 GOTO 9999
         ENDIF
         CALL ppm_alloc(sub2proc,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mkfield',   &
     &           'sub to processor assignment list SUB2PROC',__LINE__,info)
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
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',   &
     &       mesg,__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      NULLIFY(subs_bc)
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef,min_sub,max_sub,nsubs,subs_bc,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',  &
     &       'finding and defining the BC of the subs failed ',__LINE__,&
     &        info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the topology
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(topoid,min_phys,max_phys,min_sub,max_sub,subs_bc, &
      &     sub2proc,nsubs,bcdef,0._MK,isublist,nsublist,nneigh,ineigh,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',      &
     &        'Storing topology failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new mesh internally
      !-------------------------------------------------------------------------
      CALL ppm_mesh_store(topoid,meshid,ndata,istart,Nm,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkfield',      &
     &        'Storing mesh definition failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ineigh,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','ineigh',__LINE__,info)
      ENDIF
      CALL ppm_alloc(nneigh,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','nneigh',__LINE__,info)
      ENDIF
      IF (decomp.NE.ppm_param_decomp_cartesian.AND. &
      &   decomp.NE.ppm_param_decomp_user_defined) THEN
         CALL ppm_alloc(min_box,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','min_box',__LINE__,info)
         ENDIF
         CALL ppm_alloc(max_box,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','max_box',__LINE__,info)
         ENDIF
         CALL ppm_alloc(nchld,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','nchld',__LINE__,info)
         ENDIF
      ENDIF
      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','subs_bc',__LINE__,info)
      ENDIF
      CALL ppm_alloc(istart,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','istart',__LINE__,info)
      ENDIF
      CALL ppm_alloc(ndata,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','ndata',__LINE__,info)
      ENDIF
      CALL ppm_alloc(isublist,ldc,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','isublist',__LINE__,info)
      ENDIF
      IF (.NOT.PRESENT(user_minsub)) THEN
         CALL ppm_alloc(min_sub,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','min_sub',__LINE__,info)
         ENDIF
      ENDIF
      IF (.NOT.PRESENT(user_maxsub)) THEN
         CALL ppm_alloc(max_sub,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','max_sub',__LINE__,info)
         ENDIF
      ENDIF
      IF (.NOT.PRESENT(user_sub2proc)) THEN
         CALL ppm_alloc(sub2proc,ldc,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_dealloc,'ppm_topo_mkfield','sub2proc',__LINE__,info)
         ENDIF
      ENDIF
      9999 CONTINUE
      CALL substop('ppm_topo_mkfield',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mkfield',  &
     &           'Please call ppm_init first!',__LINE__,info)
             GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF(max_phys(i).LE.min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'max_phys must be > min_phys',__LINE__, info)
               GOTO 8888
            ENDIF
            IF(Nm(i) .LT. 2) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'Nm must be > 1 in all dimensions',__LINE__, info)
               GOTO 8888
            ENDIF
            IF(ighostsize(i) .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'ighostsize must be >= 0',__LINE__, info)
               GOTO 8888
            ENDIF
         ENDDO
         IF (Npart .GT. 0 .AND. ASSOCIATED(xp)) THEN
            IF(SIZE(xp,2) .LT. Npart) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'not enough particles contained in xp',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF(SIZE(xp,1) .LT. ppm_dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'leading dimension of xp too small',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF (PRESENT(pcost)) THEN
                IF(SIZE(pcost,1) .LT. Npart) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &                 'pcost does not contain costs for all particles',&
     &                 __LINE__, info)
                   GOTO 8888
                ENDIF
            ENDIF
         ENDIF
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (decomp .NE. ppm_param_decomp_user_defined) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'decomp type is set to user_defined for this assignment',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF(user_nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF (.NOT.ASSOCIATED(user_sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 8888
            ENDIF
            DO i=1,user_nsubs
                IF ((user_sub2proc(i).LT.0).OR.(user_sub2proc(i).GE. &
     &               ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 8888
                ENDIF
            ENDDO
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
            IF(user_nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &             'no subs defined in user_defined decomposition',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF ((.NOT.ASSOCIATED(user_minsub)).OR. &
     &          (.NOT.ASSOCIATED(user_maxsub))) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield', &
     &              'min_sub/max_sub must be allocated for user def. decomp',&
     &              __LINE__, info)
                GOTO 8888
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
            DO i=1,user_nsubs
               larea = (user_maxsub(1,i)-user_minsub(1,i))* &
     &                 (user_maxsub(2,i)-user_minsub(2,i))
               IF(ppm_dim.EQ.3) THEN
                  larea = larea * (user_maxsub(3,i)-user_minsub(3,i))
               END IF
               sarea = sarea + larea
            ENDDO
            IF(ABS(sarea-parea)/parea.GE.lmyeps) THEN
               !----------------------------------------------------------------
               !  Mismatch!
               !----------------------------------------------------------------
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkfield',   &
     &              'faulty subdomains defined',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDIF
 8888    CONTINUE

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
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkfield_d
#endif
