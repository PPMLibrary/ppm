      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_mkgeom
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
      SUBROUTINE ppm_topo_mkgeom_s(topoid,decomp,assig,min_phys,         &
     &              max_phys,bcdef,ghostsize,cost,info,user_minsub,      &
     &              user_maxsub,user_nsubs,user_sub2proc)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_topo_mkgeom_d(topoid,decomp,assig,min_phys,         &
     &              max_phys,bcdef,ghostsize,cost,info,user_minsub,      &
     &              user_maxsub,user_nsubs,user_sub2proc)
#endif
      !!! This routine is the topology creation routine for purely
      !!! geometry-based decompositions, i.e. without particles and without
      !!! meshes.
      !!!
      !!! [NOTE]
      !!! In the user_defined case, the user must submit existing
      !!! subdomains and all of `min_sub`, `max_sub`, `cost` and `user_nsubs`
      !!! must be provided. Subs are then mapped onto processors and
      !!! the topology is stored.

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_topo_cost
      USE ppm_module_topo_store
      USE ppm_module_define_subs_bc
      USE ppm_module_topo_subs2proc
      USE ppm_module_topo_metis_s2p
      USE ppm_module_find_neigh
      USE ppm_module_tree
      USE ppm_module_alloc
      USE ppm_module_topo_box2subs
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
      !!! * ppm_param_decomp_cuboid
      !!! * ppm_param_decomp_user_defined
      INTEGER,                            INTENT(IN   ) :: assig
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
      REAL(MK), DIMENSION(:  ),           INTENT(IN   ) :: min_phys
      !!! Minimum of physical extend of the computational domain (double)
      !!!
      !!! first index is ppm_dim
      REAL(MK), DIMENSION(:  ),           INTENT(IN   ) :: max_phys
      !!! Maximum of physical extend of the computational domain (double)
      !!!
      !!! first index is ppm_dim
      INTEGER , DIMENSION(:  ),           INTENT(IN   ) :: bcdef
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
      !!! The size (width) of the ghost layer.
      REAL(MK), DIMENSION(:  ),           POINTER       :: cost
      !!! Estimated cost associated with subdomains. Either user-defined on
      !!! input or decomposition result on output. The cost of a subdomain
      !!! is given by its volume.
      INTEGER,                            INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
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
      !!! Used when decomp is user defined.
      INTEGER,  DIMENSION(:  ), OPTIONAL, POINTER       :: user_sub2proc
      !!! Subdomain to processor assignment.
      !!! Used if assignment is user defined.
      !!!
      !!! index: subID (global).
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                           :: i,treetype,nbox,isub
      INTEGER                           :: iopt
      INTEGER, DIMENSION(:,:), POINTER  :: ineigh  => NULL()
      INTEGER, DIMENSION(:,:), POINTER  :: subs_bc => NULL()
      INTEGER, DIMENSION(1  )           :: Nmdummy,ldc
      INTEGER, DIMENSION(3,1)           :: nnodes
      INTEGER, DIMENSION(:  ), POINTER  :: nneigh  => NULL()
      INTEGER, DIMENSION(:  ), POINTER  :: nchld   => NULL()
      REAL(MK)                          :: t0,parea,sarea,larea,lmyeps
      REAL(MK), DIMENSION(ppm_dim)      :: gsvec
      REAL(MK), DIMENSION(1,1)          :: xpdummy
      LOGICAL , DIMENSION(ppm_dim)      :: fixed
      REAL(MK), DIMENSION(3,2)          :: weights
      REAL(MK), DIMENSION(:,:), POINTER :: min_box  => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_box  => NULL()
      CHARACTER(LEN=ppm_char)           :: mesg
      INTEGER                           :: nsublist, nsubs
      INTEGER , DIMENSION(  :), POINTER :: isublist => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: min_sub  => NULL()
      REAL(MK), DIMENSION(:,:), POINTER :: max_sub  => NULL()
      INTEGER,  DIMENSION(:  ), POINTER :: sub2proc => NULL()
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_mkgeom',t0,info)
#if    __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif  __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

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
      ENDIF
      IF (PRESENT(user_maxsub)) THEN
          max_sub => user_maxsub
      ENDIF
      IF (PRESENT(user_sub2proc)) THEN
          sub2proc => user_sub2proc
      ENDIF
      !-------------------------------------------------------------------------
      !  Dummy arguments for non-existing particles and meshes
      !-------------------------------------------------------------------------
      xpdummy(1,1)  = 0.0_MK
      Nmdummy(1)    = 0
      nnodes(1:3,1) = 0

      !-------------------------------------------------------------------------
      !  Recursive bisection
      !-------------------------------------------------------------------------
      SELECT CASE (decomp)
      CASE (ppm_param_decomp_bisection)
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
         &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
         &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
             &    'Bisection decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Pencils
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_xpencil, &
      &     ppm_param_decomp_ypencil, &
      &     ppm_param_decomp_zpencil)
         IF (decomp.EQ.ppm_param_decomp_zpencil.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
             &   'Cannot make z pencils in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  pencil quadrisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a quad tree, binary in 2d
         treetype         = ppm_param_tree_quad
         IF (ppm_dim .EQ. 2) treetype = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! fix the proper direction
         fixed(1:ppm_dim) = .FALSE.
         IF (decomp .EQ. ppm_param_decomp_xpencil) fixed(1) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_ypencil) fixed(2) = .TRUE.
         IF (decomp .EQ. ppm_param_decomp_zpencil) fixed(3) = .TRUE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
         &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
         &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
             &    'Pencil decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &       max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  Slabs
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_xy_slab, &
      &     ppm_param_decomp_xz_slab, &
      &     ppm_param_decomp_yz_slab)
         IF (decomp.EQ.ppm_param_decomp_xz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
             &   'Cannot make x-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         IF (decomp.EQ.ppm_param_decomp_yz_slab.AND.ppm_dim.LT.3) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',  &
             &           'Cannot make y-z slabs in 2D!',__LINE__,info)
             GOTO 9999
         ENDIF
         !-------------------------------------------------------------------
         !  slab bisection using the general ppm_tree
         !-------------------------------------------------------------------
         ! build a binary tree
         treetype         = ppm_param_tree_bin
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
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
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
         &    ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
         &    max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',    &
             &    'Slab decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
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
         ! no particles and no mesh
         weights(1,1:2)   = 0.0_MK
         weights(2,1:2)   = 0.0_MK
         ! geometry has unit weight
         weights(3,1:2)   = 1.0_MK
         ! all directions can be cut
         fixed(1:ppm_dim) = .FALSE.
         gsvec(1:ppm_dim) = ghostsize
         ! build tree
         CALL ppm_tree(xpdummy,0,Nmdummy,min_phys,max_phys,treetype,  &
         &       ppm_nproc,.FALSE.,gsvec,0.1_MK,-1.0_MK,fixed,weights,min_box, &
         &       max_box,nbox,nchld,info)
         IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
             &   'Cuboid decomposition failed',__LINE__,info)
             GOTO 9999
         ENDIF
         ! convert tree to subs
         CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,   &
         &    max_sub,nsubs,info)
         IF (info .NE. 0) GOTO 9999

      !-------------------------------------------------------------------------
      !  User provides decomposition: Do nothing
      !-------------------------------------------------------------------------
      CASE (ppm_param_decomp_user_defined)
         !Do nothing. Just take the stuff from the user and trust the guy.
         gsvec(1:ppm_dim) = ghostsize
      !-------------------------------------------------------------------------
      !  Unknown decomposition type
      !-------------------------------------------------------------------------
      CASE DEFAULT
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown decomposition type: ',decomp
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
         &    mesg,__LINE__,info)
         GOTO 9999
      END SELECT

      !-------------------------------------------------------------------------
      !  Find the neighbors of the subdomains
      !-------------------------------------------------------------------------
      CALL ppm_find_neigh(min_phys,max_phys,bcdef, &
      &    min_sub,max_sub,nsubs,nneigh,ineigh,gsvec,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
          &    'Finding neighbors failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Find the cost of each subdomain
      !-------------------------------------------------------------------------
      IF (decomp .NE. ppm_param_decomp_user_defined) THEN
          CALL ppm_topo_cost(xpdummy,0,min_sub,max_sub,nsubs,nnodes,  &
          &    cost,info)
          IF (info.NE.0) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
             &   'Computing costs failed',__LINE__,info)
             GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign the subdomains to processors
      !-------------------------------------------------------------------------
      SELECT CASE (assig)
      CASE (ppm_param_assign_internal)
         !-------------------------------------------------------------------
         !  internal assignment routine
         !-------------------------------------------------------------------
         CALL ppm_topo_subs2proc(cost,nneigh,ineigh,nsubs,sub2proc, &
         &    isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
            &   'Assigning subs to processors failed',__LINE__,info)
            GOTO 9999
         ENDIF
      CASE (ppm_param_assign_nodal_cut,  &
      &     ppm_param_assign_nodal_comm, &
      &     ppm_param_assign_dual_cut,   &
      &     ppm_param_assign_dual_comm)
         !-------------------------------------------------------------------
         !  use METIS library to do assignment
         !-------------------------------------------------------------------
         CALL ppm_topo_metis_s2p(min_sub,max_sub,nneigh,ineigh,cost,nsubs,&
         &       assig,sub2proc,isublist,nsublist,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
            &    'Assigning subs to processors using METIS failed',__LINE__,&
            &    info)
            GOTO 9999
         ENDIF
      CASE (ppm_param_assign_user_defined)
         !-------------------------------------------------------------------
         !  user defined assignment
         !-------------------------------------------------------------------
         iopt = ppm_param_alloc_fit
         ldc(1) = nsubs
         CALL ppm_alloc(isublist,ldc,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_topo_mkgeom',   &
     &           'list of local subs ISUBLIST',__LINE__,info)
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
      CASE DEFAULT
         !-------------------------------------------------------------------
         !  unknown assignment scheme
         !-------------------------------------------------------------------
         info = ppm_error_error
         WRITE(mesg,'(A,I5)') 'Unknown assignment scheme: ',assig
         CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
     &       mesg,__LINE__,info)
         GOTO 9999
      END SELECT

      !-------------------------------------------------------------------------
      !  Find and define the boundary conditions on the subs on the local
      !  processor (the routine will allocate the requried memory)
      !-------------------------------------------------------------------------
      NULLIFY(subs_bc)
      CALL ppm_define_subs_bc(min_phys,max_phys,bcdef,min_sub,max_sub, &
      &    nsubs,subs_bc,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkpart',  &
         &    'finding and defining the BC of the subs failed ', &
         &    __LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the topology internally
      !-------------------------------------------------------------------------
      CALL ppm_topo_store(topoid,min_phys,max_phys,min_sub,max_sub,subs_bc, &
     &                    sub2proc,nsubs,bcdef,ghostsize,isublist,nsublist,&
     &                    nneigh,ineigh,info)
      IF (info.NE.0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_topo_mkgeom',  &
     &        'Storing topology failed',__LINE__,info)
          GOTO 9999
      ENDIF

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

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_mkgeom',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (.NOT. ppm_initialized) THEN
             info = ppm_error_error
             CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_mkgeom',  &
     &           'Please call ppm_init first!',__LINE__,info)
             GOTO 8888
         ENDIF
         IF(ghostsize .LT. 0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &          'ghostsize must be >= 0.0',__LINE__, info)
            GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF(max_phys(i).LE.min_phys(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'max_phys must be > min_phys',__LINE__, info)
               GOTO 8888
            ENDIF
         ENDDO
         IF (assig .EQ. ppm_param_assign_user_defined) THEN
            IF (decomp .NE. ppm_param_decomp_user_defined) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'decomp type is set to user_defined for this assignment',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF(user_nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'no subs defined in user_defined assignment',&
     &             __LINE__, info)
               GOTO 8888
            ENDIF
            IF (.NOT.ASSOCIATED(sub2proc)) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &              'sub2proc must be allocated for user defined assignment',&
     &              __LINE__, info)
                GOTO 8888
            ENDIF
            DO i=1,user_nsubs
                IF ((sub2proc(i).LT.0).OR.(sub2proc(i).GE.ppm_nproc)) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &                 'invalid processor specified in sub2proc',&
     &                 __LINE__, info)
                   GOTO 8888
                ENDIF
            ENDDO
         ENDIF
         IF (decomp .EQ. ppm_param_decomp_user_defined) THEN
            IF ((.NOT.ASSOCIATED(min_sub)).OR.(.NOT.ASSOCIATED(max_sub))) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &              'min_sub/max_sub must be allocated for user def. decomp',&
     &              __LINE__, info)
                GOTO 8888
            ENDIF
            IF(user_nsubs .LE. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom', &
     &             'no subs defined in user_defined decomposition',&
     &             __LINE__, info)
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
               larea = (max_sub(1,i)-min_sub(1,i))*(max_sub(2,i)-min_sub(2,i))
               IF(ppm_dim.EQ.3) THEN
                  larea = larea * (max_sub(3,i)-min_sub(3,i))
               END IF
               sarea = sarea + larea
            ENDDO
            IF(ABS(sarea-parea)/parea.GE.lmyeps) THEN
               !----------------------------------------------------------------
               !  Mismatch!
               !----------------------------------------------------------------
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_topo_mkgeom',   &
     &              'faulty subdomains defined',__LINE__,info)
               GOTO 8888
            ENDIF

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
 8888    CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_topo_mkgeom_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_topo_mkgeom_d
#endif
