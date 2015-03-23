     !-------------------------------------------------------------------------
     !  Subroutines :               ppm_inl_vlist
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
      SUBROUTINE inl_vlist_s(topoid, xp, Np, Mp, cutoff, skin, lsymm, &
      &          ghostlayer, info, vlist, nvlist, lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE inl_vlist_d(topoid, xp, Np, Mp, cutoff, skin, lsymm, &
      &          ghostlayer, info, vlist, nvlist, lstore)
#endif
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,                           INTENT(IN   ) :: topoid
      !!! ID of the topology.
      REAL(MK), DIMENSION(:,:),          INTENT(IN   ) :: xp
      !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
      INTEGER,                           INTENT(IN   ) :: Np
      !!! Number of real particles
      INTEGER,                           INTENT(IN   ) :: Mp
      !!! Number of all particles including ghost particles
      REAL(MK), DIMENSION(:),            INTENT(IN   ) :: cutoff
      !!! Particle cutoff radii array
      REAL(MK),                          INTENT(IN   ) :: skin
      !!! Skin parameter
      LOGICAL,                           INTENT(IN   ) :: lsymm
      !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
      !!! layers only in (+) directions in all axes. Else, we have ghost
      !!! layers in all directions.
      REAL(MK), DIMENSION(2*ppm_dim),    INTENT(IN   ) :: ghostlayer
      !!! Extra area/volume over the actual domain introduced by
      !!! ghost layers.
      INTEGER,                           INTENT(  OUT) :: info
      !!! Info to be RETURNed. 0 if SUCCESSFUL.
      INTEGER,  DIMENSION(:,:),          POINTER       :: vlist
      !!! verlet lists. vlist(3, 6) is the 3rd neighbor of particle 6.
      INTEGER,  DIMENSION(:),            POINTER       :: nvlist
      !!! number of neighbors of particles. nvlist(i) is number of
      !!! neighbors particle i has.
      LOGICAL,                 OPTIONAL, INTENT(IN   ) :: lstore
      !!! OPTIONAL logical parameter to choose whether to store
      !!! vlist or not. By default, it is set to TRUE.

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo),         POINTER :: topo

      REAL(MK), DIMENSION(2*ppm_dim)    :: actual_subdomain
      REAL(MK), DIMENSION(:,:), POINTER :: xp_sub
      REAL(MK), DIMENSION(:)  , POINTER :: cutoff_sub
      REAL(MK)                          :: t0

      INTEGER,  DIMENSION(:,:), POINTER :: vlist_sub
      INTEGER,  DIMENSION(:)  , POINTER :: nvlist_sub
      INTEGER,  DIMENSION(:)  , POINTER :: p_id
      INTEGER                           :: Np_sub
      INTEGER                           :: Mp_sub
      INTEGER                           :: rank_sub
      INTEGER                           :: neigh_max
      INTEGER                           :: n_part
      INTEGER                           :: isub
      INTEGER                           :: i
      INTEGER                           :: j
      !-------------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !-------------------------------------------------------------------------
      INTEGER                           :: iopt
      INTEGER, DIMENSION(2)             :: lda

      CHARACTER(LEN=ppm_char) :: caller='inl_vlist'

      LOGICAL :: lst

      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check Arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      !---------------------------------------------------------------------
      ! In symmetric verlet lists, we need to allocate nvlist for all
      ! particles including ghost particles (Mp). Otherwise, we need to store
      ! verlet lists of real particles only, hence we allocate nvlist for
      ! real particles (Np) only.
      !---------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      lda(1) = MERGE(Mp,Np,lsymm)
      CALL ppm_alloc(nvlist, lda, iopt, info)
      or_fail_alloc('nvlist',exit_point=no,ppm_error=ppm_error_fatal)

      !---------------------------------------------------------------------
      ! As no neighbors have been found yet, maximum number of neighbors
      ! (neigh_max) is set to 0.
      !---------------------------------------------------------------------
      neigh_max = 0

      !---------------------------------------------------------------------
      ! For each subdomain
      !---------------------------------------------------------------------
      DO rank_sub = 1, topo%nsublist
         isub = topo%isublist(rank_sub)
         !-----------------------------------------------------------------
         ! Get physical extent of subdomain without ghost layers
         !-----------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
         DO i = 1, ppm_dim
            actual_subdomain(2*i-1) = topo%min_subs(i, isub)
            actual_subdomain(2*i)   = topo%max_subs(i, isub)
         ENDDO
#elif __KIND == __DOUBLE_PRECISION
         DO i = 1, ppm_dim
            actual_subdomain(2*i-1) = topo%min_subd(i, isub)
            actual_subdomain(2*i)   = topo%max_subd(i, isub)
         ENDDO
#endif

         !-----------------------------------------------------------------
         ! Get xp and cutoff arrays for the given subdomain. Also get number
         ! of real particles (Np_sub) and total number of particles (Mp_sub)
         ! of this subdomain.
         !-----------------------------------------------------------------
         NULLIFY(xp_sub,cutoff_sub,p_id)
         CALL getSubdomainParticles(xp,Np,Mp,cutoff,lsymm,   &
         &    actual_subdomain,ghostlayer,xp_sub,cutoff_sub, &
         &    Np_sub,Mp_sub,p_id)

         !-----------------------------------------------------------------
         ! Create verlet lists for particles of this subdomain which will
         ! be stored in vlist_sub and nvlist_sub.
         !-----------------------------------------------------------------
         NULLIFY(vlist_sub,nvlist_sub)

         lst = MERGE(lstore,.TRUE.,PRESENT(lstore))

         CALL create_inl_vlist(xp_sub,Np_sub,Mp_sub,cutoff_sub, &
         &    skin,lsymm,actual_subdomain,ghostlayer,info,      &
         &    vlist_sub,nvlist_sub,lst)

         n_part = MERGE(Mp_sub,Np_sub,lsymm)

         DO i = 1, n_part
            nvlist(p_id(i)) = nvlist_sub(i)
         ENDDO

         IF (lst) THEN
            IF (neigh_max.LT.MAXVAL(nvlist_sub))  THEN
               neigh_max = MAXVAL(nvlist_sub)

               iopt = ppm_param_alloc_grow_preserve
               lda(1) = neigh_max
               lda(2) = MERGE(Mp,Np,lsymm)
               CALL ppm_alloc(vlist, lda, iopt, info)
               or_fail_alloc('vlist',exit_point=no,ppm_error=ppm_error_fatal)
            ENDIF

            n_part = MERGE(Mp_sub,Np_sub,lsymm)

            DO i = 1, n_part
               DO j = 1, nvlist_sub(i)
                  vlist(j, p_id(i)) = p_id(vlist_sub(j, i))
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          LOGICAL :: valid
          IF (.NOT. ppm_initialized) THEN
             fail('Please call ppm_init first!',ppm_err_ppm_noinit,exit_point=8888)
          ENDIF
          IF (skin .LT. 0.0_MK) THEN
             fail('skin must be >= 0',exit_point=8888)
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
             fail('Geometric topology required',exit_point=8888)
          ENDIF
          IF (topoid .NE. ppm_param_topo_undefined) THEN
             CALL ppm_check_topoid(topoid,valid,info)
             IF (.NOT. valid) THEN
                fail('topoid out of range',exit_point=8888)
             ENDIF
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE inl_vlist_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE inl_vlist_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE create_inl_vlist_s(xp, Np, Mp, cutoff, skin, lsymm, &
      & actual_domain, ghostlayer, info, vlist, nvlist, lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE create_inl_vlist_d(xp, Np, Mp, cutoff, skin, lsymm, &
      & actual_domain, ghostlayer, info, vlist, nvlist, lstore)
#endif
        !!! This subroutine creates verlet lists for particles whose coordinates
        !!! and cutoff radii are provided by xp and cutoff, respectively.
        !!! Here, Np denotes the number of real particles and Mp is the total
        !!! number of particles including ghost particles. Given these inputs and
        !!! others that are required, this subroutine allocates and fills nvlist.
        !!! If the OPTIONAL parameter lstore is set to TRUE or not passed,
        !!! vlist is also allocated and filled.
#ifdef __DEBUG
        USE ppm_module_util_time
        USE ppm_module_write
#endif
        IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !---------------------------------------------------------------------
        !  Arguments
        !---------------------------------------------------------------------
        REAL(MK), DIMENSION(:,:),       INTENT(IN   ) :: xp
        !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
        INTEGER,                        INTENT(IN   ) :: Np
        !!! Number of real particles
        INTEGER,                        INTENT(IN   ) :: Mp
        !!! Number of all particles including ghost particles
        REAL(MK), DIMENSION(:),         INTENT(IN   ) :: cutoff
        !!! Particles cutoff radii
        REAL(MK),                       INTENT(IN   ) :: skin
        !!! Skin parameter
        LOGICAL,                        INTENT(IN   ) :: lsymm
        !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
        !!! layers only in (+) directions in all axes. Else, we have ghost
        !!! layers in all directions.
        REAL(MK), DIMENSION(2*ppm_dim), INTENT(IN   ) :: actual_domain
        ! Physical extent of actual domain without ghost layers.
        REAL(MK), DIMENSION(ppm_dim),   INTENT(IN   ) :: ghostlayer
        !!! Extra area/volume over the actual domain introduced by
        !!! ghost layers.
        INTEGER,                        INTENT(  OUT) :: info
        !!! Info to be RETURNed. 0 if SUCCESSFUL.
        INTEGER,  DIMENSION(:, :),      POINTER       :: vlist
        !!! verlet lists. vlist(3, 6) is the 3rd neighbor of particle 6.
        INTEGER,  DIMENSION(:),         POINTER       :: nvlist
        !!! number of neighbors of particles. nvlist(i) is number of
        !!! neighbors particle i has.
        LOGICAL,              OPTIONAL, INTENT(IN   ) :: lstore
        !!! OPTIONAL logical parameter to choose whether to store
        !!! vlist or not. By default, it is set to TRUE.

        !---------------------------------------------------------------------
        !  Local variables and counters
        !---------------------------------------------------------------------
        REAL(MK), DIMENSION(2*ppm_dim) :: whole_domain
        ! Physical extent of whole domain including ghost layers.
        REAL(MK)                       :: t0,t1,t2
        REAL(MK)                       :: max_size
        REAL(MK)                       :: size_diff

        INTEGER :: i
        INTEGER :: j

        LOGICAL :: lst

        !---------------------------------------------------------------------
        !  Variables and parameters for ppm_alloc
        !---------------------------------------------------------------------
        INTEGER               :: iopt
        INTEGER, DIMENSION(2) :: lda

        CHARACTER(LEN=ppm_char) :: caller='create_inl_vlist'

        CALL substart(caller,t0,info)

        max_size = 0._MK

        IF (lsymm)  THEN
           DO i = 1, ppm_dim
              whole_domain(2*i-1) = actual_domain(2*i-1)
              whole_domain(2*i)   = actual_domain(2*i) + ghostlayer(i)
              max_size = MAX(max_size, (whole_domain(2*i) - whole_domain(2*i-1)))
           ENDDO
        ELSE
           DO i = 1, ppm_dim
              whole_domain(2*i-1) = actual_domain(2*i-1) - ghostlayer(i)
              whole_domain(2*i)   = actual_domain(2*i)   + ghostlayer(i)
              max_size = MAX(max_size, (whole_domain(2*i) - whole_domain(2*i-1)))
           ENDDO
        ENDIF

        DO i = 1, ppm_dim
           IF ((whole_domain(2*i) - whole_domain(2*i-1)) .LE. max_size/2.0_MK) THEN
              whole_domain(2*i)=whole_domain(2*i-1) + max_size
           ENDIF
        ENDDO

        !-------------------------------------------------------------------------
        !  Create inhomogeneous cell list
        !-------------------------------------------------------------------------
#ifdef __DEBUG
        CALL ppm_util_time(t1)
#endif

        CALL ppm_create_inl_clist(xp, Np, Mp, cutoff, skin, actual_domain, &
        &    ghostlayer, lsymm, clist, info)
        or_fail('ppm_create_inl_clist')

#ifdef __DEBUG
        CALL ppm_util_time(t2)
        stdout("creation of cell list took : ", 't2-t1', " secs")
#endif
        !-------------------------------------------------------------------------
        !  Allocate own_plist array, which will be used to get list of particles
        !  that are located in the same cell. Even though its a temporary array,
        !  it will be used many times throughout the code, hence it is defined
        !  in the module and allocated here, then deallocated in the end when
        !  there is no more need for that.
        !-------------------------------------------------------------------------
        iopt   = ppm_param_alloc_fit
        lda(1) = clist%n_all_p
        CALL ppm_alloc(own_plist, lda, iopt, info)
        or_fail_alloc('own_plist',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Allocate neigh_plist array, which will be used to get list of particles
        !  that are located in neighboring cells. Even though its a temporary
        !  array, it will be used many times throughout the code as own_plist,
        !  hence it is defined in the module and allocated here, then deallocated
        !  in the end when there is no more need for that.
        !-------------------------------------------------------------------------
        CALL ppm_alloc(neigh_plist, lda, iopt, info)
        or_fail_alloc('neigh_plist',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Allocate ncells array, which contains offset directions.
        !-------------------------------------------------------------------------
        lda(1) = 3**ppm_dim
        lda(2) = ppm_dim
        CALL ppm_alloc(ncells, lda, iopt, info)
        or_fail_alloc('ncells',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Fill in ncells array such that it contains offset directions. For
        !  example, ncell(1, :) in 2D will have (-1, -1) which is the bottom-left
        !  neighbor of the reference cell. Works for nD.
        !-------------------------------------------------------------------------
        DO j = 1,ppm_dim
           DO i = 1, 3**(ppm_dim)
              ncells(i,j) = MOD((i-1)/(3**(j-1)),3) - 1
           ENDDO
        ENDDO

        !-------------------------------------------------------------------------
        !  Allocate empty_list array, which will be used to store empty cells.
        !  If not large enough, it will be regrowed in putInEmptyList subroutine.
        !-------------------------------------------------------------------------
        iopt = ppm_param_alloc_fit
        lda(1) = 10*clist%max_depth
        CALL ppm_alloc(empty_list, lda, iopt, info)
        or_fail_alloc('empty_list',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Call getVerletLists subroutine. If lstore is not present, default case
        !  which is lstore = TRUE is applied.
        !-------------------------------------------------------------------------
        lst=MERGE(lstore,.TRUE.,PRESENT(lstore))

        CALL getVerletLists(xp, cutoff, clist, skin, lsymm, whole_domain, &
        &    actual_domain, vlist, nvlist, lst, info)
        or_fail('getVerletLists')

        !-------------------------------------------------------------------------
        !  Deallocate empty_list.
        !-------------------------------------------------------------------------
        iopt = ppm_param_dealloc
        lda = 0
        CALL ppm_alloc(empty_list, lda, iopt, info)
        or_fail_dealloc('empty_list',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Deallocate ncells array.
        !-------------------------------------------------------------------------
        CALL ppm_alloc(ncells,  lda, iopt, info)
        or_fail_dealloc('ncells',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Deallocate own_plist.
        !-------------------------------------------------------------------------
        CALL ppm_alloc(own_plist, lda, iopt, info)
        or_fail_dealloc('own_plist',ppm_error=ppm_error_fatal)

        !-------------------------------------------------------------------------
        !  Deallocate neigh_plist.
        !-------------------------------------------------------------------------
        CALL ppm_alloc(neigh_plist, lda, iopt, info)
        or_fail_dealloc('neigh_plist',ppm_error=ppm_error_fatal)

        CALL ppm_destroy_inl_clist(clist,info)
        or_fail_dealloc('ppm_destroy_inl_clist',ppm_error=ppm_error_fatal)

      9999 CONTINUE
      CALL substop(caller,t0,info)
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE create_inl_vlist_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE create_inl_vlist_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getVerletLists_s(xp, cutoff, clist, skin, lsymm, whole_domain, &
      &                           actual_domain, vlist, nvlist, lstore, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getVerletLists_d(xp, cutoff, clist, skin, lsymm, whole_domain, &
      &                           actual_domain, vlist, nvlist, lstore, info)
#endif
        !!! This subroutine allocates nvlist and fills it with number of
        !!! neighbors of each particle. Then, if lstore is TRUE, it also allocates
        !!! vlist array and fills it with neighbor particles IDs for each
        !!! particle. Depending on lsymm parameter, it calls the appropriate
        !!! subroutine.
        IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !-------------------------------------------------------------------------
        !  Arguments
        !-------------------------------------------------------------------------
        REAL(MK), DIMENSION(:,:),       INTENT(IN   ) :: xp
        !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
        REAL(MK), DIMENSION(:),         INTENT(IN   ) :: cutoff
        !!! Particles cutoff radii
        REAL(MK),                       INTENT(IN   ) :: skin
        !!! Skin parameter
        TYPE(ppm_clist),                INTENT(IN   ) :: clist
        !!! cell list
        LOGICAL,                        INTENT(IN   ) :: lsymm
        !!! Logical parameter to define whether lists are symmetric or not
        !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
        !!! layers only in (+) directions in all axes. Else, we have ghost
        !!! layers in all directions.
        REAL(MK), DIMENSION(2*ppm_dim), INTENT(IN   ) :: whole_domain
        !!! Physical extent of whole domain including ghost layers.
        REAL(MK), DIMENSION(2*ppm_dim), INTENT(IN   ) :: actual_domain
        !!! Physical extent of actual domain without ghost layers.
        INTEGER,  DIMENSION(:, :),      POINTER       :: vlist
        !!! Verlet lists of particles. vlist(j, i) corresponds to jth neighbor
        !!! of particle i.
        INTEGER,  DIMENSION(: ),        POINTER       :: nvlist
        !!! Number of neighbors that particles have. nvlist(i) is the
        !!! number of neighbor particle i has.
        LOGICAL,                        INTENT(IN   ) :: lstore
        !!! Logical parameter to choose whether to store vlist or not.
        INTEGER,                        INTENT(  OUT) :: info

        !-------------------------------------------------------------------------
        !  Local variables and counters
        !-------------------------------------------------------------------------
        REAL(MK) :: t0

        INTEGER               :: p_idx
        !-------------------------------------------------------------------------
        !  Variables and parameters for ppm_alloc
        !-------------------------------------------------------------------------
        INTEGER               :: iopt
        INTEGER, DIMENSION(2) :: lda

        !-------------------------------------------------------------------------
        !  Allocate used array, which will be used as a mask for particles,
        !  to keep track of whether they were used before or not. Then,
        !  initialize it to FALSE.
        !-------------------------------------------------------------------------

        CHARACTER(LEN=ppm_char) :: caller='getVerletLists'

        CALL substart(caller,t0,info)

        iopt = ppm_param_alloc_fit
        lda(1) = clist%n_all_p
        CALL ppm_alloc(used, lda, iopt, info)
        or_fail_alloc('used',ppm_error=ppm_error_fatal)

        used = .FALSE.

        !-------------------------------------------------------------------------
        !  Set size of nvlist. If lsymm = TRUE, we also need to store number of
        !  neighbors of some ghost particles.
        !-------------------------------------------------------------------------
        lda(1) = MERGE(clist%n_all_p,clist%n_real_p,lsymm)

        !-------------------------------------------------------------------------
        !  Allocate nvlist array and initialize it to 0.
        !-------------------------------------------------------------------------
        CALL ppm_alloc(nvlist, lda, iopt, info)
        or_fail_alloc('nvlist',ppm_error=ppm_error_fatal)

        nvlist = 0

        !-------------------------------------------------------------------------
        !  Fill nvlist array with number of neighbors.
        !-------------------------------------------------------------------------
        IF (lsymm) THEN ! If lists are symmetric
           DO p_idx = 1, clist%n_all_p
              CALL count_neigh_sym(clist%rank(p_idx), clist, whole_domain, &
              &    actual_domain, xp, cutoff, skin, nvlist)
           ENDDO
        ELSE ! If lists are not symmetric

           ! Yaser: I think this has a mistake, as we are in the non symmetric
           ! and the SIZE(nvlist) = clist%n_real_p, so I corrected it from clist%n_all_p
           !   DO p_idx = 1, clist%n_all_p
           DO p_idx = 1, clist%n_real_p
              CALL count_neigh(clist%rank(p_idx), clist, whole_domain, &
              &    xp, cutoff, skin, nvlist)
           ENDDO
        ENDIF

        !-------------------------------------------------------------------------
        !  If vlist will be stored,
        !-------------------------------------------------------------------------
        IF (lstore) THEN
           !-----------------------------------------------------------------
           !  Get maximum number of neighbors, for allocation of vlist
           !-----------------------------------------------------------------
           max_nneigh = MAXVAL(nvlist)
           !-----------------------------------------------------------------
           !  Initialize nvlist to 0, since it will be used again to
           !  keep track of where to put the next neighbor on vlist array.
           !-----------------------------------------------------------------
           nvlist = 0
           !-----------------------------------------------------------------
           !  Initialize used to FALSE, since particles should be visited
           !  again.
           !-----------------------------------------------------------------
           used = .FALSE.
           !-----------------------------------------------------------------
           !  As maximum number of neighbors can be max_nneigh at most,
           !  number of rows in vlist can be this number.
           !-----------------------------------------------------------------
           lda(1) = max_nneigh
           !-----------------------------------------------------------------
           !  If we have symmetric lists, we also need to store verlet list
           !  of some ghost particles, so we allocate columns of vlist
           !  by number of all particles.
           !-----------------------------------------------------------------
           lda(2) = MERGE(clist%n_all_p,clist%n_real_p,lsymm)
           !-----------------------------------------------------------------
           !  Allocate vlist
           !-----------------------------------------------------------------
           CALL ppm_alloc(vlist, lda, iopt, info)
           or_fail_alloc('vlist',ppm_error=ppm_error_fatal)

           !-----------------------------------------------------------------
           !  Fill in verlet lists depending on whether lists will be
           !  symmetric or not.
           !-----------------------------------------------------------------
           IF (lsymm) THEN ! If symmetric
              DO p_idx = 1, clist%n_all_p
                 CALL get_neigh_sym(clist%rank(p_idx), clist, whole_domain, &
                 &    actual_domain, xp, cutoff, skin, vlist, nvlist)
              ENDDO
           ELSE            ! If not symmetric
              ! Yaser: I think this has a mistake, as we are in the non symmetric
              ! and the SIZE(nvlist) = clist%n_real_p, so I corrected it from clist%n_all_p
              !   DO p_idx = 1, clist%n_all_p
              DO p_idx = 1, clist%n_real_p
                 CALL get_neigh(clist%rank(p_idx), clist, whole_domain,  &
                 &    xp, cutoff, skin, vlist, nvlist)
              ENDDO
           ENDIF
        ENDIF

        9999 CONTINUE
        CALL substop(caller,t0,info)
        RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getVerletLists_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getVerletLists_d
#endif
