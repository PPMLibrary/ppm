     !-------------------------------------------------------------------------
     !  Subroutines :               ppm_inl_vlist
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


#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE inl_vlist_s(topoid, xp, Np, Mp, cutoff, skin, lsymm,    &
     & ghostlayer, info, vlist, nvlist, lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE inl_vlist_d(topoid, xp, Np, Mp, cutoff, skin, lsymm,    &
     & ghostlayer, info, vlist, nvlist, lstore)
#endif
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,  INTENT(IN)                       :: topoid
      !!! ID of the topology.
      REAL(MK), INTENT(IN), DIMENSION(:,:)       :: xp
      !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
      INTEGER , INTENT(IN)                       :: Np
      !!! Number of real particles
      INTEGER , INTENT(IN)                       :: Mp
      !!! Number of all particles including ghost particles
      REAL(MK), INTENT(IN), DIMENSION(:)         :: cutoff
      !!! Particle cutoff radii array
      REAL(MK), INTENT(IN)                       :: skin
      !!! Skin parameter
      LOGICAL,  INTENT(IN)                       :: lsymm
      !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
      !!! layers only in (+) directions in all axes. Else, we have ghost
      !!! layers in all directions.
      REAL(MK), INTENT(IN), DIMENSION(2*ppm_dim) :: ghostlayer
      !!! Extra area/volume over the actual domain introduced by
      !!! ghost layers.
      INTEGER , INTENT(OUT)                      :: info
      !!! Info to be RETURNed. 0 if SUCCESSFUL.
      INTEGER , POINTER,    DIMENSION(:, :)      :: vlist
      !!! verlet lists. vlist(3, 6) is the 3rd neighbor of particle 6.
      INTEGER , POINTER,    DIMENSION(:)         :: nvlist
      !!! number of neighbors of particles. nvlist(i) is number of
      !!! neighbors particle i has.
      LOGICAL,  INTENT(IN), OPTIONAL             :: lstore
      !!! OPTIONAL logical parameter to choose whether to store
      !!! vlist or not. By default, it is set to TRUE.

      !-------------------------------------------------------------------------
      !  Local variables, arrays and counters
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(2*ppm_dim)             :: actual_subdomain
      REAL(MK), DIMENSION(:,:), POINTER          :: xp_sub     => NULL()
      REAL(MK), DIMENSION(:)  , POINTER          :: cutoff_sub => NULL()
      INTEGER , DIMENSION(:,:), POINTER          :: vlist_sub  => NULL()
      INTEGER , DIMENSION(:)  , POINTER          :: nvlist_sub => NULL()
      INTEGER , DIMENSION(:)  , POINTER          :: p_id       => NULL()
      INTEGER                                    :: Np_sub
      INTEGER                                    :: Mp_sub
      INTEGER                                    :: rank_sub
      INTEGER                                    :: neigh_max
      INTEGER                                    :: n_part
      TYPE(ppm_t_topo)        , POINTER          :: topo      => NULL()
      LOGICAL                                    :: lst
      INTEGER                                    :: i
      INTEGER                                    :: isub
      INTEGER                                    :: j
      REAL(MK)                                   :: t0

      !-------------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !-------------------------------------------------------------------------
      INTEGER                               :: iopt
      INTEGER, DIMENSION(2)                 :: lda

      CALL substart('ppm_inl_vlist',t0,info)

      ! TODO: check wheter topology exists
      topo => ppm_topo(topoid)%t
      !---------------------------------------------------------------------
      ! In symmetric verlet lists, we need to allocate nvlist for all
      ! particles including ghost particles (Mp). Otherwise, we need to store
      ! verlet lists of real particles only, hence we allocate nvlist for
      ! real particles (Np) only.
      !---------------------------------------------------------------------
      IF(lsymm)   THEN
          lda(1) = Mp
      ELSE
          lda(1) = Np
      ENDIF

      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(nvlist, lda, iopt, info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',     &
    &                       'nvlist',__LINE__,info)
      END IF

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
          CALL getSubdomainParticles(xp,Np,Mp,cutoff,lsymm,actual_subdomain,&
 &                 ghostlayer, xp_sub, cutoff_sub, Np_sub,Mp_sub,p_id)

          !-----------------------------------------------------------------
          ! Create verlet lists for particles of this subdomain which will
          ! be stored in vlist_sub and nvlist_sub.
          !-----------------------------------------------------------------
          IF(present(lstore))  THEN
              lst = lstore
          ELSE
              lst = .TRUE.
          END IF

          CALL create_inl_vlist(xp_sub, Np_sub, Mp_sub, cutoff_sub,   &
     &             skin, lsymm, actual_subdomain, ghostlayer, info, vlist_sub,&
     &             nvlist_sub,lst)

          IF(lsymm)  THEN
              n_part = Mp_sub
          ELSE
              n_part = Np_sub
          END IF

          DO i = 1, n_part
              nvlist(p_id(i)) = nvlist_sub(i)
          ENDDO

          IF(lst)  THEN
              IF(neigh_max .LT. MAXVAL(nvlist_sub))  THEN
                  neigh_max = MAXVAL(nvlist_sub)
                  lda(1) = neigh_max
                  IF(lsymm)  THEN
                      lda(2) = Mp
                  ELSE
                      lda(2) = Np
                  ENDIF

                  iopt = ppm_param_alloc_grow_preserve
                  CALL ppm_alloc(vlist, lda, iopt, info)
                  IF (info.NE.0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',        &
                &                       'vlist',__LINE__,info)
                  END IF
              ENDIF

              IF(lsymm)  THEN
                  n_part = Mp_sub
              ELSE
                  n_part = Np_sub
              ENDIF

              DO i = 1, n_part
                  DO j = 1, nvlist_sub(i)
                      vlist(j, p_id(i)) = p_id(vlist_sub(j, i))
                  END DO
              ENDDO
          END IF
      ENDDO
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
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !---------------------------------------------------------------------
      !  Arguments
      !---------------------------------------------------------------------
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: xp
          !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
          INTEGER , INTENT(IN)                       :: Np
          !!! Number of real particles
          INTEGER , INTENT(IN)                       :: Mp
          !!! Number of all particles including ghost particles
          REAL(MK), INTENT(IN), DIMENSION(:)         :: cutoff
          !!! Particles cutoff radii
          REAL(MK), INTENT(IN)                       :: skin
          !!! Skin parameter
          LOGICAL,  INTENT(IN)                       :: lsymm
          !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
          !!! layers only in (+) directions in all axes. Else, we have ghost
          !!! layers in all directions.
          REAL(MK), DIMENSION(2*ppm_dim)             :: actual_domain
          ! Physical extent of actual domain without ghost layers.
          REAL(MK), INTENT(IN), DIMENSION(ppm_dim)   :: ghostlayer
          !!! Extra area/volume over the actual domain introduced by
          !!! ghost layers.
          INTEGER , INTENT(OUT)                      :: info
          !!! Info to be RETURNed. 0 if SUCCESSFUL.
          INTEGER , POINTER,    DIMENSION(:, :)      :: vlist
          !!! verlet lists. vlist(3, 6) is the 3rd neighbor of particle 6.
          INTEGER , POINTER,    DIMENSION(:)         :: nvlist
          !!! number of neighbors of particles. nvlist(i) is number of
          !!! neighbors particle i has.
          LOGICAL,  INTENT(IN), OPTIONAL             :: lstore
          !!! OPTIONAL logical parameter to choose whether to store
          !!! vlist or not. By default, it is set to TRUE.

      !---------------------------------------------------------------------
      !  Local variables and counters
      !---------------------------------------------------------------------
          REAL(MK), DIMENSION(2*ppm_dim)             :: whole_domain
          ! Physical extent of whole domain including ghost layers.
          TYPE(ppm_t_topo), POINTER                  :: topo => NULL()
          ! Topology.
          INTEGER                                    :: i
          INTEGER                                    :: j
          REAL(MK)                                   :: t0
          LOGICAL                                    :: lst

      !---------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !---------------------------------------------------------------------
          INTEGER :: iopt
          INTEGER, DIMENSION(2) :: lda

      !<<<<<<<<<<<<<<<<<<<<<<<<< Start of the code >>>>>>>>>>>>>>>>>>>>>>>>>!
          CALL substart('ppm_neighlist_clist',t0,info)

          IF(lsymm .EQV. .TRUE.)  THEN
              DO i = 1, ppm_dim
                  whole_domain(2*i-1) = actual_domain(2*i-1)
                  whole_domain(2*i)   = actual_domain(2*i) + ghostlayer(i)
              END DO
          ELSE
              DO i = 1, ppm_dim
                  whole_domain(2*i-1) = actual_domain(2*i-1) - ghostlayer(i)
                  whole_domain(2*i)   = actual_domain(2*i)   + ghostlayer(i)
              END DO
          END IF

          n_real_p = Np !Set number of real particles

      !-------------------------------------------------------------------------
      !  Create inhomogeneous cell list
      !-------------------------------------------------------------------------
          CALL ppm_create_inl_clist(xp, Mp, cutoff, skin, actual_domain, &
     & ghostlayer, lsymm, info)
          IF(info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',     &
     &                       'ppm_create_inl_clist',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Allocate own_plist array, which will be used to get list of particles
      !  that are located in the same cell. Even though its a temporary array,
      !  it will be used many times throughout the code, hence it is defined
      !  in the module and allocated here, then deallocated in the end when
      !  there is no more need for that.
      !-------------------------------------------------------------------------
          iopt   = ppm_param_alloc_fit
          lda(1) = n_all_p
          CALL ppm_alloc(own_plist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'own_plist',__LINE__,info)
             GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Allocate neigh_plist array, which will be used to get list of particles
      !  that are located in neighboring cells. Even though its a temporary
      !  array, it will be used many times throughout the code as own_plist,
      !  hence it is defined in the module and allocated here, then deallocated
      !  in the end when there is no more need for that.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(neigh_plist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'neigh_plist',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Allocate ncells array, which contains offset directions.
      !-------------------------------------------------------------------------
          lda(1) = 3**ppm_dim
          lda(2) = ppm_dim
          CALL ppm_alloc(ncells, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'ncells',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Fill in ncells array such that it contains offset directions. For
      !  example, ncell(1, :) in 2D will have (-1, -1) which is the bottom-left
      !  neighbor of the reference cell. Works for nD.
      !-------------------------------------------------------------------------
          DO j = 1,ppm_dim
              DO i = 1, 3**(ppm_dim)
                  ncells(i,j) = MOD((i-1)/(3**(j-1)),3) - 1
              END DO
          END DO

      !-------------------------------------------------------------------------
      !  Allocate empty_list array, which will be used to store empty cells.
      !  If not large enough, it will be regrowed in putInEmptyList subroutine.
      !-------------------------------------------------------------------------
          iopt = ppm_param_alloc_fit
          lda(1) = 10*max_depth
          CALL ppm_alloc(empty_list, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'empty_list',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Call getVerletLists subroutine. If lstore is not present, default case
      !  which is lstore = TRUE is applied.
      !-------------------------------------------------------------------------
          IF(present(lstore)) THEN
              lst = lstore
          ELSE
              lst = .TRUE.
          END IF

          CALL getVerletLists(xp, cutoff, skin, lsymm, whole_domain, &
     &                        actual_domain, vlist, nvlist, lst,info)
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',     &
     &                       'call to getVerletLists failed',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate empty_list.
      !-------------------------------------------------------------------------
          iopt = ppm_param_dealloc
          lda = 0
          CALL ppm_alloc(empty_list, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'empty_list',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate borders array which contains cell lists.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(borders, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'borders',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate rank array.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(rank,    lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'rank',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate ncells array.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(ncells,  lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'ncells',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate own_plist.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(own_plist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'own_plist',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate neigh_plist.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(neigh_plist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'neigh_plist',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Destroy hash table.
      !-------------------------------------------------------------------------
          CALL destroy_htable(info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &            'Could not destroy htable',__LINE__,info)
              GOTO 9999
          END IF

9999  CONTINUE
      CALL substop('create_inl_vlist',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE create_inl_vlist_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE create_inl_vlist_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE getVerletLists_s(xp, cutoff, skin, lsymm, whole_domain,    &
     & actual_domain, vlist, nvlist, lstore,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE getVerletLists_d(xp, cutoff, skin, lsymm, whole_domain,    &
     & actual_domain, vlist, nvlist, lstore,info)
#endif
      !!! This subroutine allocates nvlist and fills it with number of
      !!! neighbors of each particle. Then, if lstore is TRUE, it also allocates
      !!! vlist array and fills it with neighbor particles IDs for each
      !!! particle. Depending on lsymm parameter, it calls the appropriate
      !!! subroutine.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          REAL(MK), INTENT(IN), DIMENSION(:,:)  :: xp
          !!! Particle coordinates array. F.e., xp(1, i) is the x-coor of particle i.
          REAL(MK), INTENT(IN), DIMENSION(:)    :: cutoff
          !!! Particles cutoff radii
          REAL(MK), INTENT(IN)                  :: skin
          !!! Skin parameter
          LOGICAL,  INTENT(IN)                  :: lsymm
          !!! Logical parameter to define whether lists are symmetric or not
          !!! If lsymm = TRUE, verlet lists are symmetric and we have ghost
          !!! layers only in (+) directions in all axes. Else, we have ghost
          !!! layers in all directions.
          REAL(MK), DIMENSION(2*ppm_dim)        :: whole_domain
          !!! Physical extent of whole domain including ghost layers.
          REAL(MK), DIMENSION(2*ppm_dim)        :: actual_domain
          !!! Physical extent of actual domain without ghost layers.
          INTEGER , POINTER,    DIMENSION(:, :) :: vlist
          !!! Verlet lists of particles. vlist(j, i) corresponds to jth neighbor
          !!! of particle i.
          INTEGER , POINTER,    DIMENSION(: )   :: nvlist
          !!! Number of neighbors that particles have. nvlist(i) is the
          !!! number of neighbor particle i has.
          LOGICAL,  INTENT(IN)                  :: lstore
          !!! Logical parameter to choose whether to store vlist or not.

      !-------------------------------------------------------------------------
      !  Local variables and counters
      !-------------------------------------------------------------------------
          INTEGER                               :: i
          INTEGER                               :: p_idx

          REAL(MK)                              :: t0

      !-------------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !-------------------------------------------------------------------------
          INTEGER                               :: iopt
          INTEGER                               :: info
          INTEGER, DIMENSION(2)                 :: lda

      !-------------------------------------------------------------------------
      !  Allocate used array, which will be used as a mask for particles,
      !  to keep track of whether they were used before or not. Then,
      !  initialize it to FALSE.
      !-------------------------------------------------------------------------

      CALL substart('getVerletLists',t0,info)
          iopt = ppm_param_alloc_fit
          lda(1) = n_all_p
          CALL ppm_alloc(used, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'getVerletLists',     &
        &                       'used',__LINE__,info)
              GOTO 9999
          END IF
          used = .FALSE.

      !-------------------------------------------------------------------------
      !  Set size of nvlist. If lsymm = TRUE, we also need to store number of
      !  neighbors of some ghost particles.
      !-------------------------------------------------------------------------
          IF(lsymm .EQV. .TRUE.)  THEN
              lda(1) = n_all_p  ! Store number of neighbors also of ghost particles
          ELSE
              lda(1) = n_real_p ! Store number of neighbors of real particles only
          ENDIF

      !-------------------------------------------------------------------------
      !  Allocate nvlist array and initialize it to 0.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(nvlist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'getVerletLists',     &
        &                       'nvlist',__LINE__,info)
              GOTO 9999
          END IF
          nvlist = 0

      !-------------------------------------------------------------------------
      !  Fill nvlist array with number of neighbors.
      !-------------------------------------------------------------------------
          IF(lsymm .EQV. .TRUE.)   THEN ! If lists are symmetric
              DO i = 1, n_all_p
                  p_idx = rank(i)
                  CALL count_neigh_sym(p_idx, whole_domain,       &
                  & actual_domain, xp, cutoff, skin, vlist, nvlist)
              END DO
          ELSE                          ! If lists are not symmetric
              DO i = 1, n_all_p
                  p_idx = rank(i)
                  CALL count_neigh(p_idx, whole_domain, xp, cutoff,    &
                  & skin, vlist, nvlist)
              END DO
          END IF

      !-------------------------------------------------------------------------
      !  If vlist will be stored,
      !-------------------------------------------------------------------------
          IF(lstore) THEN
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
              IF(lsymm .EQV. .TRUE.)  THEN
                  lda(2) = n_all_p
              ELSE
                  lda(2) = n_real_p
              ENDIF
              !-----------------------------------------------------------------
              !  Allocate vlist
              !-----------------------------------------------------------------
              CALL ppm_alloc(vlist, lda, iopt, info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'getVerletLists',     &
            &                       'vlist',__LINE__,info)
                  GOTO 9999
              END IF

              !-----------------------------------------------------------------
              !  Fill in verlet lists depending on whether lists will be
              !  symmetric or not.
              !-----------------------------------------------------------------
              IF(lsymm .EQV. .TRUE.)    THEN ! If symmetric
                  DO p_idx = 1, n_all_p
                      CALL get_neigh_sym(rank(p_idx),       &
                      & whole_domain, actual_domain, xp, cutoff, skin, vlist, nvlist)
                  END DO
              ELSE                           ! If not symmetric
                  DO p_idx = 1, n_all_p
                      CALL get_neigh(rank(p_idx), whole_domain,  &
                      & xp, cutoff, skin, vlist, nvlist)
                  END DO
              END IF
          END IF

9999      CONTINUE
          CALL substop('getVerletLists',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE getVerletLists_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE getVerletLists_d
#endif