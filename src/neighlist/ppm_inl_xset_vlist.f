     !-------------------------------------------------------------------------
     !  Subroutines :               ppm_inl_xset_vlist
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
#if   __MODE == __INL
      SUBROUTINE inl_xset_vlist_s(topoid,red,nred,mred, &
      &          blue,nblue,mblue,rcblue,skin,          &
      &          ghostlayer,info,vlist,nvlist,lstore)
#elif __MODE == __HNL
      SUBROUTINE hnl_xset_vlist_s(topoid,red,nred,mred, &
      &          blue,nblue,mblue,cutoff,skin,          &
      &          ghostlayer,info,vlist,nvlist,lstore)
#endif
#elif __KIND == __DOUBLE_PRECISION
#if   __MODE == __INL
      SUBROUTINE inl_xset_vlist_d(topoid,red,nred,mred, &
      &          blue,nblue,mblue,rcblue,skin,          &
      &          ghostlayer,info,vlist,nvlist,lstore)
#elif __MODE == __HNL
      SUBROUTINE hnl_xset_vlist_d(topoid,red,nred,mred, &
      &          blue,nblue,mblue,cutoff,skin,          &
      &          ghostlayer,info,vlist,nvlist,lstore)
#endif
#endif
      !!! Inhomogeneous cross-set neighborlists. This routine provides
      !!! Verlet-list-like neighbor lists of particles of one set (blue) to a
      !!! set of other particles (red). The neighborlists are built analogous to
      !!! the inhomogenous neighborlists described in Awile2011. The red
      !!! particles are assumed to have interaction radius = domain size and we
      !!! use the same neighborhood relation as in INL.

      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER,  INTENT(IN)                       :: topoid
      !!! ID of the topology.
      REAL(MK), INTENT(IN), DIMENSION(:,:)       :: red
      !!! Coordinates array for reference particles (red)
      INTEGER , INTENT(IN)                       :: nred
      !!! Number of real red particles
      INTEGER , INTENT(IN)                       :: mred
      !!! Number of all red particles including ghost particles
      REAL(MK), INTENT(IN), DIMENSION(:,:)       :: blue
      !!! Coordinates array for neighbor particle set (blue)
      INTEGER , INTENT(IN)                       :: nblue
      !!! Number of real blue particles
      INTEGER , INTENT(IN)                       :: mblue
      !!! Number of all blue particles including ghost particles
#if   __MODE == __INL
      REAL(MK), INTENT(IN), DIMENSION(:)         :: rcblue
      !!! Blue particle cutoff radii array
#elif __MODE == __HNL
      REAL(MK), INTENT(IN)                       :: cutoff
      !!! Blue particle cutoff radii scalar
#endif
      REAL(MK), INTENT(IN)                       :: skin
      !!! Skin parameter
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
      TYPE(ppm_t_topo),         POINTER     :: topo

      REAL(MK), DIMENSION(2*ppm_dim)        :: curr_sub
      REAL(MK), DIMENSION(:,:), POINTER     :: red_sub    => NULL()
      REAL(MK), DIMENSION(:,:), POINTER     :: blue_sub   => NULL()
      REAL(MK), DIMENSION(:)  , POINTER     :: rcred_sub  => NULL()
      REAL(MK), DIMENSION(:)  , POINTER     :: rcblue_sub => NULL()
      REAL(MK)                              :: t0
      REAL(MK), DIMENSION(:),   POINTER     :: rcred      => NULL()
      REAL(MK)                              :: max_sub_size
#if   __MODE == __HNL
      REAL(MK), DIMENSION(:), POINTER, SAVE :: rcblue
#endif
      INTEGER , DIMENSION(:,:), POINTER     :: vlist_sub  => NULL()
      INTEGER , DIMENSION(:)  , POINTER     :: nvlist_sub => NULL()
      INTEGER , DIMENSION(:)  , POINTER     :: red_p_id   => NULL()
      INTEGER , DIMENSION(:)  , POINTER     :: blue_p_id  => NULL()
      INTEGER                               :: nred_sub
      INTEGER                               :: mred_sub
      INTEGER                               :: nblue_sub
      INTEGER                               :: mblue_sub
      INTEGER                               :: rank_sub
      INTEGER                               :: neigh_max
      INTEGER                               :: n_part
      INTEGER                               :: i
      INTEGER                               :: isub
      INTEGER                               :: j

      LOGICAL                               :: lst

      !-------------------------------------------------------------------------
      !  Variables and parameters for ppm_alloc
      !-------------------------------------------------------------------------
      INTEGER                               :: iopt
      INTEGER, DIMENSION(2)                 :: lda

      CALL substart('ppm_inl_vlist',t0,info)

#if   __MODE == __HNL
      !create an array of homogeneous cutoffs
      ! TODO: this is an inefficient hack
      ! (one would need to use the routines for homogeneous lists instead)
      lda(1) = mblue
      iopt = ppm_param_alloc_grow
      CALL ppm_alloc(rcblue,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',     &
    &                       'rcblue',__LINE__,info)
      END IF
      rcblue = cutoff
#endif
      ! TODO: check wheter topology exists
      topo => ppm_topo(topoid)%t

      lda(1) = nred
      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(nvlist, lda, iopt, info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',     &
    &                       'nvlist',__LINE__,info)
      END IF
      lda(1) = mred
      CALL ppm_alloc(rcred, lda, iopt, info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',     &
    &                       'rcred',__LINE__,info)
      END IF
      !---------------------------------------------------------------------
      ! set rcred to a large value to make sure all red particles are on the
      ! top level of the cell tree
      !---------------------------------------------------------------------
      max_sub_size=0._MK
      DO rank_sub = 1, topo%nsublist
          isub = topo%isublist(rank_sub)
#if   __KIND == __SINGLE_PRECISION
          max_sub_size = MAX(max_sub_size,&
              MAXVAL(topo%max_subs(1:ppm_dim, isub) - &
              topo%min_subs(1:ppm_dim, isub)))
#elif   __KIND == __DOUBLE_PRECISION
          max_sub_size = MAX(max_sub_size,&
              MAXVAL(topo%max_subd(1:ppm_dim, isub) - &
              topo%min_subd(1:ppm_dim, isub)))
#endif
      ENDDO

      rcred(1:mred) = max_sub_size
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
              curr_sub(2*i-1) = topo%min_subs(i, isub)
              curr_sub(2*i)   = topo%max_subs(i, isub)
          ENDDO
#elif __KIND == __DOUBLE_PRECISION
          DO i = 1, ppm_dim
              curr_sub(2*i-1) = topo%min_subd(i, isub)
              curr_sub(2*i)   = topo%max_subd(i, isub)
          ENDDO
#endif

          !-----------------------------------------------------------------
          ! Get xp and cutoff arrays for the given subdomain. Also get number
          ! of real particles (Np_sub) and total number of particles (Mp_sub)
          ! of this subdomain.
          !-----------------------------------------------------------------
          CALL getSubdomainParticles(red,nred,mred,rcred,.FALSE.,  &
          &                 curr_sub,ghostlayer,red_sub,rcred_sub, &
          &                 nred_sub,mred_sub,red_p_id)

          CALL getSubdomainParticles(blue,nblue,mblue,rcblue,.FALSE.,&
          &                 curr_sub,ghostlayer,blue_sub,rcblue_sub, &
          &                 nblue_sub,mblue_sub,blue_p_id)

          !-----------------------------------------------------------------
          ! Create verlet lists for particles of this subdomain which will
          ! be stored in vlist_sub and nvlist_sub.
          !-----------------------------------------------------------------
          lst=MERGE(lstore,.TRUE.,PRESENT(lstore))

          CALL create_inl_xset_vlist(red_sub,nred_sub,mred_sub,rcred_sub, &
          &                 blue_sub,nblue_sub,mblue_sub,rcblue_sub,skin, &
          &                 curr_sub,ghostlayer,info,vlist_sub,nvlist_sub,lst)

          n_part = nred_sub
          DO i = 1, n_part
              nvlist(red_p_id(i)) = nvlist_sub(i)
          ENDDO

          IF(lst)  THEN
              IF(neigh_max .LT. MAXVAL(nvlist_sub))  THEN
                  neigh_max = MAXVAL(nvlist_sub)
                  lda(1) = neigh_max
                  lda(2) = nred

                  iopt = ppm_param_alloc_grow_preserve
                  CALL ppm_alloc(vlist, lda, iopt, info)
                  IF (info.NE.0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_inl_vlist',        &
                &                       'vlist',__LINE__,info)
                  END IF
              ENDIF

              n_part = nred_sub

              DO i = 1, n_part
                  DO j = 1, nvlist_sub(i)
                      vlist(j, red_p_id(i)) = blue_p_id(vlist_sub(j, i))
                  END DO
              ENDDO
          END IF
      ENDDO
      CALL substop('ppm_inl_vlist',t0,info)
#if   __KIND == __SINGLE_PRECISION
#if   __MODE == __INL
      END SUBROUTINE inl_xset_vlist_s
#elif __MODE == __HNL
      END SUBROUTINE hnl_xset_vlist_s
#endif
#elif __KIND == __DOUBLE_PRECISION
#if   __MODE == __INL
      END SUBROUTINE inl_xset_vlist_d
#elif __MODE == __HNL
      END SUBROUTINE hnl_xset_vlist_d
#endif
#endif

      !if __MODE == _HNL -> skip till end of file
#if   __MODE == __INL
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE create_inl_xset_vlist_s(red,nred,mred,rcred,blue,nblue,  &
      &                        mblue,rcblue,skin,curr_dom,ghostlayer,info,&
      &                        vlist, nvlist, lstore)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE create_inl_xset_vlist_d(red,nred,mred,rcred,blue,nblue,  &
      &                        mblue,rcblue,skin,curr_dom,ghostlayer,info,&
      &                        vlist, nvlist, lstore)
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
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: red
          !!! Coordinates array for refernece particles (red)
          INTEGER , INTENT(IN)                       :: nred
          !!! Number of real red particles
          INTEGER , INTENT(IN)                       :: mred
          !!! Number of all red particles
          REAL(MK), INTENT(IN), DIMENSION(:)         :: rcred
          !!! Red particles cutoff radii
          REAL(MK), INTENT(IN), DIMENSION(:,:)       :: blue
          !!! Coordinates array for refernece particles (blue)
          INTEGER , INTENT(IN)                       :: nblue
          !!! Number of real blue particles
          INTEGER , INTENT(IN)                       :: mblue
          !!! Number of all blue particles
          REAL(MK), INTENT(IN), DIMENSION(:)         :: rcblue
          !!! blue particles cutoff radii
          REAL(MK), INTENT(IN)                       :: skin
          !!! Skin parameter
          REAL(MK), DIMENSION(2*ppm_dim)             :: curr_dom
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

          DO i = 1, ppm_dim
              whole_domain(2*i-1) = curr_dom(2*i-1) - ghostlayer(i)
              whole_domain(2*i)   = curr_dom(2*i)   + ghostlayer(i)
          END DO

      !-------------------------------------------------------------------------
      !  Create inhomogeneous cell list for red and for blue
      !-------------------------------------------------------------------------
          CALL ppm_create_inl_clist(red, nred, mred, rcred, skin, curr_dom, &
          &                         ghostlayer, .FALSE., red_clist, info)
          IF(info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',     &
     &                    'ppm_create_inl_clist with red',__LINE__,info)
              GOTO 9999
          END IF

          CALL ppm_create_inl_clist(blue, nblue, mblue, rcblue, skin, curr_dom, &
     & ghostlayer, .FALSE., blue_clist, info)
          IF(info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',     &
     &                    'ppm_create_inl_clist with blue',__LINE__,info)
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
          lda(1) = red_clist%n_all_p
          CALL ppm_alloc(own_red, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'own_red',__LINE__,info)
             GOTO 9999
          END IF
          CALL ppm_alloc(own_blue, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'own_blue',__LINE__,info)
             GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Allocate neigh_plist array, which will be used to get list of particles
      !  that are located in neighboring cells. Even though its a temporary
      !  array, it will be used many times throughout the code as own_plist,
      !  hence it is defined in the module and allocated here, then deallocated
      !  in the end when there is no more need for that.
      !-------------------------------------------------------------------------
          lda(1) = blue_clist%n_all_p
          CALL ppm_alloc(neigh_red, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'neigh_red',__LINE__,info)
              GOTO 9999
          END IF
          CALL ppm_alloc(neigh_blue, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'neigh_blue',__LINE__,info)
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
          lda(1) = 10*blue_clist%max_depth
          CALL ppm_alloc(empty_list, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'create_inl_vlist',     &
     &                       'empty_list',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Call get_xset_VerletLists subroutine. If lstore is not present, default case
      !  which is lstore = TRUE is applied.
      !-------------------------------------------------------------------------
          IF(PRESENT(lstore)) THEN
              lst = lstore
          ELSE
              lst = .TRUE.
          END IF
          CALL get_xset_VerletLists(red, rcred, red_clist, blue, rcblue, blue_clist, &
 &                            skin, whole_domain, &
     &                        curr_dom, vlist, nvlist, lst,info)
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',     &
     &                       'call to get_xset_VerletLists failed',__LINE__,info)
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
          CALL ppm_alloc(own_red, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'own_red',__LINE__,info)
              GOTO 9999
          END IF
          CALL ppm_alloc(own_blue, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'own_blue',__LINE__,info)
              GOTO 9999
          END IF

      !-------------------------------------------------------------------------
      !  Deallocate neigh_plist.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(neigh_red, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'neigh_red',__LINE__,info)
              GOTO 9999
          END IF
          CALL ppm_alloc(neigh_blue, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_dealloc,'create_inl_vlist',   &
     &                       'neigh_blue',__LINE__,info)
              GOTO 9999
          END IF
          CALL ppm_destroy_inl_clist(red_clist,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',   &
     &                       'ppm_destroy_inl_clist red',__LINE__,info)
              GOTO 9999
          END IF
          CALL ppm_destroy_inl_clist(blue_clist,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_sub_failed,'create_inl_vlist',   &
     &                       'ppm_destroy_inl_clist blue',__LINE__,info)
              GOTO 9999
          END IF

9999  CONTINUE
      CALL substop('create_inl_vlist',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE create_inl_xset_vlist_s
#elif   __KIND == __DOUBLE_PRECISION
      END SUBROUTINE create_inl_xset_vlist_d
#endif

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE get_xset_VerletLists_s(red, rcred, red_clist, blue, rcblue, &
 &                                blue_clist, skin, whole_domain,    &
 &                                curr_dom, vlist, nvlist, lstore,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE get_xset_VerletLists_d(red, rcred, red_clist, blue, rcblue, &
 &                                blue_clist, skin, whole_domain,    &
 &                                curr_dom, vlist, nvlist, lstore,info)
#endif
      !!! This subroutine allocates nvlist and fills it with number of
      !!! neighbors of each particle. Then, if lstore is TRUE, it also allocates
      !!! vlist array and fills it with neighbor particles IDs for each
      !!! particle.
          IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
          INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
          REAL(MK), INTENT(IN), DIMENSION(:,:)  :: red
          !!! coordinates array of red particles
          REAL(MK), INTENT(IN), DIMENSION(:)    :: rcred
          !!! Red particles cutoff radii
          TYPE(ppm_clist), INTENT(IN)           :: red_clist
          !!! Red particle cell list
          REAL(MK), INTENT(IN), DIMENSION(:,:)  :: blue
          !!! coordinates array of blue particles
          REAL(MK), INTENT(IN), DIMENSION(:)    :: rcblue
          !!! Blue particles cutoff radii
          TYPE(ppm_clist), INTENT(IN)           :: blue_clist
          !!! Blue particle cell list
          REAL(MK), INTENT(IN)                  :: skin
          !!! Skin parameter
          REAL(MK), DIMENSION(2*ppm_dim)        :: whole_domain
          !!! Physical extent of whole domain including ghost layers.
          REAL(MK), DIMENSION(2*ppm_dim)        :: curr_dom
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

      CALL substart('get_xset_VerletLists',t0,info)
          iopt = ppm_param_alloc_fit
          lda(1) = red_clist%n_all_p
          CALL ppm_alloc(used, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'get_xset_VerletLists',     &
        &                       'used',__LINE__,info)
              GOTO 9999
          END IF
          used = .FALSE.

      !-------------------------------------------------------------------------
      !  Set size of nvlist.
      !-------------------------------------------------------------------------
          lda(1) = red_clist%n_real_p ! Store number of neighbors of real particles only
      !-------------------------------------------------------------------------
      !  Allocate nvlist array and initialize it to 0.
      !-------------------------------------------------------------------------
          CALL ppm_alloc(nvlist, lda, iopt, info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'get_xset_VerletLists',     &
        &                       'nvlist',__LINE__,info)
              GOTO 9999
          END IF
          nvlist = 0


      !-------------------------------------------------------------------------
      !  Fill nvlist array with number of neighbors.
      !-------------------------------------------------------------------------
          DO i = 1, red_clist%n_real_p
              p_idx = red_clist%rank(i)
              CALL count_xset_neigh(p_idx, red_clist, blue_clist, &
 &                 whole_domain, &
 &                 red, rcred, blue, rcblue, skin, nvlist)
          END DO

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
              lda(2) = red_clist%n_real_p
              !-----------------------------------------------------------------
              !  Allocate vlist
              !-----------------------------------------------------------------
              CALL ppm_alloc(vlist, lda, iopt, info)
              IF (info.NE.0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'get_xset_VerletLists',     &
            &                       'vlist',__LINE__,info)
                  GOTO 9999
              END IF

              !-----------------------------------------------------------------
              !  Fill in verlet lists depending on whether lists will be
              !  symmetric or not.
              !-----------------------------------------------------------------
              DO i = 1, red_clist%n_real_p
                  p_idx = red_clist%rank(i)
                  CALL get_xset_neigh(p_idx, red_clist, blue_clist, &
 &                     whole_domain,  &
 &                     red, rcred, blue, rcblue, skin, vlist, nvlist)
              END DO
          END IF

9999      CONTINUE
          CALL substop('get_xset_VerletLists',t0,info)

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE get_xset_VerletLists_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE get_xset_VerletLists_d
#endif
#endif
