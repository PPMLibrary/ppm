      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_part_modify_add
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

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_part_modify_add_s(topoid,xp,Npart,Mpart,xpn,Nnew,&
              isymm,ghostsize,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_modify_add_d(topoid,xp,Npart,Mpart,xpn,Nnew,&
              isymm,ghostsize,info)
#endif 
      !!! This routine adds new particles to xp and updates the ghost mappings
      !!! and the corresponding buffers.
      !!! This routine is similar to the partial mapping routine 
      !!! (`ppm_map_part_partial`) in the sense that the ghost particles are
      !!! assumed to be located on neighbouring processors only, and thus only
      !!! require a nearest neighbour communication.
      !!!
      !!! [WARNING]
      !!! It is an error to specify a topology ID to which the particles are not
      !!! currently mapped
      !!!
      !!!
      !!! 1. determine which new particles are ghost particles,
      !!!    and which are real particles.
      !!! 2. insert new real particles between xp(:,Np) and xp(:,Np+1)
      !!! 3. append new ghost particles at the end of xp
      !!! 4. get potential ghosts within the new real particles
      !!! 5. update mappings and fill buffers
      !!!
      !!! The ghosts are found in three steps:
      !!!
      !!!  1. consider ghosts within the local processor
      !!!  2. consider ghosts on neighbouring processors
      !!!  3. consider ghosts on periodic images of neighbouring
      !!!     processors
      !!!
      !!! Comments:
      !!!
      !!! 1. Local ghost particles come in two types: ghost particles that
      !!!    exist because two sub domains touch, and ghosts across a
      !!!    periodic boundary. The first is automatically handled
      !!!    by the particles - and no copy or/and send is
      !!!    required. The second ghosts could be handled during
      !!!    the calculation of the interactions, but this would
      !!!    require a check for periodicity and no explicit
      !!!    ghosts from the processor itself. However, it seems
      !!!    more natural to have ghosts - irrespectively of their
      !!!    origin, so ghosts from periodicity (on the same
      !!!    processor) are copied here.
      !!! 2. is standard procedure
      !!! 3. is handled by copying particles that are adjacent to
      !!!    faces on a physical boundary to their image position.
      !!!
      !!! Now, item 2. and 3. can be treated within the same logic
      !!! whereas 1. require a bit of thought: to keep the source
      !!! compact, we could loop over all processors including
      !!! the local one and check for ghosts - the problem with
      !!! this procedure is, that we check for ghosts by comparing
      !!! the location of particles within the extended subs
      !!! boundaries (and not their true ghost layer - which
      !!! would result in more IFs than necessary). However,
      !!! because of this, looping over the processor itself
      !!! would find ghosts that are not really ghost - the only
      !!! ghosts that should be found in this step are those due
      !!! to periodicity.  The solution is to store the ghosts as
      !!! ighost() with nghost denoting the total number of
      !!! ghosts (excluding those due to periodicity) and
      !!! nghostplus the total number of ghosts. During the loop
      !!! over the processor itself we therefore only need to
      !!! loop from nghost+1,nghostplus.
      !!! ==============================================================


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_typedef
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_util_commopt
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! ID of current topology
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: xp
      !!! The position of the particles
      INTEGER                 , INTENT(INOUT) :: Npart
      !!! The number of particles (on the local processor)
      INTEGER                 , INTENT(INOUT) :: Mpart
      !!! The number of particles (including ghosts)
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xpn
      !!! The position of the new particles to be added
      INTEGER                 , INTENT(IN   ) :: Nnew
      !!! The number of new particles (on the local processor)
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      !!! The size of the ghost layer
      INTEGER                 , INTENT(IN   ) :: isymm
      !!! Indicator for the use of symmetry
      !!!
      !!! * isymm > 0 use symmetry
      !!! * isymm = 0 do not use symmetry
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3) :: ldu
      INTEGER               :: i,j,k,isub,lda,iadd,jadd
      INTEGER               :: nlist1,nlist2,nghost,nghostplus
      INTEGER               :: ipart,sendrank !,recvrank
      INTEGER               :: iopt,iset,ibuffer
      REAL(MK), DIMENSION(:,:), POINTER :: xt  => NULL()
      ! position of potential ghosts
      REAL(MK), DIMENSION(:,:), POINTER :: xt_offset => NULL()
      ! offset of pot. ghosts
      REAL(MK)                      :: xminf,yminf,zminf ! full domain
      REAL(MK)                      :: xmaxf,ymaxf,zmaxf ! full domain
      REAL(MK)                      :: xmini,ymini,zmini ! inner domain
      REAL(MK)                      :: xmaxi,ymaxi,zmaxi ! inner domain
      REAL(MK), DIMENSION(ppm_dim)  :: min_phys
      REAL(MK), DIMENSION(ppm_dim)  :: max_phys
      REAL(MK), DIMENSION(ppm_dim)  :: len_phys
      REAL(MK)                      :: t0
      INTEGER                       :: nbc,npbc,nsbc
      LOGICAL, DIMENSION(2*ppm_dim) :: lextra
      INTEGER, DIMENSION(2*ppm_dim) :: ibc
      LOGICAL                       :: valid
      CHARACTER(ppm_char)           :: mesg
      TYPE(ppm_t_topo),POINTER      :: topo => NULL()
      REAL(MK)                      :: eps
      INTEGER, DIMENSION(:), POINTER :: idx_xpnr
      !!! indices of the new real particles to be added (the list of potential
      !!! ghosts to be sent will be extracted from this index set).

         INTEGER             ,DIMENSION(:),POINTER :: ppm_psendbuffer_add => NULL()
         INTEGER                                        ::  ppm_nsendbuffer_add
         INTEGER                                        ::  ppm_nrecvbuffer_add
         INTEGER                                        ::  ppm_sendbufsize_add
         INTEGER                                        ::  ppm_recvbufsize_add
         INTEGER                                        ::  ppm_buffer_set_add
         INTEGER            , DIMENSION(:),POINTER :: ppm_buffer2part_add => NULL()
         INTEGER            , DIMENSION(:),POINTER :: ppm_buffer_type_add => NULL()
         INTEGER            , DIMENSION(:),POINTER :: ppm_buffer_dim_add => NULL()
         REAL(ppm_kind_single),DIMENSION(:),POINTER::ppm_ghost_offsets_add => NULL()
         REAL(ppm_kind_double),DIMENSION(:),POINTER::ppm_ghost_offsetd_add => NULL()
         REAL(ppm_kind_single),DIMENSION(:),POINTER::ppm_sendbuffers_add => NULL()
         REAL(ppm_kind_double),DIMENSION(:),POINTER::ppm_sendbufferd_add => NULL()

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_part_modify_add',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 1) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
#if __KIND == __DOUBLE_PRECISION
      eps = ppm_myepsd
#else
      eps = ppm_myepss
#endif

      !-------------------------------------------------------------------------
      !  Split new particles into real and ghosts
      !  The split is stored as 2 index arrays (real and ghosts) in modify%
      !-------------------------------------------------------------------------
      CALL ppm_part_split_compute(topoid,xpn,Nnew,info)
      IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,'ppm_part_modify_add',  &
            &          'ppm_part_split_compute failed',__LINE__,info)
      ENDIF

      CALL ppm_part_split_apply(topoid,xp,ppm_dim,Npart,Mpart,xpn,Nnew,info)
      IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,'ppm_part_modify_add',  &
            &          'ppm_part_split_apply failed',__LINE__,info)
      ENDIF

      !-------------------------------------------------------------------------
      !  Count Boundary conditions
      !-------------------------------------------------------------------------
      npbc = 0
      nsbc = 0
      nbc = 0
      ibc(:) = 0
      lextra(:) = .FALSE.
      
      DO i=1,2*ppm_dim
          SELECT CASE (topo%bcdef(i))
          CASE (ppm_param_bcdef_periodic)
              npbc = npbc + 1
              ibc(i) = ppm_param_bcdef_periodic
          CASE (ppm_param_bcdef_symmetry)
              nsbc = nsbc + 1
              lextra(i) = .TRUE.
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_antisymmetry)
              nsbc = nsbc + 1
              lextra(i) = .TRUE.
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_neumann)
              nsbc = nsbc + 1
              lextra(i) = .TRUE.
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_dirichlet)
              nsbc = nsbc + 1
              lextra(i) = .TRUE.
              ibc(i) = ppm_param_bcdef_symmetry
          END SELECT
      ENDDO
      nbc = npbc + nsbc

      ! if there is still some data left in the buffer, warn the user
      IF (ppm_buffer_set .GT. 0) THEN
        info = ppm_error_warning
        CALL ppm_error(ppm_err_map_incomp,'ppm_part_modify_add',  &
     &       'Buffer was not empty. Possible loss of data!',__LINE__,info)
      ENDIF

      !----------------------------------------------------------------------
      !  first check if the optimal communication protocol is known
      !----------------------------------------------------------------------
      IF (.NOT.topo%isoptimized) THEN
        !-------------------------------------------------------------------
        !  if not: determine it
        !-------------------------------------------------------------------
        CALL ppm_util_commopt(topoid,info)
        IF (ppm_debug .GT. 1) THEN
            DO i=1,topo%nneighproc
                WRITE(mesg,'(A,I4)') 'have neighbor: ', topo%ineighproc(i)
                CALL ppm_write(ppm_rank,'ppm_part_modify_add',mesg,info)
            ENDDO
            DO i=1,topo%ncommseq
                WRITE(mesg,'(A,I4)') 'communicate: ', topo%icommseq(i)
                CALL ppm_write(ppm_rank,'ppm_part_modify_add',mesg,info)
            ENDDO
        ENDIF
      ENDIF

      ! for now here, but this can be done better
      ! TODO: kepp this consistant throughout simulation
      ppm_sendbufsize = 0
      IF (ppm_kind.EQ.ppm_kind_double) THEN
         IF (ASSOCIATED(ppm_sendbufferd_add)) THEN
            ppm_sendbufsize = SIZE(ppm_sendbufferd_add)
         ENDIF
      ELSE
         IF (ASSOCIATED(ppm_sendbuffers_add)) THEN
            ppm_sendbufsize = SIZE(ppm_sendbuffers_add)
         ENDIF
      ENDIF
      !-------------------------------------------------------------------------
      !  Compute the size of the computational box on this topology
      !-------------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
      min_phys(:) = topo%min_physd(:)
      max_phys(:) = topo%max_physd(:)
#else
      min_phys(:) = topo%min_physs(:)
      max_phys(:) = topo%max_physs(:)
#endif
      len_phys(:) = max_phys(:) - min_phys(:)

      !-------------------------------------------------------------------------
      !  Save the map type for the subsequent calls 
      !-------------------------------------------------------------------------
      ppm_map_type = ppm_param_map_ghost_get

      !-------------------------------------------------------------------------
      !  Allocate memory for the list of particle on the local processor that
      !  may be ghosts on other processors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = Nnew
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'list1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'list2',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ighost,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'ighost',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for their positions
      !-------------------------------------------------------------------------
      ldu(1) = ppm_dim
      ldu(2) = modify%Nrnew
      CALL ppm_alloc(xt,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'xt',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(xt_offset,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'xt',__LINE__,info)
          GOTO 9999
      ENDIF
      ldu(1) = modify%Nrnew
      CALL ppm_alloc(modify%idx_ghost_send_new,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'modify%idx_ghost_send_new',__LINE__,info)
          GOTO 9999
      ENDIF
      idx_xpnr => modify%idx_ghost_send_new


      !-------------------------------------------------------------------------
      !  List ilist1() holds the particles we are currently considering 
      !  List ilist2() holds the particles that have not yet been associated 
      !  with a sub, and ghost() holds the ghosts
      !-------------------------------------------------------------------------
      nlist1  = modify%Nrnew
      DO i=1,nlist1
         ilist1(i) = modify%idx_real_new(i)
      ENDDO
      nlist2  = 0
      nghost  = 0

      !-------------------------------------------------------------------------
      !  Fill the list with particles that are within the reach of the ghost
      !  regions of other subs. We do that by looping over the subs belonging
      !  to this processor and checking if the particles are well within the
      !  sub. If this is the case, the particle will never be a ghost 
      !-------------------------------------------------------------------------
      DO k=1,topo%nsublist
         !----------------------------------------------------------------------
         !  Initialize the second list counter to zero
         !----------------------------------------------------------------------
         nlist2 = 0

         !----------------------------------------------------------------------
         !  Get the (global) id of the sub
         !----------------------------------------------------------------------
         isub = topo%isublist(k)

         !----------------------------------------------------------------------
         !  Store the full extend of the sub in ?minf and ?maxf
         !----------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
         xminf = topo%min_subd(1,isub)
         xmaxf = topo%max_subd(1,isub)

         yminf = topo%min_subd(2,isub)
         ymaxf = topo%max_subd(2,isub)

         IF (ppm_dim.EQ.3) THEN
            zminf = topo%min_subd(3,isub)
            zmaxf = topo%max_subd(3,isub)
         ENDIF 
#else
         xminf = topo%min_subs(1,isub)
         xmaxf = topo%max_subs(1,isub)

         yminf = topo%min_subs(2,isub)
         ymaxf = topo%max_subs(2,isub)
         IF (ppm_dim.EQ.3) THEN
            zminf = topo%min_subs(3,isub)
            zmaxf = topo%max_subs(3,isub)
         ENDIF
#endif

         !----------------------------------------------------------------------
         !  Compute the size of the inner region
         !----------------------------------------------------------------------
         IF (isymm.GT.0) THEN
            !-------------------------------------------------------------------
            !  if we use symmetry the upper/right part of our sub will have 
            !  ghosts and therefore the particle at the upper/right part of the 
            !  sub cannot be ghosts on other processors. Thus the ghosts must 
            !  be found at the lower/left of the sub
            !-------------------------------------------------------------------
            
            xmini = xminf + ghostsize
            IF ((ABS(xmaxf - max_phys(1)).LT. eps).AND.lextra(2)) THEN
                xmaxi = xmaxf - ghostsize
            ELSE
                xmaxi = xmaxf
            ENDIF

            ymini = yminf + ghostsize
            IF ((ABS(ymaxf - max_phys(2)).LT. eps).AND.lextra(4)) THEN
                ymaxi = ymaxf - ghostsize
            ELSE
                ymaxi = ymaxf
            ENDIF
 
            IF (ppm_dim.EQ.3) THEN
               zmini = zminf + ghostsize
               IF ((ABS(zmaxf - max_phys(3)).LT. eps).AND.lextra(6)) THEN
                   zmaxi = zmaxf - ghostsize
               ELSE
                   zmaxi = zmaxf
               ENDIF
            ENDIF 
         ELSE
            !-------------------------------------------------------------------
            !  if we do not use symmetry, the particles along the entire 
            !  boundary of the sub will be ghosts on other processors
            !-------------------------------------------------------------------
            xmini = xminf + ghostsize
            xmaxi = xmaxf - ghostsize

            ymini = yminf + ghostsize
            ymaxi = ymaxf - ghostsize

            IF (ppm_dim.EQ.3) THEN
               zmini = zminf + ghostsize
               zmaxi = zmaxf - ghostsize
            ENDIF 
         ENDIF 

         !----------------------------------------------------------------------
         !  loop over the remaining particles
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xpn(1,ipart).GE.xminf.AND.xpn(1,ipart).LT.xmaxf.AND. &
     &             xpn(2,ipart).GE.yminf.AND.xpn(2,ipart).LT.ymaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, check if the particles are in the ghost layer 
                  !-------------------------------------------------------------
                  IF (xpn(1,ipart).LT.xmini.OR.xpn(1,ipart).GT.xmaxi.OR. &
     &                xpn(2,ipart).LT.ymini.OR.xpn(2,ipart).GT.ymaxi) THEN
                     !----------------------------------------------------------
                     !  if yes the particle will be a ghost somewhere
                     !----------------------------------------------------------
                     nghost         = nghost + 1   
                     ighost(nghost) = ipart
                     xt(1,nghost)   = xpn(1,ipart)
                     xt(2,nghost)   = xpn(2,ipart)
                     xt_offset(1,nghost) = 0.0_MK
                     xt_offset(2,nghost) = 0.0_MK
                     modify%idx_ghost_send_new(nghost)=ipart
                  ENDIF
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ELSE
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xpn(1,ipart).GE.xminf.AND.xpn(1,ipart).LT.xmaxf.AND. &
     &             xpn(2,ipart).GE.yminf.AND.xpn(2,ipart).LT.ymaxf.AND. &
     &             xpn(3,ipart).GE.zminf.AND.xpn(3,ipart).LT.zmaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, check if the particles is in the ghost layer 
                  !-------------------------------------------------------------
                  IF (xpn(1,ipart).LT.xmini.OR.xpn(1,ipart).GT.xmaxi.OR. &
     &                xpn(2,ipart).LT.ymini.OR.xpn(2,ipart).GT.ymaxi.OR. &
     &                xpn(3,ipart).LT.zmini.OR.xpn(3,ipart).GT.zmaxi) THEN
                     !----------------------------------------------------------
                     !  if yes the particle will be a ghost somewhere
                     !----------------------------------------------------------
                     nghost         = nghost + 1   
                     ighost(nghost) = ipart
                     xt(1,nghost)   = xpn(1,ipart)
                     xt(2,nghost)   = xpn(2,ipart)
                     xt(3,nghost)   = xpn(3,ipart)
                     xt_offset(1,nghost) = 0.0_MK
                     xt_offset(2,nghost) = 0.0_MK
                     xt_offset(3,nghost) = 0.0_MK
                     modify%idx_ghost_send_new(nghost)=ipart
                  ENDIF
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ENDIF 

         !----------------------------------------------------------------------
         !  swap the lists
         !----------------------------------------------------------------------
         DO j=1,nlist2
            ilist1(j) = ilist2(j)
         ENDDO
         nlist1 = nlist2

      ENDDO ! end of subs on local processor


      !-------------------------------------------------------------------------
      !  At the end the nlist2 should be zero
      !-------------------------------------------------------------------------
      IF (nlist2.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_part_unass,'ppm_part_modify_add',     &
     &       'nlist2 > 0',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the total number of ghosts incl. those due to boundary
      !  conditions
      !-------------------------------------------------------------------------
      modify%Ngsendnew = nghost
      nghostplus = nghost

      !-------------------------------------------------------------------------
      ! find out whether we have any non-freespace BCs
      !----------------------------------------------------------------------
      !  handle boundary conditions in x
      !----------------------------------------------------------------------
      IF (ibc(1).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_x.inc"
      ELSEIF (ibc(1).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_lx.inc"
      ENDIF  
      
      IF (ibc(2).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_ux.inc"
      ENDIF 
      !----------------------------------------------------------------------
      !  handle boundary conditions in y
      !----------------------------------------------------------------------
      IF (ibc(3).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_y.inc"
      ELSEIF (ibc(3).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_ly.inc"
      ENDIF
      IF (ibc(4).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_uy.inc"
      ENDIF
      !----------------------------------------------------------------------
      !  handle periodicity in z (if 3D)
      !----------------------------------------------------------------------
      IF (ppm_dim.EQ.3) THEN
         !-------------------------------------------------------------------
         !  yes, we split the if in two, since we do not know in what order
         !  the compiler will check and ppm_bcdef will only be allocated to
         !  four (4) in 2D
         !-------------------------------------------------------------------
         IF (ibc(5).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_z.inc"
         ELSEIF (ibc(5).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_lz.inc"
         ENDIF
         IF (ibc(6).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_uz.inc"
         ENDIF

      ENDIF ! of 3D

      !-------------------------------------------------------------------------
      !  Ok, now we have a shorter list of potential ghosts to search 
      !  allocate memory for the lghosts 
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nghostplus
      CALL ppm_alloc(lghost,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'logical list: lghost',__LINE__,info)
          GOTO 9999
      ENDIF

      !!-------------------------------------------------------------------------
      !!  Allocate memory for the global send/recv lists
      !!-------------------------------------------------------------------------
      !ldu(1) = topo%ncommseq
      !CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
      !IF (info .NE. 0) THEN
          !info = ppm_error_fatal
          !CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     !&        'global send rank list PPM_ISENDLIST',__LINE__,info)
          !GOTO 9999
      !ENDIF
      !CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
      !IF (info .NE. 0) THEN
          !info = ppm_error_fatal
          !CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     !&        'global recv rank list PPM_IRECVLIST',__LINE__,info)
          !GOTO 9999
      !ENDIF

      !-------------------------------------------------------------------------
      !  Allocate memory for the pointers to the particles that will be send
      !  and received from and to a processor. The ppm_recvbuffer is NOT used
      !  in this routine, but initialized from the precv() array in the routine
      !  ppm_map_part_send().
      !-------------------------------------------------------------------------
      ldu(1) = topo%ncommseq + 1
      CALL ppm_alloc(ppm_psendbuffer_add,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'global send buffer pointer ppm_psendbuffer_add',__LINE__,info)
          GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Step 1:
      !  First we find the ghosts due to periodicity on the local processor
      !-------------------------------------------------------------------------
      ppm_psendbuffer_add(1) = 1
      !ppm_nsendlist      = 0
      !ppm_nrecvlist      = 0
      iset               = 0
      ibuffer            = 0
      k                  = 1

      !-------------------------------------------------------------------------
      !  Get the rank of the processor
      !-------------------------------------------------------------------------
      sendrank = topo%icommseq(k) ! should be ppm_rank !
      !recvrank = sendrank

      !-------------------------------------------------------------------------
      !  Store the processor to which we will send to
      !-------------------------------------------------------------------------
      !ppm_nsendlist                = ppm_nsendlist + 1
      !ppm_isendlist(ppm_nsendlist) = sendrank

      !-------------------------------------------------------------------------
      !  Store the processor to which we will recv from
      !-------------------------------------------------------------------------
      !ppm_nrecvlist                = ppm_nrecvlist + 1
      !ppm_irecvlist(ppm_nrecvlist) = recvrank

      !-------------------------------------------------------------------------
      !  Store the number of buffer entries (this is the first)
      !-------------------------------------------------------------------------
      ppm_buffer_set = 1

      !-------------------------------------------------------------------------
      !  Allocate memory for the field registers that holds the dimension and
      !  type of the data
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_buffer_set
      CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'buffer dimensions PPM_BUFFER_DIM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &        'buffer types PPM_BUFFER_TYPE',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store the dimension and type
      !-------------------------------------------------------------------------
      ppm_buffer_dim(ppm_buffer_set)  = ppm_dim
#if    __KIND == __SINGLE_PRECISION
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#else
      ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Well we only do the stuff for real if we have periodicity !
      !  (Re)allocate memory for the send buffer. The required size of these 
      !  arrays is not easy to compute. The first step (step 1) require at
      !  most ppm_nsublist(topoid)*(nghostplus - nghost)*ppm_dim. The arrays
      !  are resized further below during step 2 and 3.
      !-------------------------------------------------------------------------
      IF (nbc.GT.0) THEN
         !----------------------------------------------------------------------
         !  We can grow or fit the arrays - a matter of taste 
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_grow
         ldu(1) = ppm_dim*(nghostplus - nghost)*topo%nsublist
         ppm_sendbufsize = ldu(1)
         !----------------------------------------------------------------------
         !  First allocate the sendbuffer
         !----------------------------------------------------------------------
         IF (ppm_kind.EQ.ppm_kind_double) THEN
            CALL ppm_alloc(ppm_sendbufferd_add,ldu,iopt,info)
            CALL ppm_alloc(ppm_ghost_offsetd_add,ldu,iopt,info)
         ELSE
            CALL ppm_alloc(ppm_sendbuffers_add,ldu,iopt,info)
            CALL ppm_alloc(ppm_ghost_offsets_add,ldu,iopt,info)
         ENDIF
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &           'global send buffer PPM_SENDBUFFER',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  then allocate the index list: buffer2part
         !----------------------------------------------------------------------
         ldu(1) = (nghostplus - nghost)*topo%nsublist
         CALL ppm_alloc(ppm_buffer2part_add,ldu,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &           'buffer-to-particles map ppm_buffer2part_add',__LINE__,info)
             GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  flag all ghosts as not yet taken
         !----------------------------------------------------------------------
         DO j=1,nghostplus
            lghost(j) = .TRUE.
         ENDDO

         !----------------------------------------------------------------------
         !  loop over the subs on the local processor
         !----------------------------------------------------------------------
         DO j=1,topo%nsublist
            !-------------------------------------------------------------------
            !  Get the global ID of the sub
            !-------------------------------------------------------------------
            isub = topo%isublist(j)

            !-------------------------------------------------------------------
            !  Define the extended resize of this sub 
            !-------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
               xminf = topo%min_subs(1,isub)
               xmaxf = topo%max_subs(1,isub)
   
               yminf = topo%min_subs(2,isub)
               ymaxf = topo%max_subs(2,isub)
   
               IF (ppm_dim.EQ.3) THEN
                  zminf = topo%min_subs(3,isub)
                  zmaxf = topo%max_subs(3,isub)
               ENDIF 
#else
               xminf = topo%min_subd(1,isub)
               xmaxf = topo%max_subd(1,isub)
   
               yminf = topo%min_subd(2,isub)
               ymaxf = topo%max_subd(2,isub)
   
               IF (ppm_dim.EQ.3) THEN
                  zminf = topo%min_subd(3,isub)
                  zmaxf = topo%max_subd(3,isub)
               ENDIF 
#endif
            IF (isymm.GT.0) THEN
               !----------------------------------------------------------------
               !  if we use symmetry ghosts will only be present at the
               !  upper/right part of the sub
               !----------------------------------------------------------------
               xmini = xminf  
               xmaxi = xmaxf + ghostsize
   
               ymini = yminf 
               ymaxi = ymaxf + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = zminf
                  zmaxi = zmaxf + ghostsize
               ENDIF

               !----------------------------------------------------------------
               ! If this subdomain is at the physical boundary, then we need an
               ! extra ghostlayer for handling the boundary conditions
               !----------------------------------------------------------------
               IF (ABS(xmini - min_phys(1)).LT. eps .AND. lextra(1)) THEN
                   xmini = xmini - ghostsize
               ENDIF
               IF (ABS(ymini - min_phys(2)).LT. eps .AND. lextra(3)) THEN
                   ymini = ymini - ghostsize
               ENDIF
               IF (ppm_dim.EQ.3) THEN
                   IF (ABS(zmini - min_phys(3)).LT. eps .AND. lextra(5)) THEN
                       zmini = zmini - ghostsize
                   ENDIF
               ENDIF 
            ELSE
               !----------------------------------------------------------------
               !  if we do not use symmetry, we have ghost all around
               !----------------------------------------------------------------
               xmini = xminf - ghostsize
               xmaxi = xmaxf + ghostsize
   
               ymini = yminf - ghostsize
               ymaxi = ymaxf + ghostsize
   
               IF (ppm_dim.EQ.3) THEN
                  zmini = zminf - ghostsize
                  zmaxi = zmaxf + ghostsize
               ENDIF 
            ENDIF 

            !-------------------------------------------------------------------
            !  Reallocate the arrays: ppm_sendbuffers_add/d and ppm_buffer2part_add. 
            !  The required size is the current size minus the current use plus 
            !  the maximum increment (nghostplus-nghost).
            !  We perform a grow_preserve to preserve the contents but only to 
            !  reallocate if an increased size is required.
            !-------------------------------------------------------------------
            iopt   = ppm_param_alloc_grow_preserve
            IF (ppm_kind.EQ.ppm_kind_double) THEN
               IF (ASSOCIATED(ppm_sendbufferd_add)) THEN
                  ldu(1) = ppm_sendbufsize + &
     &                     (nghostplus - nghost - ibuffer)*ppm_dim
               ELSE
                  ldu(1) = (nghostplus - ibuffer)*ppm_dim
               ENDIF
               CALL ppm_alloc(ppm_sendbufferd_add,ldu,iopt,info)
               CALL ppm_alloc(ppm_ghost_offsetd_add,ldu,iopt,info)
            ELSE
               IF (ASSOCIATED(ppm_sendbuffers_add)) THEN
                  ldu(1) = ppm_sendbufsize + &
     &                     (nghostplus - nghost - ibuffer)*ppm_dim
               ELSE
                  ldu(1) = (nghostplus - nghost - ibuffer)*ppm_dim 
               ENDIF 
               CALL ppm_alloc(ppm_sendbuffers_add,ldu,iopt,info)
               CALL ppm_alloc(ppm_ghost_offsets_add,ldu,iopt,info)
            ENDIF
            ppm_sendbufsize = ldu(1) 
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &              'global send buffer PPM_SENDBUFFER',__LINE__,info)
                GOTO 9999
            ENDIF

            IF (ASSOCIATED(ppm_buffer2part_add)) THEN
               ldu(1) = SIZE(ppm_buffer2part_add) + nghostplus - nghost - iset
            ELSE
               ldu(1) = nghostplus - nghost - iset
            ENDIF 
            CALL ppm_alloc(ppm_buffer2part_add,ldu,iopt,info)
            IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',     &
     &              'buffer-to-particles map ppm_buffer2part_add',__LINE__,info)
                GOTO 9999
            ENDIF

            !-------------------------------------------------------------------
            !  loop over the potential ghost particles due to periodicity
            !-------------------------------------------------------------------
            IF (ppm_dim.EQ.2) THEN
               !----------------------------------------------------------------
               !  Two dimensions
               !----------------------------------------------------------------
               DO i=nghost+1,nghostplus
                  !-------------------------------------------------------------
                  !  and check if it is inside the ghost region 
                  !-------------------------------------------------------------
                  IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND.lghost(i)) THEN
                     !----------------------------------------------------------
                     !  mark the ghost as taken
                     !----------------------------------------------------------
                     lghost(i) = .FALSE.

                     !----------------------------------------------------------
                     !  found one - increment the buffer counter
                     !----------------------------------------------------------
                     iset                  = iset + 1
   
                     !----------------------------------------------------------
                     !  store the ID of the particles
                     !----------------------------------------------------------
                     ppm_buffer2part_add(iset) = ighost(i)
   
                     !----------------------------------------------------------
                     !  store the particle
                     !----------------------------------------------------------
                     IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = REAL(xt(1,i),        &
     &                       ppm_kind_double)
                         ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(1,i), &
     &                       ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = REAL(xt(2,i),        &
     &                       ppm_kind_double)
                         ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(2,i), &
     &                       ppm_kind_double)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = xt(1,i)
                         ppm_ghost_offsetd_add(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = xt(2,i)
                         ppm_ghost_offsetd_add(ibuffer) = xt_offset(2,i)
#endif
                     ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = xt(1,i)
                         ppm_ghost_offsets_add(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = xt(2,i)
                         ppm_ghost_offsets_add(ibuffer) = xt_offset(2,i)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = REAL(xt(1,i),        &
     &                       ppm_kind_single)
                         ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(1,i), &
     &                       ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = REAL(xt(2,i),        &
     &                       ppm_kind_single)
                         ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(2,i), &
     &                       ppm_kind_single)
#endif
                     ENDIF
                  ENDIF 
               ENDDO
            ELSE
               !----------------------------------------------------------------
               !  Three dimensions
               !----------------------------------------------------------------
               DO i=nghost+1,nghostplus
                  !-------------------------------------------------------------
                  !  and check if it is inside the ghost region 
                  !-------------------------------------------------------------
                  IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. & 
     &                xt(3,i).GT.zmini.AND.xt(3,i).LT.zmaxi.AND.lghost(i)) THEN
                     !----------------------------------------------------------
                     !  mark the ghost as taken
                     !----------------------------------------------------------
                     lghost(i) = .FALSE.

                     !----------------------------------------------------------
                     !  found one - increment the buffer counter
                     !----------------------------------------------------------
                     iset                  = iset + 1

                     !----------------------------------------------------------
                     !  store the ID of the particles
                     !----------------------------------------------------------
                     ppm_buffer2part_add(iset) = ighost(i)

                     !----------------------------------------------------------
                     !  store the particle
                     !----------------------------------------------------------
                     IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = REAL(xt(1,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(1,i), &
     &                      ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = REAL(xt(2,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(2,i), &
     &                      ppm_kind_double)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = REAL(xt(3,i),        &
     &                      ppm_kind_double)
                         ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(3,i), &
     &                      ppm_kind_double)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = xt(1,i)
                         ppm_ghost_offsetd_add(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = xt(2,i)
                         ppm_ghost_offsetd_add(ibuffer) = xt_offset(2,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbufferd_add(ibuffer)   = xt(3,i)
                         ppm_ghost_offsetd_add(ibuffer) = xt_offset(3,i)
#endif
                     ELSE
#if    __KIND == __SINGLE_PRECISION
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = xt(1,i)
                         ppm_ghost_offsets_add(ibuffer) = xt_offset(1,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = xt(2,i)
                         ppm_ghost_offsets_add(ibuffer) = xt_offset(2,i)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = xt(3,i)
                         ppm_ghost_offsets_add(ibuffer) = xt_offset(3,i)
#else
                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = REAL(xt(1,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(1,i), &
     &                      ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = REAL(xt(2,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(2,i), &
     &                      ppm_kind_single)

                         ibuffer = ibuffer + 1 
                         ppm_sendbuffers_add(ibuffer)   = REAL(xt(3,i),        &
     &                      ppm_kind_single)
                         ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(3,i), &
     &                      ppm_kind_single)
#endif
                     ENDIF
                  ENDIF 
               ENDDO
            ENDIF ! 2/3 dimension 
           
         ENDDO ! end of loop over subs on local processor
      ENDIF ! of BC ghosts

      !-------------------------------------------------------------------------
      !  Update the buffer pointer (ie the current iset ... or the number of
      !  particle we will send to the k-th entry in the icommseq list
      !-------------------------------------------------------------------------
      ppm_psendbuffer_add(k+1) = iset + 1

      !-------------------------------------------------------------------------
      !  Step 2 and 3:
      !  Notice the ppm_icommseq(1,topoid) points to the local processor, and
      !  we already considered this on step 1, so we now start at k = 2
      !-------------------------------------------------------------------------
      DO k=2,topo%ncommseq
         !----------------------------------------------------------------------
         !  for each processor in the sendrank list we need to send the ghosts
         !  only once, so we flag the lghost as true initially
         !----------------------------------------------------------------------
         DO j=1,nghostplus
            lghost(j) = .TRUE.
         ENDDO

         !----------------------------------------------------------------------
         !  Get the rank of the processor
         !----------------------------------------------------------------------
         sendrank = topo%icommseq(k)
         !recvrank = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
         !ppm_nsendlist                = ppm_nsendlist + 1
         !ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
         !ppm_nrecvlist                = ppm_nrecvlist + 1
         !ppm_irecvlist(ppm_nrecvlist) = recvrank

         !----------------------------------------------------------------------
         !  only consider non-negative sendranks
         !----------------------------------------------------------------------
         IF (sendrank.GE.0) THEN
            !-------------------------------------------------------------------
            !  loop over all the subs in the current topology 
            !  this is a bit tedious and perhaps inefficient, but we do not
            !  have the subs belonging to all processors stored for each topoid
            !  sorry ! so we have to perform the calculations each time
            !-------------------------------------------------------------------
            DO j=1,topo%nsubs
               !----------------------------------------------------------------
               !  consider only those beloging to rank
               !----------------------------------------------------------------
               IF (topo%sub2proc(j).EQ.sendrank) THEN
                  !-------------------------------------------------------------
                  !  Define the extended resize of this sub 
                  !-------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
                  xminf = topo%min_subs(1,j)
                  xmaxf = topo%max_subs(1,j)

                  yminf = topo%min_subs(2,j)
                  ymaxf = topo%max_subs(2,j)

                  IF (ppm_dim.EQ.3) THEN
                     zminf = topo%min_subs(3,j)
                     zmaxf = topo%max_subs(3,j)
                  ENDIF 
#else
                  xminf = topo%min_subd(1,j)
                  xmaxf = topo%max_subd(1,j)

                  yminf = topo%min_subd(2,j)
                  ymaxf = topo%max_subd(2,j)
 
                  IF (ppm_dim.EQ.3) THEN
                     zminf = topo%min_subd(3,j)
                     zmaxf = topo%max_subd(3,j)
                  ENDIF 
#endif
                  IF (isymm.GT.0) THEN
                     !----------------------------------------------------------
                     !  if we use symmetry ghosts will only be present at the
                     !  upper/right part of the sub
                     !----------------------------------------------------------
                     xmini = xminf
                     xmaxi = xmaxf + ghostsize

                     ymini = yminf 
                     ymaxi = ymaxf + ghostsize

                     IF (ppm_dim.EQ.3) THEN
                        zmini = zminf 
                        zmaxi = zmaxf + ghostsize
                     ENDIF 
                     !----------------------------------------------------------
                     ! If we are at the border of the physical domain we have
                     ! to have ghostlayers either way
                     !----------------------------------------------------------
                     IF (ABS(xmini - min_phys(1)).LT.eps.AND.lextra(1)) THEN
                         xmini = xmini - ghostsize
                     ENDIF
                     IF (ABS(ymini - min_phys(2)).LT.eps.AND.lextra(3)) THEN
                         ymini = ymini - ghostsize
                     ENDIF
                     IF (ppm_dim.EQ.3) THEN
                       IF (ABS(zmini - min_phys(3)).LT.eps.AND.lextra(5)) THEN
                          zmini = zmini - ghostsize
                       ENDIF
                     ENDIF

                  ELSE
                     !----------------------------------------------------------
                     !  if we do not use symmetry, we have ghost all around
                     !----------------------------------------------------------
                     xmini = xminf - ghostsize
                     xmaxi = xmaxf + ghostsize

                     ymini = yminf - ghostsize
                     ymaxi = ymaxf + ghostsize

                     IF (ppm_dim.EQ.3) THEN
                        zmini = zminf - ghostsize
                        zmaxi = zmaxf + ghostsize
                     ENDIF 
                  ENDIF 

                  !-------------------------------------------------------------
                  !  Reallocate to make sure we have enough memory in the
                  !  ppm_buffer2part_add and ppm_sendbuffers_add/d 
                  !-------------------------------------------------------------
                  iopt   = ppm_param_alloc_grow_preserve
                  ldu(1) = ibuffer + nghostplus*ppm_dim*10 ! test 20061109
                  IF (ppm_kind.EQ.ppm_kind_double) THEN
                     IF ((ibuffer + nghostplus*ppm_dim).GT.ppm_sendbufsize) THEN
                        CALL ppm_alloc(ppm_sendbufferd_add,ldu,iopt,info)
                        CALL ppm_alloc(ppm_ghost_offsetd_add,ldu,iopt,info)
                        ppm_sendbufsize = ldu(1)
                     ENDIF
                  ELSE
                     IF ((ibuffer + nghostplus*ppm_dim).GT.ppm_sendbufsize) THEN
                        CALL ppm_alloc(ppm_sendbuffers_add,ldu,iopt,info)
                        CALL ppm_alloc(ppm_ghost_offsets_add,ldu,iopt,info)
                        ppm_sendbufsize = ldu(1)
                     ENDIF
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',   &
     &                    'global send buffer PPM_SENDBUFFER',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  ldu(1) = iset + nghostplus*10 ! test 20061109
                  IF ((iset + nghostplus).GT.SIZE(ppm_buffer2part_add)) THEN
                     CALL ppm_alloc(ppm_buffer2part_add,ldu,iopt,info)
                  ENDIF
                  IF (info .NE. 0) THEN
                      info = ppm_error_fatal
                      CALL ppm_error(ppm_err_alloc,'ppm_part_modify_add',   &
     &                   'buffer2particles map ppm_buffer2part_add',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  !-------------------------------------------------------------
                  !  loop over the potential ghost particles 
                  !-------------------------------------------------------------
                  IF (ppm_dim.EQ.2) THEN
                     !----------------------------------------------------------
                     !  Two dimensions
                     !----------------------------------------------------------
                     DO i=1,nghostplus
                        !-------------------------------------------------------
                        !  and check if it is inside the ghost region 
                        !-------------------------------------------------------
                        IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                      xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. &
     &                      lghost(i)) THEN
                           !----------------------------------------------------
                           !  Mark the ghost as taken for this sendrank
                           !----------------------------------------------------
                           lghost(i) = .FALSE.

                           !----------------------------------------------------
                           !  found one - increment the buffer counter
                           !----------------------------------------------------
                           iset = iset + 1

                           !----------------------------------------------------
                           !  store the ID of the particles
                           !----------------------------------------------------
                           ppm_buffer2part_add(iset) = ighost(i)

                           !----------------------------------------------------
                           !  store the particle
                           !----------------------------------------------------
                           IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = REAL(xt(1,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = REAL(xt(2,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_double)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = xt(1,i)
                               ppm_ghost_offsetd_add(ibuffer) = xt_offset(1,i)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = xt(2,i)
                               ppm_ghost_offsetd_add(ibuffer) = xt_offset(2,i)
#endif
                           ELSE
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = xt(1,i)
                               ppm_ghost_offsets_add(ibuffer) = xt_offset(1,i)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = xt(2,i)
                               ppm_ghost_offsets_add(ibuffer) = xt_offset(2,i)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = REAL(xt(1,i),        &
     &                            ppm_kind_single)
                               ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = REAL(xt(2,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_single)
#endif
                           ENDIF 
                        ENDIF 
                     ENDDO
                  ELSE
                     !----------------------------------------------------------
                     !  Three dimensions
                     !----------------------------------------------------------
                     DO i=1,nghostplus
                        !-------------------------------------------------------
                        !  and check if it is inside the ghost region 
                        !-------------------------------------------------------
                        IF (xt(1,i).GT.xmini.AND.xt(1,i).LT.xmaxi.AND. &
     &                      xt(2,i).GT.ymini.AND.xt(2,i).LT.ymaxi.AND. & 
     &                      xt(3,i).GT.zmini.AND.xt(3,i).LT.zmaxi.AND. &
     &                      lghost(i)) THEN

                           !----------------------------------------------------
                           !  Mark the ghost as taken for this sendrank
                           !----------------------------------------------------
                           lghost(i) = .FALSE.

                           !----------------------------------------------------
                           !  found one - increment the buffer counter
                           !----------------------------------------------------
                           iset                  = iset + 1

                           !----------------------------------------------------
                           !  store the ID of the particles
                           !----------------------------------------------------
                           ppm_buffer2part_add(iset) = ighost(i)

                           !----------------------------------------------------
                           !  store the particle
                           !----------------------------------------------------
                           IF (ppm_kind.EQ.ppm_kind_double) THEN
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = REAL(xt(1,i),          &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer) = REAL(xt(2,i),          &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_double)

                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = REAL(xt(3,i),        &
     &                            ppm_kind_double)
                               ppm_ghost_offsetd_add(ibuffer) = REAL(xt_offset(3,i), &
     &                            ppm_kind_double)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = xt(1,i)
                               ppm_ghost_offsetd_add(ibuffer) = xt_offset(1,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = xt(2,i)
                               ppm_ghost_offsetd_add(ibuffer) = xt_offset(2,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbufferd_add(ibuffer)   = xt(3,i)
                               ppm_ghost_offsetd_add(ibuffer) = xt_offset(3,i)
#endif
                           ELSE
#if    __KIND == __SINGLE_PRECISION
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = xt(1,i)
                               ppm_ghost_offsets_add(ibuffer) = xt_offset(1,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = xt(2,i)
                               ppm_ghost_offsets_add(ibuffer) = xt_offset(2,i)
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = xt(3,i)
                               ppm_ghost_offsets_add(ibuffer) = xt_offset(3,i)
#else
                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = REAL(xt(1,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(1,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = REAL(xt(2,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(2,i), &
     &                            ppm_kind_single)

                               ibuffer = ibuffer + 1 
                               ppm_sendbuffers_add(ibuffer)   = REAL(xt(3,i),   &
     &                            ppm_kind_single)
                               ppm_ghost_offsets_add(ibuffer) = REAL(xt_offset(3,i), &
     &                            ppm_kind_single)
#endif
                           ENDIF 
                        ENDIF 
                     ENDDO
                  ENDIF ! 2/3 dimension 
               ENDIF ! only consider subs belonging to rank
            ENDDO ! loop over all subs in topo
         ENDIF ! skipped negative ranks
         !----------------------------------------------------------------------
         !  Update the buffer pointer (ie the current iset ... or the number of
         !  particle we will send to the k-th entry in the icommseq list
         !----------------------------------------------------------------------
         ppm_psendbuffer_add(k+1) = iset + 1

      ENDDO ! loop over all processors in commseq

      !-------------------------------------------------------------------------
      !  Store the current size of the buffer
      !-------------------------------------------------------------------------
      ppm_nsendbuffer = ibuffer

      !-------------------------------------------------------------------------
      !  Fuse the actual ghost mappings with the ones we have just computed
      !  (and which map the newly added ghost particles to the corresponding
      !   real particles on the neighbouring processors)
      !-------------------------------------------------------------------------
      iset = 0
      DO k=topo%ncommseq,1,-1

          i=ppm_psendbuffer(k)
          j=ppm_psendbuffer(k+1)
          iadd=ppm_psendbuffer_add(k)
          jadd=ppm_psendbuffer_add(k+1)

          !copy to buffers and move data within buffers
          ppm_buffer2part(i+iadd-1:j+iadd-1) = & 
              ppm_buffer2part(i:j)
          ppm_buffer2part(i+iadd-1:jadd+iadd-1) = & 
              ppm_buffer2part_add(iadd:jadd)
              
          lda=ppm_dim
          IF (ppm_kind.EQ.ppm_kind_double) THEN
              ppm_ghost_offsetd(1+lda*(i+iadd-2):1+lda*(j+iadd-2)) = &
                  ppm_ghost_offsetd(1+lda*(i-1):1+lda*(j-1))
              ppm_ghost_offsetd(1+lda*(i+iadd-2):1+lda*(j+jadd-2)) = &
                  ppm_ghost_offsetd_add(1+lda*(iadd-1):1+lda*(jadd-1))
              ppm_sendbufferd(1+lda*(i+iadd-2):1+lda*(j+iadd-2)) = &
                  ppm_sendbufferd(1+lda*(i-1):1+lda*(j-1))
              ppm_sendbufferd(1+lda*(i+iadd-2):1+lda*(j+jadd-2)) = &
                  ppm_sendbufferd_add(1+lda*(iadd-1):1+lda*(jadd-1))
          ELSE
              ppm_ghost_offsets(1+lda*(i+iadd-2):1+lda*(j+iadd-2)) = &
                  ppm_ghost_offsets(1+lda*(i-1):1+lda*(j-1))
              ppm_ghost_offsets(1+lda*(i+iadd-2):1+lda*(j+jadd-2)) = &
                  ppm_ghost_offsets_add(1+lda*(iadd-1):1+lda*(jadd-1))
              ppm_sendbuffers(1+lda*(i+iadd-2):1+lda*(j+iadd-2)) = &
                  ppm_sendbuffers(1+lda*(i-1):1+lda*(j-1))
              ppm_sendbuffers(1+lda*(i+iadd-2):1+lda*(j+jadd-2)) = &
                  ppm_sendbuffers_add(1+lda*(iadd-1):1+lda*(jadd-1))
          ENDIF

         iset = iset + ppm_psendbuffer_add(k) !current offset
         ppm_psendbuffer(k) = ppm_psendbuffer(k) + iset
      ENDDO ! loop over all processors in commseq

      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      idx_xpnr => NULL()

      iopt = ppm_param_dealloc

      CALL ppm_alloc(    xt,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_part_modify_add',     &
     &       'xt',__LINE__,info)
      ENDIF

      CALL ppm_alloc(xt_offset,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_part_modify_add',     &
     &       'xt',__LINE__,info)
      ENDIF

      CALL ppm_alloc(lghost,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_part_modify_add',     &
     &       'lghost',__LINE__,info)
      ENDIF

      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_part_modify_add',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_part_modify_add',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .EQ. ppm_param_topo_undefined) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_no_topo,'ppm_part_modify_add',  &
     &            'This routine needs a topology defined topo',__LINE__,info)
              GOTO 8888
          ENDIF
          IF (topoid .NE. ppm_param_topo_undefined) THEN
            CALL ppm_check_topoid(topoid,valid,info)
            IF (.NOT. valid) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,'ppm_part_modify_add',  &
     &               'topoid out of range',__LINE__,info)
                GOTO 8888
            ENDIF
        ENDIF
        IF (Nnew .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_modify_add',  &
     &          'Nnew must be >=0',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (ghostsize .LT. 0.0_MK) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_part_modify_add',  &
     &          'ghostsize must be >=0.0',__LINE__,info)
            GOTO 8888
        ENDIF
 8888     CONTINUE
      END SUBROUTINE check
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_modify_add_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_modify_add_d
#endif
