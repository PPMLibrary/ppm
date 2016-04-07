      !-------------------------------------------------------------------------
      !     Subroutine   :                   ppm_rmsh_comp_weights
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
      SUBROUTINE ppm_rmsh_comp_weights_s(topoid,meshid,xp,Np,kernel,info, &
     &                                   wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_comp_weights_d(topoid,meshid,xp,Np,kernel,info, &
     &                                   wx1_user,wx2_user,wx3_user)
#endif
      !!! This routine computes the weights of the particles
      !!! which need to be remeshed in 2D.
      !-------------------------------------------------------------------------
      !  INCLUDES
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      USE ppm_module_data_rmsh
      USE ppm_module_write
      USE ppm_module_util_time
      USE ppm_module_mapping_typedef
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !--------------------------------------------------------------------------
      ! Arguments
      !--------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:)       , POINTER       :: xp
      !!! Particle positions
      INTEGER                        , INTENT(IN   ) :: Np
      !!! Number of particles
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! Topology identifier of target
      INTEGER                        , INTENT(IN   ) :: meshid
      !!! ID of the mesh
      INTEGER                        , INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns 0 upon success
      REAL(MK), OPTIONAL, DIMENSION(:,:,:), POINTER  :: wx1_user
      !!! If present, this routine does nothing, because these
      !!! weights are used after for remeshing
      REAL(MK), OPTIONAL, DIMENSION(:,:,:), POINTER  :: wx2_user
      !!! If present, this routine does nothing, because these
      !!! weights are used after for remeshing
      REAL(MK), OPTIONAL, DIMENSION(:,:,:), POINTER  :: wx3_user
      !!! If present, this routine does nothing, because these
      !!! weights are used after for remeshing

      !--------------------------------------------------------------------------
      ! Local variables
      !--------------------------------------------------------------------------
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(MK), DIMENSION(:)       , POINTER :: min_phys
      REAL(MK), DIMENSION(:)       , POINTER :: max_phys
      REAL(MK),  DIMENSION(ppm_dim)           :: dxi,dx
      REAL(MK),  DIMENSION(ppm_dim)           :: len_phys
      REAL(MK)                               :: x1,x2,x3,epsilon
      INTEGER                                :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
      INTEGER,  DIMENSION(ppm_dim)           :: Nc
      INTEGER                                :: i,j,k,ii,jj,kk,iidec
      INTEGER                                :: jjdec,nb_sub,npart,ipart
      INTEGER                                :: kkdec,ip1,nbpt_z,nlist1
      INTEGER                                :: ip2,ip3,nbpt_x,nbpt_y,iface
      INTEGER                                :: isub,ifrom,ito,ip,iopt,isubl
      INTEGER                                :: max_partnumber,idom,nlist2,idoml
      INTEGER, DIMENSION(1)                  :: lda
      INTEGER, DIMENSION(ppm_dim)            :: Nm
      INTEGER                                :: nsubs
      INTEGER, DIMENSION(6)                  :: bcdef
      LOGICAL                                :: internal_weights
      ! aliases
      REAL(MK), DIMENSION(:,:),      POINTER :: min_sub
      REAL(MK), DIMENSION(:,:),      POINTER :: max_sub
      REAL(MK), DIMENSION(:,:,:)   , POINTER :: wx1
      REAL(MK), DIMENSION(:,:,:)   , POINTER :: wx2
      REAL(MK), DIMENSION(:,:,:)   , POINTER :: wx3
      REAL(MK)                               :: myeps
      REAL(MK)                               :: tim1s, tim1e
      CHARACTER(len=256)                     :: msg
      LOGICAL                                :: valid
      TYPE(ppm_t_equi_mesh), POINTER         :: p_mesh
      TYPE(ppm_t_topo)     , POINTER         :: topo


      !--------------------------------------------------------------------------
      !  Initialize
      !--------------------------------------------------------------------------


      !--------------------------------------------------------------------------
      !  This is a hack! Somehow having trouble with constructors in the
      !  ppm_module_data_rmsh module
      !--------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      NULLIFY(wx1,wx2,wx3)

      CALL substart('ppm_rmsh_comp_weights',t0,info)

      internal_weights = .FALSE.

      IF(ppm_dim.EQ.3) THEN
        IF(PRESENT(wx1_user).AND.(.NOT.(PRESENT(wx2_user))).AND. &
             & (.NOT.(PRESENT(wx3_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx2_user or wx3_user weights lacking',__LINE__,info)
           GOTO 9999
        ELSEIF(PRESENT(wx2_user).AND.(.NOT.(PRESENT(wx3_user))).AND. &
             & (.NOT.(PRESENT(wx1_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx1_user or wx3_user weights lacking',__LINE__,info)
           GOTO 9999
        ELSEIF(PRESENT(wx3_user).AND.(.NOT.(PRESENT(wx1_user))).AND. &
             & (.NOT.(PRESENT(wx2_user)))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &  ' wx1_user or wx2_user weights lacking',__LINE__,info)
           GOTO 9999
        END IF
      ELSEIF (ppm_dim.EQ.2) THEN
        IF(PRESENT(wx1_user).AND.(.NOT.PRESENT(wx2_user)).OR.&
             &PRESENT(wx2_user).AND.(.NOT.PRESENT(wx1_user))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights', &
                &  'a wx_user array is missing',__LINE__,info)
           GOTO 9999
        END IF
      END IF

     !--------------------------------------------------------------------------
     !  Check arguments
     !--------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      !--------------------------------------------------------------------------
      !  Check meshid and topoid validity
      !--------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL ppm_check_topoid(topoid,valid,info)
         IF (.NOT. valid) THEN
             info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
     &             'topoid not pointing to a valid topology',__LINE__,info)
               GOTO 9999
         ENDIF
         IF (.NOT. ppm_mesh%exists(meshid)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights', &
                 'mesh_id is invalid',__LINE__,info)
            GOTO 9999
        ENDIF
      ENDIF

      topo => ppm_topo(topoid)%t

      SELECT TYPE (t => ppm_mesh%vec(meshid)%t)
      TYPE IS (ppm_t_equi_mesh)
         p_mesh => t
      END SELECT
      !--------------------------------------------------------------------------
      !  Get istart
      !--------------------------------------------------------------------------
      istart => p_mesh%istart

      !--------------------------------------------------------------------------
      !  Assignment of the useful arrays/scalar
      !--------------------------------------------------------------------------
      Nm(1:ppm_dim) = p_mesh%Nm
      bcdef(1:(2*ppm_dim)) = topo%bcdef(1:(2*ppm_dim))
      nsubs = topo%nsublist

      !--------------------------------------------------------------------------
      !  Alloc memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !--------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Np
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights',     &
     &               'particle list 1 ILIST1',__LINE__,info)
       GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights',     &
     &               'particle list 2 ILIST2',__LINE__,info)
       GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubs
      CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights',     &
     &                'store_info allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF

     !--------------------------------------------------------------------------
     !  Mesh spacing
     !--------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys => topo%min_physs
      max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_phys => topo%min_physd
      max_phys => topo%max_physd
#endif

      DO i = 1,ppm_dim
        Nc(i)       = Nm(i) - 1
        len_phys(i) = max_phys(i) - min_phys(i)
        dx(i)       = len_phys(i)/REAL(Nc(i),MK)
      ENDDO
      dxi     = 1.0_MK/dx
      epsilon = 0.000001_mk

     !--------------------------------------------------------------------------
     !  Initialize the particle list
     !--------------------------------------------------------------------------
     nlist1     = 0
     store_info = 0

     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO
#if   __KIND == __SINGLE_PRECISION
     myeps = ppm_myepss
     min_sub => topo%min_subs
     max_sub => topo%max_subs
#elif __KIND == __DOUBLE_PRECISION
     myeps = ppm_myepsd
     min_sub => topo%min_subd
     max_sub => topo%max_subd
#endif

     nlist2 = 0
     !--------------------------------------------------------------------------
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !--------------------------------------------------------------------------
     DO idom = topo%nsublist,1,-1
        idoml = topo%isublist(idom)
        !-----------------------------------------------------------------------
        !  Loop over the remaining particles
        !-----------------------------------------------------------------------
        nlist2 = 0
        npart = 0
        DO i=1,nlist1
           ipart = ilist1(i)
           !--------------------------------------------------------------------
           !  If the particle is inside the current subdomain, assign it
           !--------------------------------------------------------------------
           IF(ppm_dim.EQ.3) THEN
              !-----------------------------------------------------------------
              ! the particle is in the closure of the subdomain
              !-----------------------------------------------------------------
              IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                   &xp(3,ipart).GE.min_sub(3,idoml) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml) .AND. &
                   &xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN

                 IF( (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                      & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                      & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
                      & (topo%subs_bc(6,idoml).EQ.1   .AND.    &
                      & bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

                    npart = npart + 1
                    store_info(idom) = npart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           ELSEIF (ppm_dim.EQ.2) THEN

              IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN

                 IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                      & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                      & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                    npart = npart + 1
                    store_info(idom) = npart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           END IF
        END DO
        !-----------------------------------------------------------------------
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !-----------------------------------------------------------------------
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF

        !-----------------------------------------------------------------------
        !  Exit if the list is empty
        !-----------------------------------------------------------------------
        IF (nlist1.EQ.0) EXIT
     END DO
     !--------------------------------------------------------------------------
     !  Check that we sold all the particles
     !--------------------------------------------------------------------------
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_rmsh_comp_weights',  &
             &       'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF

     !--------------------------------------------------------------------------
     ! select the particle for with several domain computation (french)
     !--------------------------------------------------------------------------

     max_partnumber = 0
     DO idom=1,topo%nsublist
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)
        END IF
     END DO
     iopt   = ppm_param_alloc_fit
     ldu(1) = topo%nsublist
     ldu(2) = max_partnumber

     CALL ppm_alloc(list_sub,ldu,iopt,info)
     IF(info.NE.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_rmsh_comp_weights',     &
     &        'problem in internal allocation',__LINE__,info)
        GOTO 9999
     END IF

     list_sub=0

     !--------------------------------------------------------------------------
     !  Initialize the particle list
     !--------------------------------------------------------------------------
     nlist1     = 0
     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO

     !--------------------------------------------------------------------------
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !--------------------------------------------------------------------------
     DO idom = topo%nsublist,1,-1
        idoml = topo%isublist(idom)
        !-----------------------------------------------------------------------
        !  loop over the remaining particles
        !-----------------------------------------------------------------------
        nlist2 = 0
        npart = 0

        DO i=1,nlist1
           ipart = ilist1(i)

           !--------------------------------------------------------------------
           !  If the particle is inside the current subdomain, assign it
           !--------------------------------------------------------------------
           IF(ppm_dim.EQ.3) THEN
              IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                   &xp(3,ipart).GE.min_sub(3,idoml) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml) .AND. &
                   &xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN

                 IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
                      & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                      & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
                      & (topo%subs_bc(6,idoml).EQ.1   .AND.    &
                      & bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN


                    npart = npart + 1
                    list_sub(idom,npart) = ipart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           ELSEIF (ppm_dim.EQ.2) THEN

              IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
                   &xp(2,ipart).GE.min_sub(2,idoml) .AND. &
                   &xp(1,ipart).LE.max_sub(1,idoml) .AND. &
                   &xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN

                 IF(   (xp(1,ipart).LT.max_sub(1,idom) .OR.  &
                      & (topo%subs_bc(2,idoml).EQ.1   .AND.    &
                      & bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &(xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
                      & (topo%subs_bc(4,idoml).EQ.1   .AND.    &
                      & bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                    npart = npart + 1
                    list_sub(idom,npart) = ipart
                 ELSE
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                 ENDIF
              ELSE
                 nlist2         = nlist2 + 1
                 ilist2(nlist2) = ipart
              END IF
           END IF

        END DO
        !-----------------------------------------------------------------------
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !-----------------------------------------------------------------------
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF

        !-----------------------------------------------------------------------
        !  Exit if the list is empty
        !-----------------------------------------------------------------------
        IF (nlist1.EQ.0) EXIT
     END DO

     !--------------------------------------------------------------------------
     !  Check that we sold all the particles
     !--------------------------------------------------------------------------
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_rmsh_comp_weights',  &
             &       'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF

     !--------------------------------------------------------------------------
     !  Allocate and alias the weights if we need them.
     !--------------------------------------------------------------------------
     max_partnumber = 0
     DO idom = 1,topo%nsublist
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)
        END IF
     END DO

     iopt = ppm_param_alloc_fit
     ldl(1) = 1
     ldl(2) = 1
     ldl(3) = 1 - ppm_rmsh_kernelsize(kernel)
     ldu(1) = ppm_rmsh_kernelsize(kernel)*2
     ldu(2) = topo%nsublist
     ldu(3) = max_partnumber + ppm_rmsh_kernelsize(kernel)

     IF(((PRESENT(wx1_user)).AND.PRESENT(wx2_user).AND. &
          & PRESENT(wx3_user).AND.ppm_dim.EQ.3).OR.&
          &(PRESENT(wx1_user).AND.PRESENT(wx2_user).AND.&
          ppm_dim.eq.2)) THEN
        internal_weights = .FALSE.

        IF (.NOT.(ASSOCIATED(wx1_user)).OR. &
             &   .NOT.(ASSOCIATED(wx2_user))) THEN
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx1_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx1_user => wx1_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx1_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx1_user => wx1_d
#endif
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx2_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx2_user => wx2_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx2_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx2_user => wx2_d
#endif
        END IF
        IF(PRESENT(wx3_user)) THEN
           IF(ppm_dim.EQ.3.AND.&
                &.NOT.ASSOCIATED(wx3_user)) THEN
#if   __KIND == __SINGLE_PRECISION
              CALL ppm_alloc(wx3_s,ldl,ldu,iopt,info)
              IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                      &            'pb in weights allocation',__LINE__,info)
                 GOTO 9999
              END IF

              wx3_user => wx3_s
#elif __KIND == __DOUBLE_PRECISION
              CALL ppm_alloc(wx3_d,ldl,ldu,iopt,info)
              IF (info.NE.0) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                  &            'pb in weights allocation',__LINE__,info)
                 GOTO 9999
              END IF
              wx3_user => wx3_d
#endif
           END IF
        END IF
        wx1 => wx1_user
        wx2 => wx2_user
        IF(ppm_dim.EQ.3) wx3 => wx3_user
     ELSE
#if   __KIND == __SINGLE_PRECISION
        CALL ppm_alloc(wx1_s,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
            &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx1 => wx1_s
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_alloc(wx1_d,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx1 => wx1_d
#endif
#if   __KIND == __SINGLE_PRECISION
        CALL ppm_alloc(wx2_s,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx2 => wx2_s
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_alloc(wx2_d,ldl,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &            'pb in weights allocation',__LINE__,info)
           GOTO 9999
        END IF
        wx2 => wx2_d
#endif
        IF(ppm_dim.EQ.3) THEN
#if   __KIND == __SINGLE_PRECISION
           CALL ppm_alloc(wx3_s,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx3 => wx3_s
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_alloc(wx3_d,ldl,ldu,iopt,info)
           IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &            'pb in weights allocation',__LINE__,info)
              GOTO 9999
           END IF
           wx3 => wx3_d
#endif
        END IF
     END IF
     !--------------------------------------------------------------------------
     !  Initialize Weights
     !--------------------------------------------------------------------------
     CALL ppm_util_time(tim1s)
     wx1 = 0.0_MK
     wx2 = 0.0_MK
     IF(ppm_dim.EQ.3) wx3 = 0.0_MK
     !--------------------------------------------------------------------------
     !  Beginning of the computation
     !--------------------------------------------------------------------------
     !  loop over subs
     SELECT CASE(kernel)

     CASE(ppm_param_rmsh_kernel_bsp2)

        DO isub = 1,topo%nsublist
           !  loop over particles therein
           DO ip = 1,store_info(isub)
              isubl = topo%isublist(isub)

              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1) &
     &             )*dxi(1))+2-istart(1,isubl)

              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2) &
     &             )*dxi(2))+2-istart(2,isubl)
              IF(ppm_dim.EQ.3) THEN
                 ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3) &
     &             )*dxi(3))+2-istart(3,isubl)
              END IF

              iidec = 0
              jjdec = 0
              kkdec = 0

              !-----------------------------------------------------------------
              !  Bsplines 2
              !-----------------------------------------------------------------

              DO ii    = ip1,ip1+1
                 iidec = iidec + 1
                 x1 = ABS(REAL(ii-2+istart(1,isubl),mk)-( &
     &               xp(1,list_sub(isub,ip))-min_phys(1))*dxi(1))

                 IF(x1.LE.1.0_MK) THEN
                    wx1(iidec,isub,ip) =  1.0_MK - x1
                 END IF
              END DO

              ! Loop over y direction
              DO jj    = ip2,ip2+1
                 jjdec = jjdec + 1
                 x2 = ABS(REAL(jj-2+istart(2,isubl),mk)-( &
     &               xp(2,list_sub(isub,ip))-min_phys(2))*dxi(2))

                 IF(x2.LE.1.0_MK) THEN
                    wx2(jjdec,isub,ip) =  1.0_MK - x2
                 END IF
              END DO
              IF(ppm_dim.EQ.3) THEN
                 ! Loop over z direction
                 DO kk    = ip3,ip3+1
                    kkdec = kkdec + 1
                    x3 = ABS(REAL(kk-2+istart(3,isubl),mk)- (&
     &                  xp(3,list_sub(isub,ip))-min_phys(3))*dxi(3))

                    IF(x3.LE.1.0_MK) THEN
                       wx3(kkdec,isub,ip) =  1.0_MK - x3
                    END IF
                 END DO
              END IF

           END DO

        END DO

     CASE(ppm_param_rmsh_kernel_mp4)

        !-----------------------------------------------------------------
        ! M Prime Four
        !-----------------------------------------------------------------
        DO isub = 1,topo%nsublist
           !  loop over particles therein
           DO ip = 1,store_info(isub)
              isubl = topo%isublist(isub)

              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1) &
                   &         )*dxi(1))+2-istart(1,isubl)

              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2) &
                   &         )*dxi(2))+2-istart(2,isubl)
              IF(ppm_dim.EQ.3) THEN
                 ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3) &
                   &         )*dxi(3))+2-istart(3,isubl)
              END IF

              iidec = 0
              jjdec = 0
              kkdec = 0

              DO ii = ip1-1,ip1+2
                 iidec = iidec + 1
                 x1 = ABS(REAL(ii-2+istart(1,isubl),mk)-( &
     &               xp(1,list_sub(isub,ip))-min_phys(1))*dxi(1))
                 IF(x1.LE.2.0_MK) THEN
                    IF(x1.LE.1.0_MK) THEN
                       wx1(iidec,isub,ip) = 1.0_MK - x1**2*(2.5_mk-1.5_MK*x1)
                    ELSE
                       wx1(iidec,isub,ip) = 2.0_MK + (-4.0_MK + &
     &                     (2.5_MK - 0.5_MK * x1)*x1)*x1
                    END IF
                 END IF
              END DO

              DO jj = ip2-1,ip2+2
                 jjdec = jjdec + 1
                 x2 = ABS(REAL(jj-2+istart(2,isubl),mk)-( &
     &               xp(2,list_sub(isub,ip))-min_phys(2))*dxi(2))

                 IF(x2.LE.2.0_MK) THEN
                    IF(x2.LE.1.0_MK) THEN
                       wx2(jjdec,isub,ip) = 1.0_MK - x2**2*(2.5_mk-1.5_MK*x2)
                    ELSE
                       wx2(jjdec,isub,ip) = 2.0_MK + (-4.0_MK + &
     &                     (2.5_MK - 0.5_MK * x2)*x2)*x2
                    END IF
                 END IF

              END DO

              IF(ppm_dim.EQ.3) THEN

                 DO kk    = ip3 - 1, ip3+2

                    kkdec = kkdec + 1
                    x3 = ABS(REAL(kk-2+istart(3,isubl),mk)- (&
     &                  xp(3,list_sub(isub,ip))-min_phys(3))*dxi(3))
                    IF(x3.LE.2.0_MK) THEN
                       IF(x3.LE.1.0_MK) THEN
                          wx3(kkdec,isub,ip) =  1.0_MK - x3**2*(2.5_MK - &
                               & 1.5_MK*x3)
                       ELSE
                          wx3(kkdec,isub,ip) =  2.0_MK + (-4.0_MK + &
                               & (2.5_MK - 0.5_MK*x3)*x3)*x3
                       END IF
                    END IF
                 END DO
              END IF
           END DO
        END DO
     END SELECT

     CALL ppm_util_time(tim1e)

     IF (ppm_debug.gt.0) THEN
        WRITE(msg,*) 'naked calc of weights took [ms]: ',1000.0*(tim1e-tim1s)
        CALL ppm_write(ppm_rank,'ppm_rmsh_comp_weights',msg,info)
     END IF

      !-------------------------------------------------------------------------
      ! Deallocation of the arrays....
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      lda(1) = 0
      CALL ppm_alloc(ilist1,lda,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
    &        'ppm_rmsh_comp_weights', &
    &        'pb in ilist1 deallocation',__LINE__,info)
         GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,lda,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &       'ppm_rmsh_comp_weights',    &
     &       'pb in ilist2 deallocation',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Now why on earth would you want to deallocate them?
      !-------------------------------------------------------------------------
#ifdef __COMPLETELY_LOCO
      IF (internal_weights) THEN
         CALL ppm_alloc(wx1,lda,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc, &
     &           'ppm_rmsh_comp_weights',   &
     &           'pb in weight deallocation',__LINE__,info)
            GOTO 9999
         END IF
         CALL ppm_alloc(wx2,lda,iopt,info)
         IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc, &
     &           'ppm_rmsh_comp_weights',   &
     &           'pb in weight deallocation',__LINE__,info)
            GOTO 9999
         END IF
         IF(ppm_dim.EQ.3) THEN
            CALL ppm_alloc(wx3,lda,iopt,info)
            IF (info.NE.0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_dealloc, &
     &             'ppm_rmsh_comp_weights',    &
     &             'pb in weight deallocation',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDIF
      ENDIF
#endif
      IF(.NOT.internal_weights) THEN
         NULLIFY(wx1_s,wx1_d)
         NULLIFY(wx2_s,wx2_d)
         NULLIFY(wx3_s,wx3_d)
      END IF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999  CONTINUE
      CALL substop('ppm_rmsh_comp_weights',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'not enough particles contained in xp',__LINE__,info)
              GOTO 8888
           ENDIF
           IF (SIZE(xp,1) .LT. ppm_dim) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'leading dimension of xp insufficient',__LINE__,info)
              GOTO 8888
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                   &  'particles must be specified',__LINE__,info)
              GOTO 8888
           END IF
           GOTO 8888
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &     'wrong kernel definition',__LINE__,info)
           GOTO 8888
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
             &  .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_rmsh_comp_weights',  &
                &     'wrong kernel support',__LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_rmsh_comp_weights_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_rmsh_comp_weights_d
#endif


