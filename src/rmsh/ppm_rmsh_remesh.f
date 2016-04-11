      !-------------------------------------------------------------------------
      !  Subroutine   :                     ppm_rmsh_remesh
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

      !-------------------------------------------------------------------------
      ! TODO:
      !  Have to modify field_up to be optional. if its not present then an
      !  internal field is allocated; Unfortunately a little bit of a pain cause
      !  the overloading is not going to get a grip on whats going on then.
      !-------------------------------------------------------------------------

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ss_2d(topoid,meshid,xp,Np,up,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ds_2d(topoid,meshid,xp,Np,up,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_sv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_dv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ss_3d(topoid,meshid,xp,Np,up,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_ds_3d(topoid,meshid,xp,Np,up,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_sv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_rmsh_remesh_dv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &     ghostsize,info,field_up,wx1_user,wx2_user,wx3_user)
#endif
#endif
#endif
      !!! This routine recompute the new particles values, thanks to
      !!! the precomputed weights...
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_map
      USE ppm_module_data
      USE ppm_module_data_rmsh
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: xp
      !!! Position of the particles
      INTEGER  , DIMENSION(:)         , INTENT(IN   )  :: ghostsize
      !!! Number of ghost particles
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      !!! The array to be remeshed
      REAL(MK)                                         :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
      !!! Mesh on which particles have been distributed (if field not
      !!! present)
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)     :: lda
      !!! Leading dimension
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      !!! The array to be remeshed
      REAL(MK) , DIMENSION(lda)                        :: tup
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif
      !!! Mesh on which particles have been distributed (if field not
      !!! present)
#endif
      INTEGER                         , INTENT(IN   )  :: Np
      !!! Number of particles (on processor)
      INTEGER                         , INTENT(IN   )  :: topoid
      !!! Topology identifier of target
      INTEGER                         , INTENT(IN   )  :: meshid
      !!! ID of the mesh
      INTEGER                         , INTENT(IN   )  :: kernel
      !!! Choice of the kernel used to compute the weights
      INTEGER                         , INTENT(  OUT)  :: info
      !!! Returns 0 upon success
      REAL(MK), OPTIONAL, DIMENSION(:,:,:),  POINTER   :: wx1_user
      !!! If present, these weights are used, instead of
      !!! internal weights
      REAL(MK), OPTIONAL, DIMENSION(:,:,:),  POINTER   :: wx2_user
      !!! If present, these weights are used, instead of
      !!! internal weights
      REAL(MK), OPTIONAL, DIMENSION(:,:,:),  POINTER   :: wx3_user
      !!! If present, these weights are used, instead of
      !!! internal weights

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)             :: len_phys
      REAL(MK), DIMENSION(:,:,:)   , POINTER   :: wx1
      REAL(MK), DIMENSION(:,:,:)   , POINTER   :: wx2
      REAL(MK), DIMENSION(:,:,:)   , POINTER   :: wx3
      REAL(MK), DIMENSION(:)       , POINTER   :: min_phys
      REAL(MK), DIMENSION(:)       , POINTER   :: max_phys
      REAL(MK), DIMENSION(ppm_dim)             :: dxi,dx
      REAL(MK)                                 :: dv1,dv2,dv3,epsilon
      INTEGER,  DIMENSION(ppm_dim)             :: Nc
      INTEGER,  DIMENSION(ppm_dim+2)           :: ldu,ldl
      REAL(MK)                                 :: kernel_support
      INTEGER                                  :: i,j,k,ii,jj,kk,iidec,l
      INTEGER                                  :: kkdec,ip1,ip2,ip3
      INTEGER                                  :: nb_sub,isub,dim,npart,isubl
      INTEGER                                  :: ifrom,ito,ip,max_partnumber
      INTEGER                                  :: ipart,nlist2,nlist1
      INTEGER                                  :: jjdec,iopt,idom,kkk
      INTEGER                                  :: ndata1_max,ndata2_max
      INTEGER                                  :: ndata3_max,nsubs,idim
      INTEGER,  DIMENSION(ppm_dim)             :: Nm
      INTEGER,  DIMENSION(:,:),      POINTER   :: ndata
      INTEGER,  DIMENSION(:,:),      POINTER   :: istart
      INTEGER,  DIMENSION(6)                   :: bcdef
      LOGICAL                                  :: consistent
      INTEGER                                  :: maptype
      TYPE(ppm_t_equi_mesh), POINTER           :: p_mesh
      TYPE(ppm_t_topo)     , POINTER           :: topo
      LOGICAL                                  :: valid
#if __MODE == __SCA
      INTEGER                                  :: lda = 1
#endif
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_rmsh_remesh',t0,info)

      !-------------------------------------------------------------------------
      !  This is a hack! Somehow having trouble with constructors in the
      !  ppm_module_data_rmsh module
      !-------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      dim = ppm_dim
      consistent = .TRUE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t

      SELECT TYPE (t => ppm_mesh%vec(meshid)%t)
      TYPE IS (ppm_t_equi_mesh)
         p_mesh => t
      END SELECT

      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------
      istart => p_mesh%istart

      IF(ppm_debug.GT.0) THEN
         !----------------------------------------------------------------------
         !  Check weights again
         !----------------------------------------------------------------------
         SELECT CASE(ppm_dim)

         CASE(2)
            IF(PRESENT(wx1_user)) THEN
               IF(   .NOT.ASSOCIATED(wx1_user).OR.&
     &               .NOT.ASSOCIATED(wx2_user) ) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &                'cannot remesh with empty weights',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDIF

         CASE(3)
            IF(PRESENT(wx1_user)) THEN
               IF(   .NOT.ASSOCIATED(wx1_user).OR.&
     &               .NOT.ASSOCIATED(wx2_user).OR.&
     &               .NOT.ASSOCIATED(wx3_user)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &                'cannot remesh with empty weights',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDIF

         END SELECT

         IF(.NOT.PRESENT(wx1_user)) THEN
#if   __KIND == __SINGLE_PRECISION
            IF(.NOT.ASSOCIATED(wx1_s)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'compute weights before remeshing',__LINE__,info)
               GOTO 9999
            ENDIF
#elif __KIND == __DOUBLE_PRECISION
            IF(.NOT.ASSOCIATED(wx1_d)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'compute weights before remeshing',__LINE__,info)
               GOTO 9999
            ENDIF
#endif
         ENDIF
      ENDIF

      IF(np.EQ.0) GOTO 5555
      !-------------------------------------------------------------------------
      !  Done with Argument checks, do aliasing
      !-------------------------------------------------------------------------
      IF(PRESENT(wx1_user)) THEN
         wx1 => wx1_user
         wx2 => wx2_user
         IF(ppm_dim.EQ.3) wx3 => wx3_user
      ELSE
#if   __KIND == __SINGLE_PRECISION
         wx1 => wx1_s
         wx2 => wx2_s
         IF(ppm_dim.EQ.3) wx3 => wx3_s
#elif __KIND == __DOUBLE_PRECISION
         wx1 => wx1_d
         wx2 => wx2_d
         IF(ppm_dim.EQ.3) wx3 => wx3_d
#endif
      ENDIF

 5555 CONTINUE
      !-------------------------------------------------------------------------
      !  Some more aliasing.
      !-------------------------------------------------------------------------
      nsubs                = topo%nsublist
      ! ndata(1:dim,1:nsubs) = ppm_cart_mesh(topoid,meshid)%nnodes
      ndata                => p_mesh%nnodes
      Nm(1:dim)            = p_mesh%Nm
      bcdef(1:2*ppm_dim)   = topo%bcdef(1:2*ppm_dim)

      !-------------------------------------------------------------------------
      !  Allocate memory for field if necessary
      !-------------------------------------------------------------------------
      ndata1_max = MAXVAL(ndata(1,topo%isublist(1:nsubs)))
      ndata2_max = MAXVAL(ndata(2,topo%isublist(1:nsubs)))
      IF(ppm_dim.EQ.3) ndata3_max=MAXVAL(ndata(3,topo%isublist(1:nsubs)))
      ! WARNING: distinguish between sca and vector allocation
      ! WARNING: distinguish between 2d and 3d allocation!!
      IF(ASSOCIATED(field_up)) THEN
         iopt    = ppm_param_alloc_grow
      ELSE
         iopt    = ppm_param_alloc_fit
      ENDIF
#if   __DIME == __3D
#if   __MODE == __VEC
      !-------------------------------------------------------------------------
      !  Allocation for 3D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1
      ldl(2)     = 1 - ghostsize(1)
      ldl(3)     = 1 - ghostsize(2)
      ldl(4)     = 1 - ghostsize(3)
      ldl(5)     = 1
      ldu(1)     = lda
      ldu(2)     = ndata1_max + ghostsize(1)
      ldu(3)     = ndata2_max + ghostsize(2)
      ldu(4)     = ndata3_max + ghostsize(3)
      ldu(5)     = nsubs
#elif __MODE == __SCA
      !-------------------------------------------------------------------------
      !  Allocation for 3D scalar
      !-------------------------------------------------------------------------
      ldl(1)     = 1 - ghostsize(1)
      ldl(2)     = 1 - ghostsize(2)
      ldl(3)     = 1 - ghostsize(3)
      ldl(4)     = 1
      ldu(1)     = ndata1_max + ghostsize(1)
      ldu(2)     = ndata2_max + ghostsize(2)
      ldu(3)     = ndata3_max + ghostsize(3)
      ldu(4)     = nsubs
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
      !-------------------------------------------------------------------------
      !  Allocation for 2D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1
      ldl(2)     = 1 - ghostsize(1)
      ldl(3)     = 1 - ghostsize(2)
      ldl(4)     = 1
      ldu(1)     = lda
      ldu(2)     = ndata1_max + ghostsize(1)
      ldu(3)     = ndata2_max + ghostsize(2)
      ldu(4)     = nsubs
#elif __MODE == __SCA
      !-------------------------------------------------------------------------
      !  Allocation for 2D vector
      !-------------------------------------------------------------------------
      ldl(1)     = 1 - ghostsize(1)
      ldl(2)     = 1 - ghostsize(2)
      ldl(3)     = 1
      ldu(1)     = ndata1_max + ghostsize(1)
      ldu(2)     = ndata2_max + ghostsize(2)
      ldu(3)     = nsubs
#endif
#endif
      CALL ppm_alloc(field_up,ldl,ldu,iopt,info)
      IF(info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_rmsh_remesh', &
     &       'field_up',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the field
      !-------------------------------------------------------------------------
      field_up = 0.0_MK
      IF(np.EQ.0) GOTO 9999
      !-------------------------------------------------------------------------
      !  Recover particle lists
      !-------------------------------------------------------------------------
      !  dont need to do this as they are in ppm_module_data_rmsh

      !-------------------------------------------------------------------------
      !  Pre computing operations
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys => topo%min_physs
      max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_phys => topo%min_physd
      max_phys => topo%max_physd
#endif
      DO i = 1,ppm_dim
         Nc(i)       = Nm(i)-1
         len_phys(i) = max_phys(i) - min_phys(i)
         dx(i)       = len_phys(i)/REAL(Nc(i),MK)
      ENDDO

      dxi = 1.0_MK/dx
      dv1 = dx(1)
      dv2 = dx(2)
      if(ppm_dim.EQ.3) dv3 = dx(3)

      !-------------------------------------------------------------------------
      !  Start computation, taking out the select case for speed
      !-------------------------------------------------------------------------
      SELECT CASE(kernel)
         !----------------------------------------------------------------------
         ! M prime 4
         !----------------------------------------------------------------------
      CASE(ppm_param_rmsh_kernel_mp4)
         DO isub = 1,nsubs
            isubl = topo%isublist(isub)
            DO ip=1,store_info(isub)
               ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1) &
     &                    )*dxi(1))+2-istart(1,isubl)
               ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2) &
     &                    )*dxi(2))+2-istart(2,isubl)
#if   __DIME == __3D
               ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3) &
     &                    )*dxi(3))+2-istart(3,isubl)
#endif
#if   __MODE == __VEC
               tup(1:lda) = up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
               tup        = up(list_sub(isub,ip))
#endif
               !----------------------------------------------------------------
               !  flipped upside down for ppm conformance
               !----------------------------------------------------------------
               if(lda.EQ.3.AND.ppm_dim.EQ.3) THEN
#if     __DIME == __3D
               kkdec = 0
               DO kk=ip3-1,ip3+2
                  kkdec = kkdec + 1
#endif
                  jjdec = 0
                  DO jj=ip2-1,ip2+2
                     jjdec = jjdec + 1

                     iidec = 0
                     DO ii=ip1-1,ip1+2
                        iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
!                       field_up(1:lda,ii,jj,kk,isub)&
!                            &=field_up(1:lda,ii,jj,kk,isub)&
!                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
!                            &*wx3(kkdec,isub,ip)*tup(1:lda)
                        field_up(1,ii,jj,kk,isub) &
     &                      =field_up(1,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(1)
                        field_up(2,ii,jj,kk,isub) &
     &                      =field_up(2,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(2)
                        field_up(3,ii,jj,kk,isub) &
     &                      =field_up(3,ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup(3)


#elif __MODE == __SCA
                        field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *wx3(kkdec,isub,ip)*tup
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
                        field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *tup(1:lda)
#elif __MODE == __SCA
                        field_up(ii,jj,isub)=field_up(ii,jj,isub) &
     &                      +wx1(iidec,isub,ip)*wx2(jjdec,isub,ip) &
     &                      *tup
#endif
#endif
                     ENDDO
                  ENDDO
#if  __DIME == __3D
               ENDDO
#endif
            ELSEIF (lda.EQ.4.AND.ppm_dim.EQ.3) THEN
#if     __DIME == __3D
               kkdec = 0
               DO kk=ip3-1,ip3+2
                  kkdec = kkdec + 1
#endif
                  jjdec = 0
                  DO jj=ip2-1,ip2+2
                     jjdec = jjdec + 1

                     iidec = 0
                     DO ii=ip1-1,ip1+2
                        iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
!                       field_up(1:lda,ii,jj,kk,isub)&
!                            &=field_up(1:lda,ii,jj,kk,isub)&
!                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
!                            &*wx3(kkdec,isub,ip)*tup(1:lda)
                        field_up(1,ii,jj,kk,isub)&
                             &=field_up(1,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(1)
                        field_up(2,ii,jj,kk,isub)&
                             &=field_up(2,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(2)
                        field_up(3,ii,jj,kk,isub)&
                             &=field_up(3,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(3)
                        field_up(4,ii,jj,kk,isub)&
                             &=field_up(4,ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup(4)


#elif __MODE == __SCA
                        field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*wx3(kkdec,isub,ip)*tup
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
                        field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                          &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                          &*tup(1:lda)
#elif __MODE == __SCA
                        field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                             &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                             &*tup
#endif
#endif
                     ENDDO
                  ENDDO
#if  __DIME == __3D
               ENDDO
#endif
            ELSE
#if     __DIME == __3D
              kkdec = 0
              DO kk=ip3-1,ip3+2
                 kkdec = kkdec + 1

#endif
                 jjdec = 0
                 DO jj=ip2-1,ip2+2
                    jjdec = jjdec + 1

                    iidec = 0
                    DO ii=ip1-1,ip1+2
                       iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
                       DO l=1,lda
                          field_up(l,ii,jj,kk,isub)&
                               &=field_up(l,ii,jj,kk,isub)&
                               &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                               &*wx3(kkdec,isub,ip)*tup(l)
                       ENDDO

#elif __MODE == __SCA
                       field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*tup
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
                       field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                         &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                         &*tup(1:lda)
#elif __MODE == __SCA
                       field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*tup
#endif
#endif
                    ENDDO
                 ENDDO
#if  __DIME == __3D
              ENDDO
#endif
           ENDIF
           ENDDO
        ENDDO
        !-----------------------------------------------------------------------
        ! BS2
        !-----------------------------------------------------------------------
     CASE(ppm_param_rmsh_kernel_bsp2)
        DO isub = 1,nsubs
           isubl = topo%isublist(isub)
           DO ip=1,store_info(isub)
              ip1 = INT((xp(1,list_sub(isub,ip))-min_phys(1)&
                   &)*dxi(1))+2-istart(1,isubl)
              ip2 = INT((xp(2,list_sub(isub,ip))-min_phys(2)&
                   &)*dxi(2))+2-istart(2,isubl)
#if   __DIME == __3D
              ip3 = INT((xp(3,list_sub(isub,ip))-min_phys(3)&
                   &)*dxi(3))+2-istart(3,isubl)
#endif
              !-----------------------------------------------------------------
              !  flipped upside down for ppm conformance
              !-----------------------------------------------------------------
#if     __DIME == __3D
              kkdec = 0
              DO kk=ip3,ip3+1
                 kkdec = kkdec + 1
#endif
                 jjdec = 0
                 DO jj=ip2,ip2+1
                    jjdec = jjdec + 1

                    iidec = 0
                    DO ii=ip1,ip1+1
                       iidec = iidec + 1
#if __DIME == __3D
#if   __MODE == __VEC
                       field_up(1:lda,ii,jj,kk,isub)&
                            &=field_up(1:lda,ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
                       field_up(ii,jj,kk,isub)=field_up(ii,jj,kk,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*wx3(kkdec,isub,ip)*up(list_sub(isub,ip))
#endif
#elif __DIME == __2D
#if   __MODE == __VEC
                       field_up(1:lda,ii,jj,isub)=field_up(1:lda,ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*up(1:lda,list_sub(isub,ip))
#elif __MODE == __SCA
                       field_up(ii,jj,isub)=field_up(ii,jj,isub)&
                            &+wx1(iidec,isub,ip)*wx2(jjdec,isub,ip)&
                            &*up(list_sub(isub,ip))
#endif
#endif
                    ENDDO
                 ENDDO

#if  __DIME == __3D
              ENDDO
#endif
           ENDDO
        ENDDO
      END SELECT

      !-------------------------------------------------------------------------
      !  Now map the ghosts in order to get consistent values at the border of
      !  the subdomains.
      !-------------------------------------------------------------------------

      CALL ppm_map_field_ghost_put(topoid,meshid,ghostsize,info)
      IF (info .NE. 0) GOTO 9999

#if   __MODE == __SCA
      CALL ppm_map_field_push(topoid,meshid,field_up,info)
#elif __MODE == __VEC
      CALL ppm_map_field_push(topoid,meshid,field_up,lda,info)
#endif
      IF (info .NE. 0) GOTO 9999

      CALL ppm_map_field_send(info)
      IF (info .NE. 0) GOTO 9999

#if   __MODE == __SCA
      CALL ppm_map_field_pop(topoid,meshid,field_up,ghostsize,info)
#elif __MODE == __VEC
      CALL ppm_map_field_pop(topoid,meshid,field_up,lda,ghostsize,info)
#endif

     !--------------------------------------------------------------------------
     !  Return
     !--------------------------------------------------------------------------
     9999 CONTINUE
     CALL substop('ppm_rmsh_remesh',t0,info)
     RETURN
      CONTAINS
      SUBROUTINE check
        IF(ppm_dim.EQ.3) THEN
            IF(   (PRESENT(wx1_user)).AND.&
     &            (PRESENT(wx2_user)).AND.&
     &            (PRESENT(wx3_user))) THEN
               consistent = .TRUE.
            ELSEIF(PRESENT(wx1_user).AND.(.NOT.(PRESENT(wx2_user))).AND. &
     &          (.NOT.(PRESENT(wx3_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx2_user or wx3_user weights lacking',__LINE__,info)
               GOTO 8888
            ELSEIF(PRESENT(wx2_user).AND.(.NOT.(PRESENT(wx3_user))).AND. &
     &         (.NOT.(PRESENT(wx1_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx1_user or wx3_user weights lacking',__LINE__,info)
               GOTO 8888
            ELSEIF(PRESENT(wx3_user).AND.(.NOT.(PRESENT(wx1_user))).AND. &
     &         (.NOT.(PRESENT(wx2_user)))) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             ' wx1_user or wx2_user weights lacking',__LINE__,info)
               GOTO 8888
            ENDIF
         ELSE
           IF(PRESENT(wx3_user)) THEN
               !----------------------------------------------------------------
               !  We dont care, maybe issue a warning..
               !----------------------------------------------------------------
            ELSEIF (PRESENT(wx1_user).AND..NOT.PRESENT(wx2_user)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'wx2_user not supplied',__LINE__,info)
               GOTO 8888
            ELSEIF (PRESENT(wx2_user).AND..NOT.PRESENT(wx1_user)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh', &
     &             'wx1_user not supplied',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (Np .GT. 0) THEN
            IF (SIZE(xp,2) .LT. Np) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'not enough particles contained in xp',__LINE__,info)
               GOTO 8888
            ENDIF
            IF (SIZE(xp,1) .LT. ppm_dim) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'leading dimension of xp insufficient',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDIF
         IF (Np .LE. 0) THEN
            IF (Np .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'particles must be specified',__LINE__,info)
               GOTO 8888
            ENDIF
            GOTO 8888
         ENDIF

         !----------------------------------------------------------------------
         !  Check if kernel is implemented
         !----------------------------------------------------------------------
         IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &          'wrong kernel definition',__LINE__,info)
            GOTO 8888
         ENDIF

         !----------------------------------------------------------------------
         !  This can never happen, if it does then the library is sick
         !----------------------------------------------------------------------
         kernel_support = ppm_rmsh_kernelsize(kernel)*2
        !kernel_support = ppm_rmsh_kernelsize(kernel)*2
        !IF (.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
        !     & .OR.(kernel_support.EQ.6))) THEN
        !   info = ppm_error_error
        !   CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
        !        &     'wrong kernel support',__LINE__,info)
        !   GOTO 8888
        !ENDIF

         IF (SIZE(up).LT.Np) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &         'wrong up array size',__LINE__,info)
            GOTO 8888
         ENDIF
         !----------------------------------------------------------------------
         !  This is also wrong: Ghostsize can be bigger...
         !----------------------------------------------------------------------
         !DO idim=1,dim
         !   IF (ghostsize(idim).NE.ppm_rmsh_kernelsize(kernel)) THEN
         !      info = ppm_error_error
         !      CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
         !           &     'wrong size of ghostsize',__LINE__,info)
         !      GOTO 8888
         !   ENDIF
         !ENDDO
         DO idim=1,dim
            IF (ghostsize(idim).LT.(ppm_rmsh_kernelsize(kernel)-1)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'ghostsize too small for this kernel',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
         CALL ppm_check_topoid(topoid,valid,info)
         IF (.NOT. valid) THEN
             info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_rmsh_remesh',  &
     &             'topoid not pointing to a valid topology',__LINE__,info)
               GOTO 8888
         ENDIF
         IF (.NOT.ppm_mesh%exists(meshid)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_rmsh,remesh', &
                 'meshid is invalid',__LINE__,info)
            GOTO 8888
        ENDIF

 8888   CONTINUE
      END SUBROUTINE check
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ss_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_sv_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ss_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_sv_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_rmsh_remesh_dv_3d
#endif
#endif
#endif
