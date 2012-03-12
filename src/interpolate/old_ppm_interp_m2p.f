      !-------------------------------------------------------------------------
      !     Subroutine   :                   ppm_interp_m2p
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


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_ss_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_ds_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_sv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_dv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_ss_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_ds_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_sv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_dv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif
#endif
      !!! This subroutine carries out mesh to particle interpolation.
      !!! 
      !!! Currently 2 interpolation schemes are supported:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      !!!
      !!! [WARNING]
      !!! This routine assumes that `up` is already allocated.
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
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_write
      USE ppm_module_check_id
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#endif
      !!! field from which to interpolate
#elif __MODE == __VEC
      INTEGER                         , INTENT(IN   )  :: lda
      !!! leading dimension of up
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER        :: field_up
#endif
      !!! field from which to interpolate
#endif
      REAL(MK), DIMENSION(:,:)       , INTENT(IN   ) :: xp
      !!! particle positions
      INTEGER , DIMENSION(:  )       , INTENT(IN   ) :: ghostsize
      !!! ghost size
      INTEGER                        , INTENT(IN   ) :: Np
      !!! number of particles.
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! topology identifier of target
      INTEGER                        , INTENT(IN   ) :: meshid
      !!! id of the mesh (user)
      INTEGER                        , INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.                    +
      !!! One of:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(MK),  DIMENSION(ppm_dim)          :: len_phys
      INTEGER                                :: kernel_support
      INTEGER                                :: nsubs
      INTEGER, DIMENSION(6)                  :: bcdef
      ! aliases
      REAL(mk), DIMENSION(:,:),      POINTER :: min_sub => NULL()
      REAL(mk), DIMENSION(:,:),      POINTER :: max_sub => NULL()
      REAL(mk)                               :: myeps
      REAL(mk)                               :: tim1s, tim1e
      TYPE(ppm_t_equi_mesh), POINTER         :: p_mesh => NULL()
      TYPE(ppm_t_topo)     , POINTER         :: topo   => NULL()


      ! Some arrays and variables defined as workaround tools for now
      REAL(mk), DIMENSION(:),        POINTER :: p_coor => NULL()
      REAL(mk), DIMENSION(:),        POINTER :: p_mash => NULL()
      REAL(mk), DIMENSION(:),        POINTER :: mesh_mash => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: minphys_with_ghost => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: maxphys_with_ghost => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: mesh_size => NULL()

      INTEGER                                :: nmesh
      INTEGER                                :: ndim_in
      INTEGER                                :: nmass_in
      INTEGER                                :: interp_kernel

      INTEGER                                :: i
      INTEGER                                :: j
      INTEGER                                :: k
      INTEGER                                :: dimn 

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_interp_m2p',t0,info)

#ifdef __GPU

      !-------------------------------------------------------------------------
      !  Some compilers have problems initializing constants that are declared
      ! in the module therefore we do it here:
      !-------------------------------------------------------------------------
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      dim = ppm_dim
      internal_weights = .FALSE.

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF

      IF(Np.EQ.0) GOTO 9999

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Get the mesh
      !-------------------------------------------------------------------------
      p_mesh => topo%mesh(meshid)
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------

      istart => p_mesh%istart

      !-------------------------------------------------------------------------
      !  Assignment of the useful arrays/scalar
      !-------------------------------------------------------------------------
      Nm(1:dim) = p_mesh%Nm
      bcdef(1:(2*dim)) = topo%bcdef(1:(2*dim))
      nsubs = topo%nsublist

      !-------------------------------------------------------------------------
      !  Allocate memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = Np
      IF(nsubs.NE.1) CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'particle list 1 ILIST1',__LINE__,info)
         GOTO 9999
      END IF
      IF(nsubs.NE.1) CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'particle list 2 ILIST2',__LINE__,info)
         GOTO 9999
      END IF

      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubs
      IF(nsubs.NE.1) CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                       'store_info allocation : problem',__LINE__,info)
         GOTO 9999
      END IF

      !-------------------------------------------------------------------------
      !  Mesh spacing
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_phys => topo%min_physs
      max_phys => topo%max_physs
#elif __KIND == __DOUBLE_PRECISION
      min_phys => topo%min_physd
      max_phys => topo%max_physd
#endif

      DO i = 1,dim
         Nc(i)       = Nm(i) - 1
         len_phys(i) = max_phys(i) - min_phys(i)
         dx(i)       = len_phys(i)/REAL(Nc(i),MK)
      END DO
      dxi = 1.0_mk/dx

      !-------------------------------------------------------------------------
      !  Initialize the particle list
      !-------------------------------------------------------------------------
      IF(nsubs.NE.1) THEN
         nlist1     = 0
         store_info = 0
         DO ipart=1,Np
            nlist1         = nlist1 + 1
            ilist1(nlist1) = ipart
         END DO
      END IF
#if   __KIND == __SINGLE_PRECISION
      myeps = ppm_myepss
      min_sub => topo%min_subs
      max_sub => topo%max_subs
#elif __KIND == __DOUBLE_PRECISION
      myeps = ppm_myepsd
      min_sub => topo%min_subd
      max_sub => topo%max_subd
#endif

      !-------------------------------------------------------------------------
      !  Loop over the subdomains (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !-------------------------------------------------------------------------
      IF(nsubs.NE.1) THEN
         DO idom = topo%nsublist,1,-1
            idoml = topo%isublist(idom)
            !-------------------------------------------------------------------
            !  Loop over the remaining particles
            !-------------------------------------------------------------------
            nlist2 = 0
            npart = 0
            DO i=1,nlist1
               ipart = ilist1(i)

               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
               IF(ppm_dim.EQ.3) THEN
                  !-------------------------------------------------------------
                  !  The particle is in the closure of the subdomain
                  !-------------------------------------------------------------
                  IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
     &                       xp(2,ipart).GE.min_sub(2,idoml) .AND. &
     &                       xp(3,ipart).GE.min_sub(3,idoml) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml) .AND. &
     &                       xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN

                     IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
     &                           (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
     &                           (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
     &                           (topo%subs_bc(6,idoml).EQ.1   .AND.    &
     &                           bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
     &                       xp(2,ipart).GE.min_sub(2,idoml) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN
                     IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
     &                           (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
     &                           (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

                        npart = npart + 1
                        store_info(idom) = npart
                     ELSE
                        nlist2         = nlist2 + 1
                        ilist2(nlist2) = ipart
                     END IF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  END IF
               END IF
            END DO ! end loop over remaining parts in ilist1
            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               END DO
            ENDIF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         END DO ! end loop over subs
         !----------------------------------------------------------------------
         !  Check that we sold all the particles
         !----------------------------------------------------------------------
         IF (nlist2.GT.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_part_unass,'ppm_interp_m2p',  &
     &                        'MAJOR PROBLEM',__LINE__,info)
            GOTO 9999
         END IF

         !----------------------------------------------------------------------
         !  Whats the maximum number of particles per subdomain
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO idom=1,topo%nsublist
            IF(store_info(idom).GE.max_partnumber) THEN
               max_partnumber = store_info(idom)
            END IF
         END DO
         iopt   = ppm_param_alloc_fit
         ldu(1) = topo%nsublist
         ldu(2) = max_partnumber

         !----------------------------------------------------------------------
         !  Allocate particle list
         !----------------------------------------------------------------------
         CALL ppm_alloc(list_sub,ldu,iopt,info)
         IF(info.NE.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF

         list_sub = 0
         !----------------------------------------------------------------------
         !  Initialize the particle list
         !----------------------------------------------------------------------
         nlist1     = 0
         DO ipart=1,Np
            nlist1         = nlist1 + 1
            ilist1(nlist1) = ipart
         ENDDO

         !----------------------------------------------------------------------
         !  Loop over the subdomains (since the first domains are most likely
         !  to be empty, we look backwards to reduce the number of elements in
         !  nlist2 as fast as possible)
         !----------------------------------------------------------------------
         DO idom = topo%nsublist,1,-1
            idoml = topo%isublist(idom)
            !-------------------------------------------------------------------
            !  loop over the remaining particles
            !-------------------------------------------------------------------
            nlist2 = 0
            npart = 0
            DO i=1,nlist1
               ipart = ilist1(i)
               !----------------------------------------------------------------
               !  If the particle is inside the current subdomain, assign it
               !----------------------------------------------------------------
               IF(ppm_dim.EQ.3) THEN
                  IF( ( xp(1,ipart).GE.min_sub(1,idoml) .AND. &
     &                       xp(2,ipart).GE.min_sub(2,idoml) .AND. &
     &                       xp(3,ipart).GE.min_sub(3,idoml) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml) .AND. &
     &                       xp(3,ipart).LE.max_sub(3,idoml) ) ) THEN

                     IF(   (xp(1,ipart).LT.max_sub(1,idoml) .OR.  &
     &                           (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
     &                           (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(3,ipart).LT.max_sub(3,idoml) .OR.  &
     &                           (topo%subs_bc(6,idoml).EQ.1   .AND.    &
     &                           bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
     &                       xp(2,ipart).GE.min_sub(2,idoml) .AND. &
     &                       xp(1,ipart).LE.max_sub(1,idoml) .AND. &
     &                       xp(2,ipart).LE.max_sub(2,idoml) ) ) THEN
                     IF(   (xp(1,ipart).LT.max_sub(1,idom) .OR.  &
     &                           (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                           bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                          (xp(2,ipart).LT.max_sub(2,idoml) .OR.  &
     &                           (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                           bcdef(4).NE. ppm_param_bcdef_periodic))) THEN
                        npart = npart + 1
                        list_sub(idom,npart) = ipart
                     ELSE
                        nlist2         = nlist2 + 1
                        ilist2(nlist2) = ipart
                     END IF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  END IF
               END IF
            END DO ! end loop over subs
            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               END DO
            END IF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT
         END DO

         !----------------------------------------------------------------------
         !  Check that we sold all the particles
         !----------------------------------------------------------------------
         IF (nlist2.GT.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_part_unass,'ppm_interp_m2p_3d',  &
     &                        'MAJOR PROBLEM',__LINE__,info)
            GOTO 9999
         ENDIF

         !----------------------------------------------------------------------
         !  Allocate and alias the weights if we need them.
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO idom = 1,topo%nsublist
            IF(store_info(idom).GE.max_partnumber) THEN
               max_partnumber = store_info(idom)
            END IF
         END DO
      ELSE ! now the case of one subdomain on this processor
         max_partnumber = np
         iopt = ppm_param_alloc_fit
         ldu(1) = 1
         CALL ppm_alloc(store_info,ldu,iopt,info)
         IF(info.NE.0) THEN
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF
         store_info(1) = max_partnumber
         ldu(1) = 1
         ldu(2) = max_partnumber
         CALL ppm_alloc(list_sub,ldu,iopt,info)
         IF(info.NE.0) THEN
            CALL ppm_error(ppm_err_alloc,'ppm_interp_m2p_3d',     &
     &                         'problem in internal allocation',__LINE__,info)
            GOTO 9999
         END IF
         DO i=1,max_partnumber
            list_sub(1,i) = i
         END DO
      END IF

      ! Set the arrays and variables accordingly.
      ndim_in  = ppm_dim
      nmass_in = lda

      ALLOCATE(minphys_with_ghost(ndim_in))
      ALLOCATE(maxphys_with_ghost(ndim_in))
      ALLOCATE(mesh_size(ndim_in))

      ! set interp_kernel parameter
      IF(kernel .EQ. ppm_param_rmsh_kernel_mp4) THEN
          interp_kernel = -1;
      ELSEIF(kernel .EQ. ppm_param_rmsh_kernel_bsp2) THEN
          interp_kernel = -2;
      ENDIF
      
      nmesh = nmass_in
      DO i = 1,ndim_in
          mesh_size(i) = Nm(i) + 2*ghostsize(i)
          nmesh = nmesh*mesh_size(i)
          minphys_with_ghost(i) = min_sub(i) - ghostsize(i)*dx(i)
          maxphys_with_ghost(i) = max_sub(i) + ghostsize(i)*dx(i)
      ENDDO

      ALLOCATE(mesh_mass(nmesh)) 
      ALLOCATE(p_coor(Np*ndim_in))
      ALLOCATE(p_mass(Np*nmass_in))

      IF(ndim_in .EQ. 2) THEN
          DO j = 1,mesh_size(2)
              DO i = 1,mesh_size(1)
                  mesh_mass(i + j*mesh_size(1)) = field_up(i,j,1)
              ENDDO
          ENDDO
      ELSEIF(ndim_in .EQ. 3) THEN
          DO k = 1,mesh_size(3)
              DO j = 1,mesh_size(2)
                  DO i = 1,mesh_size(1)
                      mesh_mass(i + j*mesh_size(1) + k*mesh_size(1)*   &
     &                      mesh_size(2)) = field_up(i,j,k,1)
                  ENDDO
              ENDDO
          ENDDO
      ENDIF

      DO i = 1, Np
          DO dimn = 1,ndim_in
              p_coor(ndim_in*i + dimn) = xp(dimn, i)
          ENDDO
      ENDDO

      ! Now we call the GPU code
!  ppm_gpu_interpolate_m2p_mp4(p_coor, p_mass, mesh_mass,        & 
!     &       minphys_with_ghost, maxphys_with_ghost, mesh_size,      &
!     &       Np, ndim_in, nmass_in, interp_kernel)

#endif

#ifndef __GPU
      cpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_interp_m2p',t0,info)
      RETURN
#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_dv_3d
#endif
#endif
#endif



