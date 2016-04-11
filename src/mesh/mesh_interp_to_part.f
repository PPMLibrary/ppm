      !-------------------------------------------------------------------------
      !     Subroutine   :                   ppm_interp_m2p
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


      SUBROUTINE equi_mesh_m2p(this,Part,Field,kernel,info)
      !!! This subroutine carries out mesh to particle interpolation.
      !!!
      !!! Currently 2 interpolation schemes are supported:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      !!!
      !-------------------------------------------------------------------------
      !  INCLUDES
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_topo_typedef
      USE ppm_module_interp_m2p
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh)                   :: this
      !!! Mesh on which the field is already discretized
      CLASS(ppm_t_particles_d_), INTENT(INOUT) :: Part
      !!! Particle set on which the field should be interpolated
      CLASS(ppm_t_field_),       INTENT(INOUT) :: Field
      !!! Field that we want to interpolate
      INTEGER,                   INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.                    +
      !!! One of:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      INTEGER,                   INTENT(  OUT)   :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(MK) , DIMENSION(:,:)    , POINTER :: xp => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(MK)                               :: x1,x2,x3,epsilon
      INTEGER                                :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
      INTEGER                                :: i,j,k,ii,jj,kk,iidec,maptype,Np
      INTEGER                                :: jjdec,nb_sub,npart,ipart
      INTEGER                                :: kkdec,ip1,nbpt_z,nlist1
      INTEGER                                :: ip2,ip3,nbpt_x,nbpt_y,iface
      INTEGER                                :: isub,ifrom,ito,ip,dim,iopt,isubl
      INTEGER                                :: max_partnumber,nlist2
      INTEGER                                :: nsubpatch,ipatch
      LOGICAL                                :: internal_weights,lok
      ! aliases
      REAL(MK)                               :: myeps
      REAL(MK)                               :: tim1s, tim1e
      TYPE(ppm_t_topo)     , POINTER         :: topo
      CLASS(ppm_t_subpatch_),POINTER         :: p

      REAL(MK) , DIMENSION(:      ) , POINTER :: up_1d
      REAL(MK) , DIMENSION(:,:    ) , POINTER :: up_2d
      REAL(MK) , DIMENSION(:,:    ) , POINTER :: dummy_2d
      REAL(MK) , DIMENSION(:,:,:  ) , POINTER :: dummy_3d
      REAL(MK) , DIMENSION(:,:,:,:) , POINTER :: dummy_4d

      start_subroutine("m2p")

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
      check_associated(<#Part%xp#>,"Position array xp not initialized")

      Np = Part%Npart
      CALL Part%get_xp(xp,info)
      or_fail("Get_xp failed")

      IF (Np.EQ.0) GOTO 9999

      check_associated(<#this%subpatch#>,&
      "Mesh not allocated. Call Mesh%create() first?")

      topo => ppm_topo(this%topoid)%t

      !-------------------------------------------------------------------------
      !  Number of subpatches for this processor
      !-------------------------------------------------------------------------
      nsubpatch = this%subpatch%nb

      !-------------------------------------------------------------------------
      !  Allocate memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !-------------------------------------------------------------------------
      IF (nsubpatch.GE.1) THEN
          iopt   = ppm_param_alloc_fit
          ldu(1) = Np
          CALL ppm_alloc(ilist1,ldu,iopt,info)
          or_fail_alloc("particle list 1 ILIST1")
          CALL ppm_alloc(ilist2,ldu,iopt,info)
          or_fail_alloc("particle list 2 ILIST2")

          iopt   = ppm_param_alloc_fit
          ldu(1) = nsubpatch
          CALL ppm_alloc(store_info,ldu,iopt,info)
          or_fail_alloc("store_info")
      ENDIF

      !-------------------------------------------------------------------------
      !  Initialize the particle list
      !-------------------------------------------------------------------------
      IF (nsubpatch.GE.1) THEN
         nlist1     = 0
         store_info = 0
         DO ipart=1,Np
            nlist1         = nlist1 + 1
            ilist1(nlist1) = ipart
         ENDDO
      ENDIF
      nlist2=0

      !-------------------------------------------------------------------------
      !  Loop over the subpatches in each subdomain
      !  (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !-------------------------------------------------------------------------
      IF (nsubpatch.GE.1) THEN
          p => this%subpatch%begin()
          ipatch = 1
          DO WHILE (ASSOCIATED(p))
             !-----------------------------------------------------------------
             !  Loop over the remaining particles
             !-----------------------------------------------------------------
             nlist2 = 0
             npart = 0
             DO i=1,nlist1
                ipart = ilist1(i)

                !--------------------------------------------------------------
                !  If the particle is inside the current subdomain, assign it
                !--------------------------------------------------------------
                SELECT CASE (ppm_dim)
                CASE (3)
                   !---------------------------------------------------------
                   !  The particle is in the region of influence of
                   !  the subpatch (its closure reduced by a ghostlayer)
                   !---------------------------------------------------------
                   IF (xp(1,ipart).GE.p%start_red(1) .AND. &
                   &   xp(2,ipart).GE.p%start_red(2) .AND. &
                   &   xp(3,ipart).GE.p%start_red(3) .AND. &
                   &   xp(1,ipart).LE.p%end_red(1)   .AND. &
                   &   xp(2,ipart).LE.p%end_red(2)   .AND. &
                   &   xp(3,ipart).LE.p%end_red(3)           ) THEN

                       IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                       &   (p%bc(2).GE.0   .AND.    &
                       &   p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &   (xp(2,ipart).LT.p%end(2) .OR.  &
                       &   (p%bc(4).GE.0   .AND.    &
                       &   p%bc(4).NE. ppm_param_bcdef_periodic)).AND.&
                       &   (xp(3,ipart).LT.p%end(3) .OR.  &
                       &   (p%bc(6).GE.0   .AND.    &
                       &   p%bc(6).NE. ppm_param_bcdef_periodic))   ) THEN

                           npart = npart + 1
                           store_info(ipatch) = npart
                       ELSE
                           nlist2         = nlist2 + 1
                           ilist2(nlist2) = ipart
                       ENDIF
                   ELSE
                       nlist2         = nlist2 + 1
                       ilist2(nlist2) = ipart
                   ENDIF

                CASE (2)
                   IF ((xp(1,ipart).GE.p%start_red(1) .AND. &
                   &   xp(2,ipart).GE.p%start_red(2) .AND. &
                   &   xp(1,ipart).LE.p%end_red(1) .AND. &
                   &   xp(2,ipart).LE.p%end_red(2) ) ) THEN
                       IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                       &   (p%bc(2).GE.0   .AND.    &
                       &   p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                       &   (xp(2,ipart).LT.p%end(2) .OR.  &
                       &   (p%bc(4).GE.0   .AND.    &
                       &   p%bc(4).NE. ppm_param_bcdef_periodic))) THEN

                           npart = npart + 1
                           store_info(ipatch) = npart
                       ELSE
                           nlist2         = nlist2 + 1
                           ilist2(nlist2) = ipart
                       ENDIF
                   ELSE
                       nlist2         = nlist2 + 1
                       ilist2(nlist2) = ipart
                   ENDIF

                END SELECT

             ENDDO ! end loop over remaining parts in ilist1
             !-------------------------------------------------------------------
             !  Copy the lists (well, only if nlist2 changed - decreased)
             !-------------------------------------------------------------------
             IF (nlist2.NE.nlist1) THEN
                nlist1 = nlist2
                DO i=1,nlist1
                   ilist1(i) = ilist2(i)
                ENDDO
             ENDIF

             !-------------------------------------------------------------------
             !  Exit if the list is empty
             !-------------------------------------------------------------------
             IF (nlist1.EQ.0) EXIT
             p => this%subpatch%next()
             ipatch = ipatch + 1

          ENDDO !subpatch
         !----------------------------------------------------------------------
         !  Check whether we sold all the particles
         !----------------------------------------------------------------------
         IF (ppm_debug.GT.1) THEN
             IF (nlist2.GT.0) THEN
                stdout("Some particles seem to be outside from all subpatches",&
                &      " They will not take part in the m2p interpolation")
             ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Whats the maximum number of particles per subdomain
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO ipatch=1,nsubpatch
            IF (store_info(ipatch).GE.max_partnumber) THEN
               max_partnumber = store_info(ipatch)
            ENDIF
         ENDDO

         !----------------------------------------------------------------------
         !  Allocate particle list
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit
         ldu(1) = nsubpatch
         ldu(2) = max_partnumber
         CALL ppm_alloc(list_sub,ldu,iopt,info)
         or_fail_alloc("list_sub")

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
         !  Loop over the subpatches (since the first domains are most likely
         !  to be empty, we look backwards to reduce the number of elements in
         !  nlist2 as fast as possible)
         !----------------------------------------------------------------------
         p => this%subpatch%begin()
         ipatch = 1
         DO WHILE (ASSOCIATED(p))
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
               SELECT CASE (ppm_dim)
               CASE (3)
                  IF ((xp(1,ipart).GE.p%start_red(1) .AND. &
                  &    xp(2,ipart).GE.p%start_red(2) .AND. &
                  &    xp(3,ipart).GE.p%start_red(3) .AND. &
                  &    xp(1,ipart).LE.p%end_red(1) .AND. &
                  &    xp(2,ipart).LE.p%end_red(2) .AND. &
                  &    xp(3,ipart).LE.p%end_red(3) ) ) THEN

                     IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                     &   (p%bc(2).GE.0   .AND.    &
                     &   p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                     &   (xp(2,ipart).LT.p%end(2) .OR.  &
                     &   (p%bc(4).GE.0   .AND.    &
                     &   p%bc(4).NE. ppm_param_bcdef_periodic)).AND.&
                     &   (xp(3,ipart).LT.p%end(3) .OR.  &
                     &   (p%bc(6).GE.0   .AND.    &
                     &   p%bc(6).NE. ppm_param_bcdef_periodic))   ) THEN

                         npart = npart + 1
                         list_sub(ipatch,npart) = ipart
                     ELSE
                         nlist2         = nlist2 + 1
                         ilist2(nlist2) = ipart
                     ENDIF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  ENDIF

               CASE (2)
                  IF ((xp(1,ipart).GE.p%start_red(1) .AND. &
                  &   xp(2,ipart).GE.p%start_red(2) .AND. &
                  &   xp(1,ipart).LE.p%end_red(1) .AND. &
                  &   xp(2,ipart).LE.p%end_red(2) ) ) THEN
                     IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                     &   (p%bc(2).GE.0   .AND.    &
                     &   p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                     &   (xp(2,ipart).LT.p%end(2) .OR.  &
                     &   (p%bc(4).GE.0   .AND.    &
                     &   p%bc(4).NE. ppm_param_bcdef_periodic))) THEN
                        npart = npart + 1
                        list_sub(ipatch,npart) = ipart
                     ELSE
                        nlist2         = nlist2 + 1
                        ilist2(nlist2) = ipart
                     ENDIF
                  ELSE
                     nlist2         = nlist2 + 1
                     ilist2(nlist2) = ipart
                  ENDIF
               END SELECT
            ENDDO ! i
            !-------------------------------------------------------------------
            !  Copy the lists (well, only if nlist2 changed - decreased)
            !-------------------------------------------------------------------
            IF (nlist2.NE.nlist1) THEN
               nlist1 = nlist2
               DO i=1,nlist1
                  ilist1(i) = ilist2(i)
               ENDDO
            ENDIF

            !-------------------------------------------------------------------
            !  Exit if the list is empty
            !-------------------------------------------------------------------
            IF (nlist1.EQ.0) EXIT

            p => this%subpatch%next()
            ipatch = ipatch + 1

         ENDDO !supatch

         !----------------------------------------------------------------------
         !  Check if we sold all the particles
         !----------------------------------------------------------------------
         IF (ppm_debug.GT.1) THEN
             IF (nlist2.GT.0) THEN
                stdout("Some particles seem to be outside from all subpatches",&
                &      " They will not take part in the m2p interpolation")
             ENDIF
         ENDIF

         !----------------------------------------------------------------------
         !  Allocate and alias the weights if we need them.
         !----------------------------------------------------------------------
         max_partnumber = 0
         DO ipatch = 1,nsubpatch
            IF(store_info(ipatch).GE.max_partnumber) THEN
               max_partnumber = store_info(ipatch)
            ENDIF
         ENDDO
      ELSE ! now the case of ZERO subpatch on this processor
          !nothing to do
          IF (ppm_debug.GT.1) THEN
              stdout("No subpatch on this processor, skipping m2p interpolation")
              GOTO 9999
          ENDIF

      ENDIF

      !-------------------------------------------------------------------------
      !  Create the discretization on the particle set if it doesnt exist yet
      !-------------------------------------------------------------------------
      IF (.NOT.Field%is_discretized_on(Part)) THEN
          CALL Field%discretize_on(Part,info)
          or_fail("Could not discretize Field on Particle set")
      ENDIF
      !-------------------------------------------------------------------------
      !  Get a pointer to the data array and reset quantity to zero
      !-------------------------------------------------------------------------
      IF (Field%lda.EQ.1) THEN
         NULLIFY(up_1d)
         CALL Part%get(Field,up_1d,info)
         or_fail("Part%get_field")
         DO ip=1,np
            up_1d(ip) = 0.0_MK
         ENDDO
      ELSE
         NULLIFY(up_2d)
         CALL Part%get(Field,up_2d,info)
         or_fail("Part%get_field")

         SELECT CASE(Field%lda)
         CASE (1)
            DO ip=1,np
               up_2d(1,ip) = 0.0_MK
            ENDDO

         CASE (2)
            DO ip=1,np
               up_2d(1,ip) = 0.0_MK
               up_2d(2,ip) = 0.0_MK
            ENDDO

         CASE (3)
            DO ip=1,np
               up_2d(1,ip) = 0.0_MK
               up_2d(2,ip) = 0.0_MK
               up_2d(3,ip) = 0.0_MK
            ENDDO

         CASE (4)
            DO ip=1,np
               up_2d(1,ip) = 0.0_MK
               up_2d(2,ip) = 0.0_MK
               up_2d(3,ip) = 0.0_MK
               up_2d(4,ip) = 0.0_MK
            ENDDO

         CASE (5)
            DO ip=1,np
               up_2d(1,ip) = 0.0_MK
               up_2d(2,ip) = 0.0_MK
               up_2d(3,ip) = 0.0_MK
               up_2d(4,ip) = 0.0_MK
               up_2d(5,ip) = 0.0_MK
            ENDDO

         CASE DEFAULT
            up_2d = 0.0_MK
         END SELECT
      ENDIF

      !-------------------------------------------------------------------------
      !  Beginning of the computation
      !-------------------------------------------------------------------------

      SELECT CASE (kernel)
      CASE (ppm_param_rmsh_kernel_mp4)
         IF (Field%lda.EQ.1) THEN
            IF (ppm_dim .EQ. 2) THEN
               NULLIFY(dummy_2d)
               CALL m2p_interp_mp4(this,Field,dummy_2d,xp,up_1d,info)
            ELSE
               NULLIFY(dummy_3d)
               CALL m2p_interp_mp4(this,Field,dummy_3d,xp,up_1d,info)
            ENDIF
         ELSE
            IF (ppm_dim .EQ. 2) THEN
               NULLIFY(dummy_3d)
               CALL m2p_interp_mp4(this,Field,dummy_3d,Field%lda,xp,up_2d,info)
            ELSE
               NULLIFY(dummy_4d)
               CALL m2p_interp_mp4(this,Field,dummy_4d,Field%lda,xp,up_2d,info)
            ENDIF
         ENDIF
         or_fail("m2p_interp_mp4")

      CASE (ppm_param_rmsh_kernel_bsp2)
         IF (Field%lda.EQ.1) THEN
            IF (ppm_dim .EQ. 2) THEN
               NULLIFY(dummy_2d)
               CALL m2p_interp_bsp2(this,Field,dummy_2d,xp,up_1d,info)
            ELSE
               NULLIFY(dummy_3d)
               CALL m2p_interp_bsp2(this,Field,dummy_3d,xp,up_1d,info)
            ENDIF
         ELSE
            IF (ppm_dim .EQ. 2) THEN
               NULLIFY(dummy_3d)
               CALL m2p_interp_bsp2(this,Field,dummy_3d,Field%lda,xp,up_2d,info)
            ELSE
               NULLIFY(dummy_4d)
               CALL m2p_interp_bsp2(this,Field,dummy_4d,Field%lda,xp,up_2d,info)
            ENDIF
         ENDIF
         or_fail("m2p_interp_bsp2")

      CASE DEFAULT
         fail("Only Mp4 and BSp2 are avail. Use ppm_rmsh_remesh for other kernels.", &
         & ppm_err_argument)

      END SELECT ! kernel type

      !----------------------------------------------------------
      !  Dont deallocate if something went wrong, as
      !  the arrays may not even be associated
      !----------------------------------------------------------
      IF (info.NE.0) GOTO 9999
      !----------------------------------------------------------
      !  Deallocation of the arrays....
      !----------------------------------------------------------
      iopt = ppm_param_dealloc
      ldu(1) = 0
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      or_fail_dealloc("ilist1")
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      or_fail_dealloc("ilist2")
      CALL ppm_alloc(store_info,ldu,iopt,info)
      or_fail_dealloc("store_info")
      CALL ppm_alloc(list_sub,ldu,iopt,info)
      or_fail_dealloc("list_sub")

      !Updating internal state variables
      CALL Part%set_xp(xp,info,read_only=.TRUE.)
      or_fail("Set_xp failed")

      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
           fail("Please call ppm_init first!",ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           fail("Wrong kernel definition",ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF (.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
        &   .OR.(kernel_support.EQ.6))) THEN
           fail("Wrong kernel support",ppm_err_argument,exit_point=8888)
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE equi_mesh_m2p


