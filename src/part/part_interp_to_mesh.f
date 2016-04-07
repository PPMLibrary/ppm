      !------------------------------------------------------------------------!
      !     Subroutine   :                 ppm_interp_p2m
      !------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------!

      SUBROUTINE DTYPE(part_p2m)(this,Mesh,Field,kernel,info,p2m_bcdef)
      !!! This subroutine carries out particle to mesh interpolation.
      !!!
      !!! Currently 2 interpolation schemes are supported:
      !!!
      !!! * ppm_param_rmsh_kernel_bsp2
      !!! * ppm_param_rmsh_kernel_mp4
      !!!
      !!! [WARNING]
      !!! This routine assumes that `field_up` is already allocated.
      !!!
      !!! [TIP]
      !!! There is no need to perform a `ghost_get` before calling this routine
      !!! as the routine calls itself a `ghost_put` to the field after
      !!! interpolating from particles to the field.
      !------------------------------------------------------------------------!
      !  INCLUDES
      !------------------------------------------------------------------------!

      !------------------------------------------------------------------------!
      !  Modules
      !------------------------------------------------------------------------!
      USE ppm_module_interp_p2m
      IMPLICIT NONE

      DEFINE_MK()
      !-------------------------------------------------------------------------!
      ! Arguments
      !-------------------------------------------------------------------------!
      CLASS(DTYPE(ppm_t_particles))                  :: this
      CLASS(ppm_t_equi_mesh_),         INTENT(INOUT) :: Mesh
      CLASS(ppm_t_field_),             INTENT(INOUT) :: Field
      INTEGER,                         INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.
      INTEGER,                         INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      INTEGER, DIMENSION(:), OPTIONAL, POINTER       :: p2m_bcdef
      !!! Boundary conditions used for this interpolation routine, they may be
      !!! different than the boundary conditions set in the topology.
      !!! The values in this array can be one of:
      !!!
      !!! - ppm_param_bcdef_symmetry
      !!! - ppm_param_bcdef_antisymmetry
      !-------------------------------------------------------------------------!
      ! Local variables
      !-------------------------------------------------------------------------!
      CLASS(ppm_t_subpatch_), POINTER :: p

      CLASS(ppm_t_discr_data), POINTER :: prop

      REAL(MK), DIMENSION(:      ), POINTER :: up_1d
      REAL(MK), DIMENSION(:,:    ), POINTER :: up_2d
      REAL(MK), DIMENSION(:,:    ), POINTER :: dummy_2d
      REAL(MK), DIMENSION(:,:,:  ), POINTER :: dummy_3d
      REAL(MK), DIMENSION(:,:,:,:), POINTER :: dummy_4d
      REAL(MK), DIMENSION(:,:    ), POINTER :: xp

      INTEGER, DIMENSION(:), ALLOCATABLE :: ilist1,ilist2
      INTEGER                            :: kernel_support
      INTEGER, DIMENSION(ppm_dim+2)      :: ldu
      INTEGER                            :: i,j,k,l,nsubpatch,ipatch
      INTEGER                            :: npart,ipart
      INTEGER                            :: nlist1,nlist2
      INTEGER                            :: dim,iopt,max_partnumber

      LOGICAL :: internal_weights
      ! aliases

      start_subroutine("ppm_interp_p2m")

      !-------------------------------------------------------------------------!
      !  This is a hack! Somehow having trouble with constructors in the
      !  ppm_module_data_rmsh module
      !-------------------------------------------------------------------------!
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      dim = ppm_dim
      internal_weights = .FALSE.

      !-------------------------------------------------------------------------!
      !  Check arguments
      !-------------------------------------------------------------------------!
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      NULLIFY(xp)
      CALL this%get_xp(xp,info)
      or_fail("Get_xp failed")

      !-------------------------------------------------------------------------!
      !  If there is nothing to do, do nearly nothing
      !-------------------------------------------------------------------------!
      IF (this%Npart.EQ.0) GOTO 9998

      !-------------------------------------------------------------------------
      !  Number of subpatches for this processor
      !-------------------------------------------------------------------------
      nsubpatch = Mesh%subpatch%nb

      !-------------------------------------------------------------------------!
      !  Alloc memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !-------------------------------------------------------------------------!
      ldu(1) = this%Npart
      ALLOCATE(ilist1(ldu(1)),ilist2(ldu(1)),STAT=info)
      or_fail_alloc("ilist1 & ilist2")


      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubpatch
      CALL ppm_alloc(store_info,ldu,iopt,info)
      or_fail_alloc("store_info")

      !-------------------------------------------------------------------------!
      !  Initialize the particle list
      !-------------------------------------------------------------------------!
      store_info = 0

      FORALL (ipart=1:this%Npart) ilist1(ipart)=ipart

      nlist1=this%Npart
      nlist2=0

      !-------------------------------------------------------------------------!
      !  Loop over the subpatches in each subdomain
      !  (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !-------------------------------------------------------------------------!
      p => Mesh%subpatch%begin()
      ipatch = 1
      DO WHILE (ASSOCIATED(p))
        !----------------------------------------------------------------------!
        !  Loop over the remaining particles
        !----------------------------------------------------------------------!
        nlist2 = 0
        npart = 0
        DO i=1,nlist1
            ipart=ilist1(i)

            !-------------------------------------------------------------------!
            !  If the particle is inside the current subdomain, assign it
            !-------------------------------------------------------------------!
            IF (ppm_dim.EQ.3) THEN
                !-------------------------------------------------------------
                !  The particle is in the region of influence of the subpatch
                !  (its closure reduced by a ghostlayer)
                !-------------------------------------------------------------
                IF ( xp(1,ipart).GE.p%start_ext(1) .AND. &
                &    xp(2,ipart).GE.p%start_ext(2) .AND. &
                &    xp(3,ipart).GE.p%start_ext(3) .AND. &
                &    xp(1,ipart).LT.p%end_ext(1)   .AND. &
                &    xp(2,ipart).LT.p%end_ext(2)   .AND. &
                &    xp(3,ipart).LT.p%end_ext(3) ) THEN

                    IF ( (xp(1,ipart).LT.p%end(1)               .OR.  &
                    &    (p%bc(2).GE.0                          .AND. &
                    &     p%bc(2).NE. ppm_param_bcdef_periodic)).AND. &
                    &    (xp(2,ipart).LT.p%end(2)               .OR.  &
                    &    (p%bc(4).GE.0                          .AND. &
                    &     p%bc(4).NE. ppm_param_bcdef_periodic)).AND. &
                    &    (xp(3,ipart).LT.p%end(3)               .OR.  &
                    &    (p%bc(6).GE.0                          .AND. &
                    &     p%bc(6).NE. ppm_param_bcdef_periodic)) ) THEN

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
            ELSEIF (ppm_dim.EQ.2) THEN
                IF( ( xp(1,ipart).GE.p%start_ext(1) .AND. &
    &                      xp(2,ipart).GE.p%start_ext(2) .AND. &
    &                      xp(1,ipart).LT.p%end_ext(1) .AND. &
    &                      xp(2,ipart).LT.p%end_ext(2) ) ) THEN
                         IF(   (xp(1,ipart).LT.p%end(1) .OR.  &
    &                           (p%bc(2).GE.0   .AND.    &
    &                           p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
    &                          (xp(2,ipart).LT.p%end(2) .OR.  &
    &                           (p%bc(4).GE.0   .AND.    &
    &                           p%bc(4).NE. ppm_param_bcdef_periodic))) THEN

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
            ENDIF
        ENDDO ! end loop over remaining parts in ilist1
        !----------------------------------------------------------------------
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !----------------------------------------------------------------------
        IF (nlist2.NE.nlist1) THEN
           nlist1 = nlist2
           DO i=1,nlist1
              ilist1(i) = ilist2(i)
           ENDDO
        ENDIF

        !----------------------------------------------------------------------!
        !  Exit if the list is empty
        !----------------------------------------------------------------------!
        IF (nlist1.EQ.0) EXIT
        p => Mesh%subpatch%next()
        ipatch = ipatch + 1

      ENDDO !subpatch

      !-------------------------------------------------------------------------!
      !  Check whether we sold all the particles
      !-------------------------------------------------------------------------!
      IF (ppm_debug.GT.1) THEN
         IF (nlist2.GT.0) THEN
            stdout("Some particles seem to be well outside of all subpatches",&
                    " They will not take part in the p2m interpolation")
         ENDIF
      ENDIF

      max_partnumber = 0
      DO ipatch=1,nsubpatch
         IF(store_info(ipatch).GE.max_partnumber) THEN
            max_partnumber = store_info(ipatch)
         END IF
      ENDDO

      !----------------------------------------------------------------------
      !  Allocate particle list
      !----------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubpatch
      ldu(2) = max_partnumber

      CALL ppm_alloc(list_sub,ldu,iopt,info)
      or_fail_alloc("list_sub")

      list_sub=0

      !-------------------------------------------------------------------------!
      !  Initialize the particle list
      !-------------------------------------------------------------------------!
      DO ipart=1,this%Npart
         ilist1(ipart) = ipart
      ENDDO
      nlist1 = this%Npart

      !----------------------------------------------------------------------
      !  Loop over the subpatches (since the first domains are most likely
      !  to be empty, we look backwards to reduce the number of elements in
      !  nlist2 as fast as possible)
      !----------------------------------------------------------------------
      p => Mesh%subpatch%begin()
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
            IF (ppm_dim.EQ.3) THEN
               IF ((xp(1,ipart).GE.p%start_ext(1) .AND. &
               &    xp(2,ipart).GE.p%start_ext(2) .AND. &
               &    xp(3,ipart).GE.p%start_ext(3) .AND. &
               &    xp(1,ipart).LT.p%end_ext(1) .AND. &
               &    xp(2,ipart).LT.p%end_ext(2) .AND. &
               &    xp(3,ipart).LT.p%end_ext(3) ) ) THEN

                      IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                      &   (p%bc(2).GE.0   .AND.    &
                      &    p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                      &   (xp(2,ipart).LT.p%end(2) .OR.  &
                      &   (p%bc(4).GE.0   .AND.    &
                      &    p%bc(4).NE. ppm_param_bcdef_periodic)).AND.&
                      &   (xp(3,ipart).LT.p%end(3) .OR.  &
                      &   (p%bc(6).GE.0   .AND.    &
                      &    p%bc(6).NE. ppm_param_bcdef_periodic))   ) THEN

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
            ELSEIF (ppm_dim.EQ.2) THEN
               IF ((xp(1,ipart).GE.p%start_ext(1) .AND. &
               &    xp(2,ipart).GE.p%start_ext(2) .AND. &
               &    xp(1,ipart).LT.p%end_ext(1) .AND. &
               &    xp(2,ipart).LT.p%end_ext(2) ) ) THEN
                  IF ((xp(1,ipart).LT.p%end(1) .OR.  &
                  &   (p%bc(2).GE.0   .AND.    &
                  &    p%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
                  &   (xp(2,ipart).LT.p%end(2) .OR.  &
                  &   (p%bc(4).GE.0   .AND.    &
                  &    p%bc(4).NE. ppm_param_bcdef_periodic))) THEN
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
            ENDIF
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

         p => Mesh%subpatch%next()
         ipatch = ipatch + 1
      ENDDO !supatch

      !----------------------------------------------------------------------
      !  Check if we sold all the particles
      !----------------------------------------------------------------------
      IF (ppm_debug.GT.1) THEN
         IF (nlist2.GT.0) THEN
            stdout("Some particles seem to be outside from all subpatches",&
                    " They will not take part in the m2p interpolation")
         ENDIF
      ENDIF


      !-------------------------------------------------------------------------!
      !  Allocate and alias the weights if we need them.
      !-------------------------------------------------------------------------!
      max_partnumber = 0
      DO ipatch = 1,nsubpatch
         IF (store_info(ipatch).GE.max_partnumber) THEN
            max_partnumber = store_info(ipatch)
         ENDIF
      ENDDO

      9998 CONTINUE
      !-------------------------------------------------------------------------
      !  Create the discretization on the Mesh if it doesnt exist yet
      !-------------------------------------------------------------------------
      IF (.NOT.Field%is_discretized_on(Mesh)) THEN
          CALL Field%discretize_on(Mesh,info)
          or_fail("Could not discretize Field on Mesh")
      ENDIF

      CALL Mesh%zero(Field,info)
      or_fail("this%zero(): Failed to zero the field on this Mesh")

      IF (this%Npart.EQ.0) GOTO 9997

      !Get a pointer to the data on the particles
      IF (Field%lda.EQ.1) THEN
         NULLIFY(up_1d)
         CALL this%get_field(Field,up_1d,info,read_only=.TRUE.)
      ELSE
         NULLIFY(up_2d)
         CALL this%get_field(Field,up_2d,info,read_only=.TRUE.)
      ENDIF
      or_fail("Could not get pointer to discretized data on the particles")

      !For checking only - display warnings, but should in principle fail
      !gracefully
      NULLIFY(prop)
      CALL Field%get_discr(this,prop,info)
      or_fail("Could not get pointer to discretization object on the particles")

      IF (.NOT. prop%flags(ppm_ppt_partial)) THEN
         stdout("WARNING: property has not been partial-mapped. This can produce garbage data")
      ENDIF
      !

      SELECT CASE(kernel)
      CASE (ppm_param_rmsh_kernel_mp4)
         IF (Field%lda.EQ.1) THEN
            IF (ppm_dim.EQ.2) THEN
               NULLIFY(dummy_2d)
               CALL p2m_interp_mp4(Mesh,Field,dummy_2d,xp,up_1d,info)
            ELSE
               NULLIFY(dummy_3d)
               CALL p2m_interp_mp4(Mesh,Field,dummy_3d,xp,up_1d,info)
            ENDIF
         ELSE
            IF (ppm_dim.EQ.2) THEN
               NULLIFY(dummy_3d)
               CALL p2m_interp_mp4(Mesh,Field,dummy_3d,Field%lda,xp,up_2d,info)
            ELSE
               NULLIFY(dummy_4d)
               CALL p2m_interp_mp4(Mesh,Field,dummy_4d,Field%lda,xp,up_2d,info)
            ENDIF
         ENDIF
         or_fail("p2m_interp_mp4 failed")

      CASE (ppm_param_rmsh_kernel_bsp2)
          !IF (Field%lda.EQ.1) THEN
              !IF (ppm_dim.EQ.2) THEN
                 !NULLIFY(dummy_2d)
                 !CALL p2m_interp_bsp2(Mesh,Field,dummy_2d,xp,up_1d,info)
              !ELSE
                 !NULLIFY(dummy_3d)
                 !CALL p2m_interp_bsp2(Mesh,Field,dummy_3d,xp,up_1d,info)
              !ENDIF
          !ELSE
              !IF (ppm_dim.EQ.2) THEN
                 !NULLIFY(dummy_3d)
                 !CALL p2m_interp_bsp2(Mesh,Field,dummy_3d,Field%lda,xp,up_2d,info)
              !ELSE
                 !NULLIFY(dummy_4d)
                 !CALL p2m_interp_bsp2(Mesh,Field,dummy_4d,Field%lda,xp,up_2d,info)
              !ENDIF
          !ENDIF
          fail("BSp2 has not been implemented yet - but thats quick to do...")

      CASE DEFAULT
         info = ppm_error_error
         fail("This scheme is not available. Use ppm_rmsh_remesh for other kernels")

      END SELECT         ! kernel type

      9997 CONTINUE
      !-------------------------------------------------------------------------!
      !  Now map the ghosts in order to get consistent values at the border of
      !  the subdomains.
      !-------------------------------------------------------------------------!
      CALL Mesh%map_ghost_put(info)
      or_fail("map_ghost_put")

      CALL Field%map_ghost_push(Mesh,info)
      or_fail("map_ghost_push")
      !non-blocking send
      CALL Mesh%map_isend(info)
      or_fail("map_isend")
      CALL Field%map_ghost_pop(Mesh,info)
      or_fail("map_ghost_pop")

      IF (PRESENT(p2m_bcdef)) THEN
         IF (Field%lda.EQ.1) THEN
            IF (ppm_dim.EQ.2) THEN
               CALL p2m_interp_bc(Mesh,Field,dummy_2d,p2m_bcdef,info)
            ELSE
               CALL p2m_interp_bc(Mesh,Field,dummy_3d,p2m_bcdef,info)
            ENDIF
         ELSE
            IF (ppm_dim.EQ.2) THEN
               CALL p2m_interp_bc(Mesh,Field,dummy_3d,Field%lda,p2m_bcdef,info)
            ELSE
               CALL p2m_interp_bc(Mesh,Field,dummy_4d,Field%lda,p2m_bcdef,info)
            ENDIF
         ENDIF
         or_fail("p2m_interp_bc failed")
      END IF

      IF (this%Npart.EQ.0) GOTO 9999
      !-------------------------------------------------------------------------!
      ! Deallocation of the arrays....
      !-------------------------------------------------------------------------!
      DEALLOCATE(ilist1,ilist2,STAT=info)
      or_fail_dealloc("ilist1 & ilist2")

      iopt = ppm_param_dealloc
      ldu(1) = 0

      CALL ppm_alloc(store_info,ldu,iopt,info)
      or_fail_dealloc("store_info")

      CALL ppm_alloc(list_sub,ldu,iopt,info)
      or_fail_dealloc("list_sub")

      !Updating internal state variables
      CALL this%set_xp(xp,info,read_only=.TRUE.)
      or_fail("Set_xp failed")

      end_subroutine()
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            fail("Please call ppm_init first!",&
            ppm_err_ppm_noinit,exit_point=8888)
        ENDIF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
            fail("Wrong kernel definition",&
                ppm_err_ppm_noinit,exit_point=8888)
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF (.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
        &   .OR.(kernel_support.EQ.6))) THEN
            fail("Wrong kernel support",ppm_err_argument,exit_point=8888)
        END IF
        IF (this%Npart .GT. 0) THEN
           IF (SIZE(this%xp,2) .LT. this%Npart) THEN
            fail("not enough particles contained in xp",exit_point=8888)
           ENDIF
           IF (SIZE(this%xp,1) .LT.dim) THEN
            fail("leading dimension of xp insufficient",exit_point=8888)
           ENDIF
        ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE DTYPE(part_p2m)

