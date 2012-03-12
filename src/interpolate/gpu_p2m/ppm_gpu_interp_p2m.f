      !------------------------------------------------------------------------!
      !     Subroutine   :                 ppm_gpu_interp_p2m
      !------------------------------------------------------------------------!
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
      !------------------------------------------------------------------------!


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE gpu_p2m_ss_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE gpu_p2m_ds_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE gpu_p2m_sv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE gpu_p2m_dv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE gpu_p2m_ss_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE gpu_p2m_ds_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE gpu_p2m_sv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE gpu_p2m_dv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#endif
#endif
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
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_write
      USE ppm_module_map
      USE ppm_module_check_id
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
     !-------------------------------------------------------------------------!
     ! Arguments
     !-------------------------------------------------------------------------!
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , INTENT(IN)   :: up
      !!! particle weights from which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#elif __MODE == __VEC
      INTEGER                         , INTENT(in)   :: lda
      !!! leading dimension of up
      REAL(MK) , DIMENSION(:,:)       , INTENT(IN)   :: up
      !!! particle weights from which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#endif
      REAL(MK), DIMENSION(:,:)       , INTENT(IN)    :: xp
      !!! particle positions
      INTEGER , DIMENSION(:  )       , INTENT(IN)    :: ghostsize
      !!! The size (width) of the ghost layer
      INTEGER                        , INTENT(IN   ) :: Np
      !!! number of particles
      INTEGER                        , INTENT(IN   ) :: meshid
      !!! id of the mesh (user)
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! topology identifier (user numbering) of target
      INTEGER                        , INTENT(IN   ) :: kernel
      !!! Choice of the kernel used to compute the weights.
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
#if   __MODE == __SCA
     INTEGER, DIMENSION(:  )        , POINTER, OPTIONAL :: p2m_bcdef
#elif __MODE == __VEC
     INTEGER, DIMENSION(:,:)        , POINTER, OPTIONAL :: p2m_bcdef
#endif
     !!! Boundary conditions used for this interpolation routine, they may be
     !!! different than the boundary conditions set in the topology.
     !!! The values in this array can be one of:
     !!!
     !!! - ppm_param_bcdef_symmetry
     !!! - ppm_param_bcdef_antisymmetry
     !-------------------------------------------------------------------------!
     ! Local variables
     !-------------------------------------------------------------------------!
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart   => NULL()
      INTEGER,  DIMENSION(:,:)     , POINTER :: ndata    => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(mk), DIMENSION(ppm_dim)           :: dxi,dx
      REAL(mk)                               :: dxx,dxy,dxz,dxxi,dxyi,dxzi
      REAL(mk), DIMENSION(ppm_dim)           :: len_phys
      REAL(mk)                               :: x1,x2,x3
      INTEGER                                :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
      INTEGER,  DIMENSION(ppm_dim)           :: Nc
      INTEGER                                :: i,j,k,l,ii,jj,kk,iidec,maptype
      INTEGER                                :: xlo,xhi,ylo,yhi,zlo,zhi
      INTEGER                                :: jjdec,nb_sub,npart,ipart
      INTEGER                                :: kkdec,ip1,nlist1
      INTEGER                                :: ip2,ip3,iface
      INTEGER                                :: isub,ifrom,ito,ip,dim,iopt,isubl
      INTEGER                                :: max_partnumber,idom,nlist2,idoml
      INTEGER, DIMENSION(ppm_dim)            :: Nm
      INTEGER                                :: nsubs
      INTEGER, DIMENSION(6)                  :: bcdef
      INTEGER                                :: iq
      LOGICAL                                :: internal_weights,lok
      ! aliases
      REAL(mk), DIMENSION(:,:),      POINTER :: min_sub => NULL()
      REAL(mk), DIMENSION(:,:),      POINTER :: max_sub => NULL()
      REAL(mk)                               :: myeps
      REAL(mk)                               :: tim1s, tim1e
      REAL(mk)                               :: xp1,xp2,xp3
      REAL(mk)                               :: wx1,wx2,wx3
      REAL(mk), DIMENSION(ppm_dim)           :: x0
      REAL(mk)                               :: x01,x02,x03
      INTEGER                                :: ldn
      CHARACTER(len=256)                     :: msg
      TYPE(ppm_t_equi_mesh), POINTER         :: p_mesh => NULL()
      TYPE(ppm_t_topo)     , POINTER         :: topo   => NULL()

     !-------------------------------------------------------------------------!
     ! Local variables: Temporary arrays and variables defined for the GPU code 
     !-------------------------------------------------------------------------!
      REAL(mk), DIMENSION(:),        POINTER :: p_coor => NULL()
      ! 1D array keeping particle arrays to be passed to GPU
      REAL(mk), DIMENSION(:),        POINTER :: p_prop => NULL()
      ! 1D array keeping particle properties to be passed to GPU
      REAL(mk), DIMENSION(:),        POINTER :: mesh_prop => NULL()
      ! 1D array keeping particle properties to be passed to GPU
      INTEGER,  DIMENSION(ppm_dim)           :: mesh_size
      ! Mesh size of the given subdomain
      INTEGER,  DIMENSION(ppm_dim)           :: mesh_coor
      ! mesh coordinates with the offset by istart
      INTEGER                                :: nmesh
      ! Number of mesh points in the current subdomain
      INTEGER                                :: nprop
      ! Number of properties
      INTEGER                                :: interp_kernel
      ! Interpolation kernel
      INTEGER                                :: dimn 
      ! Counter for the do-loop iterating over dimension
      INTEGER                                :: prop
      ! Counter for the do-loop iterating over properties
      REAL(mk), DIMENSION(:,:),      POINTER :: min_sub_g => NULL()
      ! Minimum physical extents of dubdomain including ghost layers
      REAL(mk), DIMENSION(:,:),      POINTER :: max_sub_g => NULL()
      ! Maximum physical extents of dubdomain including ghost layers

     !-------------------------------------------------------------------------!
     !  Initialise
     !-------------------------------------------------------------------------!

     !-------------------------------------------------------------------------!
     !  This is a hack! Somehow having trouble with constructors in the
     !  ppm_module_data_rmsh module
     !-------------------------------------------------------------------------!
      ppm_rmsh_kernelsize = (/1,2,2,4/)

      CALL substart('ppm_interp_p2m',t0,info)

      dim = ppm_dim
      internal_weights = .FALSE.

     !-------------------------------------------------------------------------!
     !  Check arguments
     !-------------------------------------------------------------------------!
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Get the meshid
      !-------------------------------------------------------------------------
      p_mesh => topo%mesh(meshid)
      !-------------------------------------------------------------------------
      !  Get istart
      istart => p_mesh%istart

     !-------------------------------------------------------------------------!
     !  Assignment of the useful arrays/scalar
     !-------------------------------------------------------------------------!
      Nm(1:dim) = p_mesh%Nm(1:dim)
      bcdef(1:(2*dim)) = topo%bcdef(1:(2*dim))
      nsubs = topo%nsublist
      !-------------------------------------------------------------------------!
      !  If there is nothing to do, do nearly nothing
      !-------------------------------------------------------------------------!
      IF(Np.EQ.0) GOTO 9998

      !-------------------------------------------------------------------------!
      !  Alloc memory for particle lists
      !  The awesome ppm_alloc will (re)allocate them, so we dont need an init
      !  or a finalize routine.
      !-------------------------------------------------------------------------!
      iopt   = ppm_param_alloc_fit
      ldu(1) = Np
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                    'particle list 1 ILIST1',__LINE__,info)
        GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info .NE. 0) THEN
       info = ppm_error_fatal
       CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                    'particle list 2 ILIST2',__LINE__,info)
       GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = nsubs
      CALL ppm_alloc(store_info,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'store_info allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF

     !-------------------------------------------------------------------------!
     !  Mesh spacing
     !-------------------------------------------------------------------------!
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
      ENDDO

     !-------------------------------------------------------------------------!
     !  Initialize the particle list
     !-------------------------------------------------------------------------!
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

     !-------------------------------------------------------------------------!
     !  ALLOCATE min_sub_g and max_sub_g which are the subdomain extents       !
     !  including ghost layers.                                                !
     !-------------------------------------------------------------------------!
      iopt   = ppm_param_alloc_fit
      ldu(1) = size(min_sub,1)
      ldu(2) = size(min_sub,2)
      CALL ppm_alloc(min_sub_g,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'min_sub_h allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF

      ! sizes will be the same as min_sub_g
      CALL ppm_alloc(max_sub_g,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'max_sub_h allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF
     !-------------------------------------------------------------------------!
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !-------------------------------------------------------------------------!
      DO idom = topo%nsublist,1,-1
        idoml = topo%isublist(idom)
        DO dimn = 1, ppm_dim
            min_sub_g(dimn, idoml) = min_sub(dimn, idoml) - dx(dimn)*ghostsize(dimn)
            max_sub_g(dimn, idoml) = max_sub(dimn, idoml) + dx(dimn)*ghostsize(dimn)
        ENDDO
        !----------------------------------------------------------------------!
        !  Loop over the remaining particles
        !----------------------------------------------------------------------!
        nlist2 = 0
        npart = 0
        DO i=1,nlist1
           ipart = ilist1(i)

           !-------------------------------------------------------------------!
           !  If the particle is inside the current subdomain, assign it
           !-------------------------------------------------------------------!
#if __DIME == __3D
           !----------------------------------------------------------------!
           ! the particle is in the closure of the subdomain
           !----------------------------------------------------------------!
           IF( ( xp(1,ipart).GE.min_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).GE.min_sub_g(2,idoml) .AND. &
     &                xp(3,ipart).GE.min_sub_g(3,idoml) .AND. &
     &                xp(1,ipart).LE.max_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).LE.max_sub_g(2,idoml) .AND. &
     &                xp(3,ipart).LE.max_sub_g(3,idoml) ) ) THEN

              IF(   (xp(1,ipart).LT.max_sub_g(1,idoml) .OR.  &
     &                    (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                    bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(2,ipart).LT.max_sub_g(2,idoml) .OR.  &
     &                    (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                    bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(3,ipart).LT.max_sub_g(3,idoml) .OR.  &
     &                    (topo%subs_bc(6,idoml).EQ.1   .AND.    &
     &                    bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN

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

#elif __DIME == __2D
           IF( ( xp(1,ipart).GE.min_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).GE.min_sub_g(2,idoml) .AND. &
     &                xp(1,ipart).LE.max_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).LE.max_sub_g(2,idoml) ) ) THEN

              IF(   (xp(1,ipart).LT.max_sub_g(1,idoml) .OR.  &
     &                    (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                    bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(2,ipart).LT.max_sub_g(2,idoml) .OR.  &
     &                    (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                    bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

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
#endif
        END DO
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
     END DO

     !-------------------------------------------------------------------------!
     !  Check that we sold all the particles
     !-------------------------------------------------------------------------!
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m',  &
     &                    'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF

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
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'problem in internal allocation',__LINE__,info)
        GOTO 9999
     END IF

     list_sub=0

     !-------------------------------------------------------------------------!
     !  Initialize the particle list
     !-------------------------------------------------------------------------!
     nlist1     = 0
     DO ipart=1,Np
        nlist1         = nlist1 + 1
        ilist1(nlist1) = ipart
     ENDDO

     !-------------------------------------------------------------------------!
     !  Loop over the subdomains (since the first domains are most likely
     !  to be empty, we look backwards to reduce the number of elements in
     !  nlist2 as fast as possible)
     !-------------------------------------------------------------------------!
     DO idom = topo%nsublist,1,-1
        idoml = topo%isublist(idom)
        !----------------------------------------------------------------------!
        !  loop over the remaining particles
        !----------------------------------------------------------------------!
        nlist2 = 0
        npart = 0

        DO i=1,nlist1
           ipart = ilist1(i)
           !-------------------------------------------------------------------!
           !  If the particle is inside the current subdomain, assign it
           !-------------------------------------------------------------------!
#if __DIME == __3D
           IF( ( xp(1,ipart).GE.min_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).GE.min_sub_g(2,idoml) .AND. &
     &                xp(3,ipart).GE.min_sub_g(3,idoml) .AND. &
     &                xp(1,ipart).LE.max_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).LE.max_sub_g(2,idoml) .AND. &
     &                xp(3,ipart).LE.max_sub_g(3,idoml) ) ) THEN

              IF(   (xp(1,ipart).LT.max_sub_g(1,idoml) .OR.  &
     &                    (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                    bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(2,ipart).LT.max_sub_g(2,idoml) .OR.  &
     &                    (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                    bcdef(4).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(3,ipart).LT.max_sub_g(3,idoml) .OR.  &
     &                    (topo%subs_bc(6,idoml).EQ.1   .AND.    &
     &                    bcdef(6).NE. ppm_param_bcdef_periodic))   ) THEN


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
#elif __DIME == __2D
           IF( ( xp(1,ipart).GE.min_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).GE.min_sub_g(2,idoml) .AND. &
     &                xp(1,ipart).LE.max_sub_g(1,idoml) .AND. &
     &                xp(2,ipart).LE.max_sub_g(2,idoml) ) ) THEN

              IF(   (xp(1,ipart).LT.max_sub_g(1,idoml) .OR.  &
     &                    (topo%subs_bc(2,idoml).EQ.1   .AND.    &
     &                    bcdef(2).NE. ppm_param_bcdef_periodic)).AND.&
     &                    (xp(2,ipart).LT.max_sub_g(2,idoml) .OR.  &
     &                    (topo%subs_bc(4,idoml).EQ.1   .AND.    &
     &                    bcdef(4).NE. ppm_param_bcdef_periodic))) THEN

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
#endif
        END DO
        !----------------------------------------------------------------------!
        !  Copy the lists (well, only if nlist2 changed - decreased)
        !----------------------------------------------------------------------!
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
     END DO

     !-------------------------------------------------------------------------!
     !  Check that we sold all the particles
     !-------------------------------------------------------------------------!
     IF (nlist2.GT.0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_part_unass,'ppm_interp_p2m_3d',  &
     &                    'MAJOR PROBLEM',__LINE__,info)
        GOTO 9999
     ENDIF

     !-------------------------------------------------------------------------!
     !  Allocate and alias the weights if we need them.
     !-------------------------------------------------------------------------!
     max_partnumber = 0
     DO idom = 1,topo%nsublist
        IF(store_info(idom).GE.max_partnumber) THEN
           max_partnumber = store_info(idom)
        END IF
     END DO

9998 CONTINUE
      IF(np.EQ.0) GOTO 9997

      SELECT CASE(kernel)

      CASE(ppm_param_rmsh_kernel_mp4)
      CASE(ppm_param_rmsh_kernel_bsp2)
      CASE DEFAULT
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',    &
     &                    'This scheme is currently not available. &
	 &                    Use ppm_rmsh_remesh for other kernels.', &
     &                    __LINE__,info)
     END SELECT         ! kernel type

     !-------------------------------------------------------------------------!
     ! Here, we run the interpolation over subdomains
     !-------------------------------------------------------------------------!
      ! Set number of properties
#if   __MODE == __SCA
      nprop = 1 
#elif __MODE == __VEC
      nprop = lda
#endif

      ! Allocate arrays defining physical domain extent and etc.
      iopt   = ppm_param_alloc_fit
      ldu(1) = max_partnumber*ppm_dim
      CALL ppm_alloc(p_coor,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'p_coor allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF

      iopt   = ppm_param_alloc_fit
      ldu(1) = max_partnumber*nprop
      CALL ppm_alloc(p_prop,ldu,iopt,info)
      IF (info .NE. 0) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'p_prop allocation : problem',__LINE__,info)
        GOTO 9999
      ENDIF

      ! set interp_kernel parameter
      IF(kernel .EQ. ppm_param_rmsh_kernel_mp4) THEN
          interp_kernel = -1;
      ELSEIF(kernel .EQ. ppm_param_rmsh_kernel_bsp2) THEN
          interp_kernel = -2;
      ENDIF

      DO idom = topo%nsublist,1,-1
          isub = topo%isublist(idom)
 
          IF(store_info(isub) .EQ. 0) THEN
              CYCLE ! Skip this iteration, move onto next isub
          ENDIF

#if __DIME == __2D
          mesh_size(1) = p_mesh%nnodes(1,isub)
          mesh_size(2) = p_mesh%nnodes(2,isub)
          nmesh        = nprop*mesh_size(1)*mesh_size(2)
#elif __DIME == __3D
          mesh_size(1) = p_mesh%nnodes(1,isub)
          mesh_size(2) = p_mesh%nnodes(2,isub)
          mesh_size(3) = p_mesh%nnodes(3,isub)
          nmesh        = nprop*mesh_size(1)*mesh_size(2)*mesh_size(3)
#endif

          iopt   = ppm_param_alloc_fit
          ldu(1) = nmesh
          CALL ppm_alloc(mesh_prop,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_interp_p2m_3d',     &
     &                     'mesh_prop allocation : problem',__LINE__,info)
              GOTO 9999
          ENDIF

          DO ip  = 1, store_info(isub)
              iq = list_sub(isub, ip)
              DO dimn = 1, ppm_dim
                  p_coor(ppm_dim*(ip - 1) + dimn) = xp(dimn, iq)
              ENDDO
#if   __MODE == __SCA
              p_prop(ip) = up(iq)
#elif __MODE == __VEC
              DO prop = 1,nprop
                  p_prop(nprop*(ip-1) + prop) = up(prop, iq)
              ENDDO
#endif
          ENDDO

#if __KIND   == __SINGLE_PRECISION
          CALL ppm_gpu_p2m_interpolation_s(p_coor, p_prop, mesh_prop,  &
     &              min_sub(:,isub), max_sub(:,isub), mesh_size,     &
     &              store_info(idom), ppm_dim, nprop, interp_kernel)
#elif __KIND == __DOUBLE_PRECISION
          CALL ppm_gpu_p2m_interpolation_d(p_coor, p_prop, mesh_prop,  &
     &              min_sub(:,isub), max_sub(:,isub), mesh_size,     &
     &              store_info(idom), ppm_dim, nprop, interp_kernel)
#endif

#if __DIME == __2D
          DO j = 0,mesh_size(2)-1
              DO i = 0,mesh_size(1)-1
                  mesh_coor(1) = i + istart(1, isub) 
                  mesh_coor(2) = j + istart(2, isub) 
#if __MODE == __SCA
                  field_up(mesh_coor(1),mesh_coor(2),1) =           &
     &                     mesh_prop(i + j*mesh_size(1) + 1)
#elif __MODE == __VEC
                  DO prop = 1, nprop
                      field_up(prop,mesh_coor(1),mesh_coor(2),1) = &
     &                   mesh_prop(nprop*(i + j*mesh_size(1)) + prop)
                  ENDDO
#endif
              ENDDO
          ENDDO
#elif __DIME == __3D
          DO k = 0,mesh_size(3)-1
              DO j = 0,mesh_size(2)-1
                  DO i = 0,mesh_size(1)-1
                      mesh_coor(1) = i + istart(1, isub) 
                      mesh_coor(2) = j + istart(2, isub) 
                      mesh_coor(3) = k + istart(3, isub) 
#if __MODE == __SCA
                      field_up(mesh_coor(1),mesh_coor(2),mesh_coor(3),1)&
     & = mesh_prop(i + j*mesh_size(1) + k*mesh_size(1)*mesh_size(2) + 1) 
#elif __MODE == __VEC
                      DO prop = 1, nprop
                         field_up(prop,mesh_coor(1),mesh_coor(2),mesh_coor(3),1) &
      &                    = mesh_prop(nprop*(i + j*mesh_size(1) +               &
      &                   k*mesh_size(1)*mesh_size(2)) + prop) 
                      ENDDO
#endif
                  ENDDO
              ENDDO
          ENDDO
#endif

          iopt = ppm_param_dealloc
          ldu(1) = 0
          CALL ppm_alloc(mesh_prop,ldu,iopt,info)
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'mesh_prop deallocation',__LINE__,info)
              GOTO 9999
          END IF
      ENDDO ! end of looping over subs

9997 CONTINUE
     IF(np.EQ.0) GOTO 9999
     !-------------------------------------------------------------------------!
     ! Deallocation of the arrays....
     !-------------------------------------------------------------------------!
      iopt = ppm_param_dealloc
      ldu(1) = 0
      CALL ppm_alloc(p_coor,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'p_coor deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(p_prop,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'p_prop deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(min_sub_g,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'min_sub_g deallocation',__LINE__,info)
         GOTO 9999
      END IF
      CALL ppm_alloc(max_sub_g,ldu,iopt,info)
      IF (info.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'max_sub_g deallocation',__LINE__,info)
         GOTO 9999
      END IF
     CALL ppm_alloc(ilist1,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m', &
     &              'pb in ilist1 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(ilist2,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(store_info,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF
     CALL ppm_alloc(list_sub,ldu,iopt,info)
     IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc, &
     &              'ppm_interp_p2m',  &
     &              'pb in ilist2 deallocation',__LINE__,info)
        GOTO 9999
     END IF

     !-------------------------------------------------------------------------!
     !  Return
     !-------------------------------------------------------------------------!
 9999 CONTINUE
      CALL substop('ppm_interp_p2m',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_interp_p2m',  &
     &               'Please call ppm_init first!',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'not enough particles contained in xp',__LINE__,info)
            GOTO 8888
           ENDIF
           IF (SIZE(xp,1) .LT.dim) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'leading dimension of xp insufficient',__LINE__,info)
            GOTO 8888
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'particles must be specified',__LINE__,info)
            GOTO 8888
           END IF
           GOTO 8888
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'wrong kernel definition',__LINE__,info)
           GOTO 8888
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
     &               .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                     'wrong kernel support',__LINE__,info)
           GOTO 8888
        END IF
        CALL ppm_check_topoid(topoid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                 'topo_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
        CALL ppm_check_meshid(topoid,meshid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',  &
     &                 'mesh_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE gpu_p2m_ss_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE gpu_p2m_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE gpu_p2m_sv_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE gpu_p2m_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE gpu_p2m_ss_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE gpu_p2m_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE gpu_p2m_sv_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE gpu_p2m_dv_3d
#endif
#endif
#endif

