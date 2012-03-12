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
     &           ghostsize,field_up,info,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_ds_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,device)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_sv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_dv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,device)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_ss_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_ds_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,device)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_sv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_dv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,device)
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
      USE ppm_module_map
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
      REAL(MK) , DIMENSION(:)                     , POINTER :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:    ), INTENT(IN),  POINTER :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:  ), INTENT(IN),  POINTER :: field_up
#endif
      !!! field from which to interpolate
#elif __MODE == __VEC
      INTEGER                         , INTENT(IN   )       :: lda
      !!! leading dimension of up
      REAL(MK) , DIMENSION(:,:)       ,             POINTER :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , INTENT(IN), POINTER :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , INTENT(IN), POINTER :: field_up
#endif
      !!! field from which to interpolate
#endif
      REAL(MK), DIMENSION(:,:)        , INTENT(IN), POINTER :: xp
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
      INTEGER                        , OPTIONAL      :: device
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(MK),  DIMENSION(ppm_dim)          :: dxi,dx
      REAL(MK),  DIMENSION(ppm_dim)          :: len_phys
      REAL(MK)                               :: x1,x2,x3,epsilon
      INTEGER                                :: kernel_support
      INTEGER,  DIMENSION(ppm_dim+2)         :: ldu,ldl
      INTEGER,  DIMENSION(ppm_dim)           :: Nc
      INTEGER                                :: i,j,k,ii,jj,kk,iidec,maptype
      INTEGER                                :: jjdec,nb_sub,npart,ipart
      INTEGER                                :: kkdec,ip1,nbpt_z,nlist1
      INTEGER                                :: ip2,ip3,nbpt_x,nbpt_y,iface
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
      TYPE(ppm_t_equi_mesh), POINTER         :: p_mesh => NULL()
      TYPE(ppm_t_topo)     , POINTER         :: topo   => NULL()

      ! Some arrays and variables defined as workaround tools for now
      REAL(mk), DIMENSION(:),        POINTER :: p_coor => NULL()
      REAL(mk), DIMENSION(:),        POINTER :: p_mass => NULL()
      REAL(mk), DIMENSION(:),        POINTER :: mesh_mass => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: minphys_with_ghost => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: maxphys_with_ghost => NULL()
      INTEGER,  DIMENSION(:),        POINTER :: mesh_size => NULL()

      INTEGER,  DIMENSION(ppm_dim)            :: mesh_coor

      INTEGER                                :: mass
      INTEGER                                :: nmesh
      INTEGER                                :: ndim_in
      INTEGER                                :: nmass_in
      INTEGER                                :: interp_kernel
      INTEGER                                :: dimn 

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_interp_m2p',t0,info)

      IF(PRESENT(device))  THEN
          IF(device .EQ. ppm_param_device_gpu)  THEN
#ifdef __GPU
! Calling GPU subroutine
#if   __MODE == __SCA
      call gpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __MODE == __VEC
      call gpu_m2p(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#else
! ----- if user chooses GPU even though GPU is NOT available ----- !
      CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',    &
     &  'PPM library was not compiled with GPU support.   &
	 &   Switching to CPU version.', &
     &           __LINE__,info)
#if   __MODE == __SCA
      call cpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __MODE == __VEC
      call cpu_m2p(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif 
          ELSE !User explicitly chooses to run CPU version
#if   __MODE == __SCA
      call cpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __MODE == __VEC
      call cpu_m2p(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
          ENDIF
      ELSE ! if user does not choose the device explicitly...
     !-------------------------------------------------------------------------!
     !  Call the p2m interpolation subroutine, depending on availability of GPU
     !-------------------------------------------------------------------------!
! ----- if GPU is available ----- !
#ifdef __GPU
! Calling GPU subroutine
#if   __MODE == __SCA
      call gpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __MODE == __VEC
      call gpu_m2p(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#else
! ----- if GPU is NOT available ----- !
#if   __MODE == __SCA
      call cpu_m2p(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info)
#elif __MODE == __VEC
      call cpu_m2p(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info)
#endif
#endif 
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_interp_m2p',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
        IF (.NOT. ppm_initialized) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_ppm_noinit,'ppm_interp_m2p',  &
     &               'Please call ppm_init first!',__LINE__,info)
            GOTO 8888
        ENDIF
        IF (Np .GT. 0) THEN
           IF (SIZE(xp,2) .LT. Np) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                     'not enough particles contained in xp',__LINE__,info)
            GOTO 8888
           ENDIF
           IF (SIZE(xp,1) .LT.dim) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                     'leading dimension of xp insufficient',__LINE__,info)
            GOTO 8888
           ENDIF
        ENDIF
        IF (Np .LE. 0) THEN
           IF (Np .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                     'particles must be specified',__LINE__,info)
            GOTO 8888
           END IF
           GOTO 8888
        END IF
        IF ((kernel.LT.1).OR.(kernel.GT.4)) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                     'wrong kernel definition',__LINE__,info)
           GOTO 8888
        END IF
        kernel_support = ppm_rmsh_kernelsize(kernel)*2
        IF(.NOT.((kernel_support.EQ.2).OR.(kernel_support.EQ.4) &
     &               .OR.(kernel_support.EQ.6))) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                     'wrong kernel support',__LINE__,info)
           GOTO 8888
        END IF
        CALL ppm_check_topoid(topoid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                 'topo_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
        CALL ppm_check_meshid(topoid,meshid,lok,info)
        IF (.NOT.lok) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_argument,'ppm_interp_m2p',  &
     &                 'mesh_id is invalid!',__LINE__,info)
           GOTO 8888
        ENDIF
 8888   CONTINUE
      END SUBROUTINE check

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



