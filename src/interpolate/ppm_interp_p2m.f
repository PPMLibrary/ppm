      !------------------------------------------------------------------------!
      !     Subroutine   :                 ppm_interp_p2m
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
      SUBROUTINE p2m_ss_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_ds_2d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_sv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_dv_2d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_ss_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_ds_3d(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_sv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_dv_3d(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef,device)
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
      REAL(MK) , DIMENSION(:), INTENT(IN), POINTER   :: up
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
      REAL(MK) , DIMENSION(:,:) ,INTENT(IN), POINTER :: up
      !!! particle weights from which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:,:  ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:,:) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#endif
      REAL(MK), DIMENSION(:,:), INTENT(IN), POINTER  :: xp
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
     INTEGER, DIMENSION(:),   POINTER  , OPTIONAL      :: p2m_bcdef
#elif __MODE == __VEC
     INTEGER, DIMENSION(:,:), POINTER  , OPTIONAL      :: p2m_bcdef
#endif
      !!  Device choice, optional
     INTEGER                         , OPTIONAL      :: device
     !!! Boundary conditions used for this interpolation routine, they may be
     !!! different than the boundary conditions set in the topology.
     !!! The values in this array can be one of:
     !!!
     !!! - ppm_param_bcdef_symmetry
     !!! - ppm_param_bcdef_antisymmetry
     !-------------------------------------------------------------------------!
     ! Local variables
     !-------------------------------------------------------------------------!
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(mk), DIMENSION(ppm_dim)           :: len_phys
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

      INTEGER                                :: i
      INTEGER                                :: Mp
      REAL(MK), DIMENSION(:),        POINTER :: ghostsize_f
      REAL(MK), DIMENSION(:),        POINTER :: dxi,dx
      INTEGER,  DIMENSION(ppm_dim)           :: Nc
      INTEGER, DIMENSION(ppm_dim)            :: Nm

     !-------------------------------------------------------------------------!
     !  Initialise
     !-------------------------------------------------------------------------!
      CALL substart('ppm_interp_p2m',t0,info)

      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Get the meshid
      !-------------------------------------------------------------------------
      p_mesh => topo%mesh(meshid)
      !-------------------------------------------------------------------------
      !  Get istart
      !-------------------------------------------------------------------------

#if   __MODE == __SCA
     ldn = 1
#elif __MODE == __VEC
     ldn = lda
#endif

     ALLOCATE(ghostsize_f(ppm_dim), dx(ppm_dim), dxi(ppm_dim))
     Nm(1:ppm_dim) = ppm_topo(topoid)%t%mesh(meshid)%Nm(1:ppm_dim)

#if   __KIND == __SINGLE_PRECISION
     min_phys => ppm_topo(topoid)%t%min_physs
     max_phys => ppm_topo(topoid)%t%max_physs
#elif __KIND == __DOUBLE_PRECISION
     min_phys => ppm_topo(topoid)%t%min_physd
     max_phys => ppm_topo(topoid)%t%max_physd
#endif

      DO i = 1,ppm_dim
        Nc(i)          = Nm(i) - 1
        len_phys(i)    = max_phys(i) - min_phys(i)
        dx(i)          = len_phys(i)/REAL(Nc(i),MK)
        ghostsize_f(i) = ghostsize(i)*dx(i)
      ENDDO

      IF(PRESENT(device))  THEN
          IF(device .EQ. ppm_param_device_gpu)  THEN
#ifdef __GPU
      CALL ppm_map_part_ghost_get(topoid,xp,ppm_dim,Np,0,ghostsize_f(1),info)
#if   __MODE == __SCA
      CALL ppm_map_part_push(up,Np,info)
#elif __MODE == __VEC
      CALL ppm_map_part_push(up,ldn,Np,info)
#endif
      CALL ppm_map_part_send(Np,Mp,info)
#if   __MODE == __SCA
      CALL ppm_map_part_pop(up,Np,Mp,info)
#elif __MODE == __VEC
      CALL ppm_map_part_pop(up,ldn,Np,Mp,info)
#endif
      CALL ppm_map_part_pop(xp,ppm_dim,Np,Mp,info)

#if   __MODE == __SCA
      CALL gpu_p2m(topoid,meshid,xp,Mp,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __MODE == __VEC
      CALL gpu_p2m(topoid,meshid,xp,Mp,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#else
! ----- if user chooses GPU even though GPU is NOT available ----- !
      CALL ppm_error(ppm_err_argument,'ppm_interp_p2m',    &
     &  'PPM library was not compiled with GPU support.   &
	 &   Switching to CPU version.', &
     &           __LINE__,info)
#if   __MODE == __SCA
      CALL cpu_p2m(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif   __MODE == __VEC
      CALL cpu_p2m(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#endif 
          ELSE !User explicitly chooses to run CPU version
#if   __MODE == __SCA
      CALL cpu_p2m(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif   __MODE == __VEC
      CALL cpu_p2m(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
          ENDIF
      ELSE ! if user does not choose the device explicitly...
     !-------------------------------------------------------------------------!
     !  Call the p2m interpolation subroutine, depending on availability of GPU
     !-------------------------------------------------------------------------!
! ----- if GPU is available ----- !
#ifdef __GPU
      CALL ppm_map_part_ghost_get(topoid,xp,ppm_dim,Np,0,ghostsize_f(1),info)
#if   __MODE == __SCA
      CALL ppm_map_part_push(up,Np,info)
#elif __MODE == __VEC
      CALL ppm_map_part_push(up,ldn,Np,info)
#endif
      CALL ppm_map_part_send(Np,Mp,info)
#if   __MODE == __SCA
      CALL ppm_map_part_pop(up,Np,Mp,info)
#elif __MODE == __VEC
      CALL ppm_map_part_pop(up,ldn,Np,Mp,info)
#endif
      CALL ppm_map_part_pop(xp,ppm_dim,Np,Mp,info)

#if   __MODE == __SCA
      CALL gpu_p2m(topoid,meshid,xp,Mp,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif __MODE == __VEC
      CALL gpu_p2m(topoid,meshid,xp,Mp,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif

#else
! ----- if GPU is NOT available ----- !
#if   __MODE == __SCA
      CALL cpu_p2m(topoid,meshid,xp,Np,up,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#elif   __MODE == __VEC
      CALL cpu_p2m(topoid,meshid,xp,Np,up,lda,kernel, &
     &           ghostsize,field_up,info,p2m_bcdef)
#endif
#endif 
      ENDIF

     !-------------------------------------------------------------------------!
     !  Return
     !-------------------------------------------------------------------------!
 9999 CONTINUE
      CALL substop('ppm_interp_p2m',t0,info)
      RETURN

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE p2m_ss_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE p2m_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE p2m_sv_2d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE p2m_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE p2m_ss_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE p2m_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
     END SUBROUTINE p2m_sv_3d
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE p2m_dv_3d
#endif
#endif
#endif

