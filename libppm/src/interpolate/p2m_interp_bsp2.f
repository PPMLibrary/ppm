      !-------------------------------------------------------------------------
      !     Subroutine   :                   p2m_interp_bsp2
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
      SUBROUTINE p2m_interp_bsp2_ss_2d(topoid,xp,up,field_up,ghostsize,dx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_ds_2d(topoid,xp,up,field_up,ghostsize,dx,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_sv_2d(topoid,xp,up,field_up,lda,ghostsize,dx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_dv_2d(topoid,xp,up,field_up,lda,ghostsize,dx,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_ss_3d(topoid,xp,up,field_up,ghostsize,dx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_ds_3d(topoid,xp,up,field_up,ghostsize,dx,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_sv_3d(topoid,xp,up,field_up,lda,ghostsize,dx,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_bsp2_dv_3d(topoid,xp,up,field_up,lda,ghostsize,dx,info)
#endif
#endif
#endif
     !!! Particle to mesh interpolation following the 2nd order B-spline scheme.
     !!!
     !!! The interpolation scheme is only implemented for 3D spaces. To
     !!! increase performance the inner loops over the number of properties to
     !!! be interpolated are unrolled for 2,3,4 and 5-vectors.
     !!! 
     !!! [NOTE]
     !!! This routine only performs the actual interpolation. It should not be
     !!! called directly from routines but instead the `ppm_interp_m2p`
     !!! routine should be used with the kernel argument set to 
     !!! `ppm_param_rmsh_kernel_bsp2`.



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


      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !------------------------------------------------------------------------!
      !  Arguments
      !------------------------------------------------------------------------!
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! topology identifier of target
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
      INTEGER                         , INTENT(IN)   :: lda
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
      REAL(MK), DIMENSION(:,:)       , INTENT(IN   ) :: xp
      !!! particle positions
      INTEGER , DIMENSION(:  )       , INTENT(IN   ) :: ghostsize
      !!! The size (width) of the ghost layer
      REAL(mk), DIMENSION(ppm_dim)   , INTENT(IN   ) :: dx
      !!! mesh spacing
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !------------------------------------------------------------------------!
      !  Local variables
      !------------------------------------------------------------------------!
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart   => NULL()
      INTEGER,  DIMENSION(:,:)     , POINTER :: ndata    => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist1   => NULL()
      INTEGER,  DIMENSION(:)       , POINTER :: ilist2   => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys => NULL()
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys => NULL()
      REAL(mk), DIMENSION(ppm_dim)           :: dxi
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
      INTEGER                                :: isub,ifrom,ito,ip,iopt,isubl
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
      CHARACTER(len=256)                     :: msg
      TYPE(ppm_t_topo)     , POINTER         :: topo => NULL()
      !------------------------------------------------------------------------!
      !  Variables for unrolled versions
      !------------------------------------------------------------------------!
      REAL(mk) :: x10,x11,x12,x13,x20,x21,x22,x23,x30,x31,x32,x33
      REAL(mk) :: a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33
      INTEGER  :: ip10,ip11,ip12,ip13,ip20,ip21,ip22,ip23,ip30,ip31,ip32,ip33
      INTEGER  :: ldn
      REAL(mk) :: a10a20a30
      REAL(mk) :: a10a20a31
      REAL(mk) :: a10a20a32
      REAL(mk) :: a10a20a33
      REAL(mk) :: a10a21a30
      REAL(mk) :: a10a21a31
      REAL(mk) :: a10a21a32
      REAL(mk) :: a10a21a33
      REAL(mk) :: a10a22a30
      REAL(mk) :: a10a22a31
      REAL(mk) :: a10a22a32
      REAL(mk) :: a10a22a33
      REAL(mk) :: a10a23a30
      REAL(mk) :: a10a23a31
      REAL(mk) :: a10a23a32
      REAL(mk) :: a10a23a33
      REAL(mk) :: a11a20a30
      REAL(mk) :: a11a20a31
      REAL(mk) :: a11a20a32
      REAL(mk) :: a11a20a33
      REAL(mk) :: a11a21a30
      REAL(mk) :: a11a21a31
      REAL(mk) :: a11a21a32
      REAL(mk) :: a11a21a33
      REAL(mk) :: a11a22a30
      REAL(mk) :: a11a22a31
      REAL(mk) :: a11a22a32
      REAL(mk) :: a11a22a33
      REAL(mk) :: a11a23a30
      REAL(mk) :: a11a23a31
      REAL(mk) :: a11a23a32
      REAL(mk) :: a11a23a33
      REAL(mk) :: a12a20a30
      REAL(mk) :: a12a20a31
      REAL(mk) :: a12a20a32
      REAL(mk) :: a12a20a33
      REAL(mk) :: a12a21a30
      REAL(mk) :: a12a21a31
      REAL(mk) :: a12a21a32
      REAL(mk) :: a12a21a33
      REAL(mk) :: a12a22a30
      REAL(mk) :: a12a22a31
      REAL(mk) :: a12a22a32
      REAL(mk) :: a12a22a33
      REAL(mk) :: a12a23a30
      REAL(mk) :: a12a23a31
      REAL(mk) :: a12a23a32
      REAL(mk) :: a12a23a33
      REAL(mk) :: a13a20a30
      REAL(mk) :: a13a20a31
      REAL(mk) :: a13a20a32
      REAL(mk) :: a13a20a33
      REAL(mk) :: a13a21a30
      REAL(mk) :: a13a21a31
      REAL(mk) :: a13a21a32
      REAL(mk) :: a13a21a33
      REAL(mk) :: a13a22a30
      REAL(mk) :: a13a22a31
      REAL(mk) :: a13a22a32
      REAL(mk) :: a13a22a33
      REAL(mk) :: a13a23a30
      REAL(mk) :: a13a23a31
      REAL(mk) :: a13a23a32
      REAL(mk) :: a13a23a33

     !-------------------------------------------------------------------------!
     !  Initialise
     !-------------------------------------------------------------------------!

      CALL substart('p2m_interp_bsp2',t0,info)

      topo => ppm_topo(topoid)%t

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

#if   __KIND == __SINGLE_PRECISION
      myeps = ppm_myepss
      min_sub => topo%min_subs
      max_sub => topo%max_subs
#elif __KIND == __DOUBLE_PRECISION
      myeps = ppm_myepsd
      min_sub => topo%min_subd
      max_sub => topo%max_subd
#endif


      dxxi = 1.0_mk/dx(1)
      dxyi = 1.0_mk/dx(2)
      IF(ppm_dim.EQ.3) dxzi = 1.0_mk/dx(3)


        !----------------------------------------------------------------------!
        !  B-spline 2 (Witch hat)
        !----------------------------------------------------------------------!
#if __DIME == __3D
        DO isub = 1,topo%nsublist

#if __MODE == __VEC

           !----------------------------------------------------------------!
           !  Unrolled versions for 3-vectors
           !----------------------------------------------------------------!
           IF (lda .EQ. 3) THEN
              isubl = topo%isublist(isub)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31


                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31

#ifdef __NOMICROINSTRUCTIONS
                field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1,iq)
                field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(2,iq)
                field_up(3,ip10,ip20,ip30,isub)=field_up(3,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(3,iq)
                field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1,iq)
                field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(2,iq)
                field_up(3,ip11,ip20,ip30,isub)=field_up(3,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(3,iq)
                field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1,iq)
                field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(2,iq)
                field_up(3,ip10,ip21,ip30,isub)=field_up(3,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(3,iq)
                field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1,iq)
                field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(2,iq)
                field_up(3,ip11,ip21,ip30,isub)=field_up(3,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(3,iq)
                field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1,iq)
                field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(2,iq)
                field_up(3,ip10,ip20,ip31,isub)=field_up(3,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(3,iq)
                field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1,iq)
                field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(2,iq)
                field_up(3,ip11,ip20,ip31,isub)=field_up(3,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(3,iq)
                field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1,iq)
                field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(2,iq)
                field_up(3,ip10,ip21,ip31,isub)=field_up(3,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(3,iq)
                field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1,iq)
                field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(2,iq)
                field_up(3,ip11,ip21,ip31,isub)=field_up(3,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(3,iq)
#else
                field_up(1:3,ip10,ip20,ip30,isub)=field_up(1:3,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1:3,iq)
                field_up(1:3,ip11,ip20,ip30,isub)=field_up(1:3,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1:3,iq)
                field_up(1:3,ip10,ip21,ip30,isub)=field_up(1:3,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1:3,iq)
                field_up(1:3,ip11,ip21,ip30,isub)=field_up(1:3,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1:3,iq)
                field_up(1:3,ip10,ip20,ip31,isub)=field_up(1:3,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1:3,iq)
                field_up(1:3,ip11,ip20,ip31,isub)=field_up(1:3,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1:3,iq)
                field_up(1:3,ip10,ip21,ip31,isub)=field_up(1:3,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1:3,iq)
                field_up(1:3,ip11,ip21,ip31,isub)=field_up(1:3,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1:3,iq)
#endif
              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 2-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 2) THEN
              isubl = topo%isublist(isub)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31

#ifdef __NOMICROINSTRUCTIONS
                field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1,iq)
                field_up(2,ip10,ip20,ip30,isub)=field_up(2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(2,iq)
                field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1,iq)
                field_up(2,ip11,ip20,ip30,isub)=field_up(2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(2,iq)
                field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1,iq)
                field_up(2,ip10,ip21,ip30,isub)=field_up(2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(2,iq)
                field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1,iq)
                field_up(2,ip11,ip21,ip30,isub)=field_up(2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(2,iq)
                field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1,iq)
                field_up(2,ip10,ip20,ip31,isub)=field_up(2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(2,iq)
                field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1,iq)
                field_up(2,ip11,ip20,ip31,isub)=field_up(2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(2,iq)
                field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1,iq)
                field_up(2,ip10,ip21,ip31,isub)=field_up(2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(2,iq)
                field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1,iq)
                field_up(2,ip11,ip21,ip31,isub)=field_up(2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(2,iq)
#else
                field_up(1:2,ip10,ip20,ip30,isub)=field_up(1:2,ip10,ip20,ip30,isub)+&
                    &a10a20a30*up(1:2,iq)
                field_up(1:2,ip11,ip20,ip30,isub)=field_up(1:2,ip11,ip20,ip30,isub)+&
                    &a11a20a30*up(1:2,iq)
                field_up(1:2,ip10,ip21,ip30,isub)=field_up(1:2,ip10,ip21,ip30,isub)+&
                    &a10a21a30*up(1:2,iq)
                field_up(1:2,ip11,ip21,ip30,isub)=field_up(1:2,ip11,ip21,ip30,isub)+&
                    &a11a21a30*up(1:2,iq)
                field_up(1:2,ip10,ip20,ip31,isub)=field_up(1:2,ip10,ip20,ip31,isub)+&
                    &a10a20a31*up(1:2,iq)
                field_up(1:2,ip11,ip20,ip31,isub)=field_up(1:2,ip11,ip20,ip31,isub)+&
                    &a11a20a31*up(1:2,iq)
                field_up(1:2,ip10,ip21,ip31,isub)=field_up(1:2,ip10,ip21,ip31,isub)+&
                    &a10a21a31*up(1:2,iq)
                field_up(1:2,ip11,ip21,ip31,isub)=field_up(1:2,ip11,ip21,ip31,isub)+&
                    &a11a21a31*up(1:2,iq)
#endif
              ENDDO
              !----------------------------------------------------------------!
              !  Unrolled versions for 1-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 1) THEN
              isubl = topo%isublist(isub)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31

             field_up(1,ip10,ip20,ip30,isub)=field_up(1,ip10,ip20,ip30,isub)+&
     &                  a10a20a30*up(1,iq)
             field_up(1,ip11,ip20,ip30,isub)=field_up(1,ip11,ip20,ip30,isub)+&
     &                  a11a20a30*up(1,iq)
             field_up(1,ip10,ip21,ip30,isub)=field_up(1,ip10,ip21,ip30,isub)+&
     &                  a10a21a30*up(1,iq)
             field_up(1,ip11,ip21,ip30,isub)=field_up(1,ip11,ip21,ip30,isub)+&
     &                  a11a21a30*up(1,iq)
             field_up(1,ip10,ip20,ip31,isub)=field_up(1,ip10,ip20,ip31,isub)+&
     &                  a10a20a31*up(1,iq)
             field_up(1,ip11,ip20,ip31,isub)=field_up(1,ip11,ip20,ip31,isub)+&
     &                  a11a20a31*up(1,iq)
             field_up(1,ip10,ip21,ip31,isub)=field_up(1,ip10,ip21,ip31,isub)+&
     &                  a10a21a31*up(1,iq)
             field_up(1,ip11,ip21,ip31,isub)=field_up(1,ip11,ip21,ip31,isub)+&
     &                  a11a21a31*up(1,iq)

              ENDDO
              !----------------------------------------------------------------!
              !  All other lda are NOT UNROLLED. Vectorization will be over lda!
              !----------------------------------------------------------------!
           ELSE
              isubl = topo%isublist(isub)
              DO ip = 1,store_info(isub)

                 iq    = list_sub(isub,ip)

                 x01 = (xp(1,iq)-min_sub(1,isubl))*dxxi
                 x02 = (xp(2,iq)-min_sub(2,isubl))*dxyi
                 x03 = (xp(3,iq)-min_sub(3,isubl))*dxzi

                 ip10 = INT(x01) + 1
                 ip20 = INT(x02) + 1
                 ip30 = INT(x03) + 1

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 xp1 = x01-REAL(ip10-1,mk)
                 xp2 = x02-REAL(ip20-1,mk)
                 xp3 = x03-REAL(ip30-1,mk)

                 x10 = xp1
                 x11 = x10 - 1.0_mk

                 x20 = xp2
                 x21 = x20 - 1.0_mk

                 x30 = xp3
                 x31 = x30 - 1.0_mk

                 a10 = 1.0_mk - x10
                 a20 = 1.0_mk - x20
                 a30 = 1.0_mk - x30

                 a11 = 1.0_mk + x11
                 a21 = 1.0_mk + x21
                 a31 = 1.0_mk + x31

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 DO ldn=1,lda

          field_up(ldn,ip10,ip20,ip30,isub)=field_up(ldn,ip10,ip20,ip30,isub)+&
     &               a10a20a30*up(ldn,iq)
          field_up(ldn,ip11,ip20,ip30,isub)=field_up(ldn,ip11,ip20,ip30,isub)+&
     &               a11a20a30*up(ldn,iq)
          field_up(ldn,ip10,ip21,ip30,isub)=field_up(ldn,ip10,ip21,ip30,isub)+&
     &               a10a21a30*up(ldn,iq)
          field_up(ldn,ip11,ip21,ip30,isub)=field_up(ldn,ip11,ip21,ip30,isub)+&
     &               a11a21a30*up(ldn,iq)
          field_up(ldn,ip10,ip20,ip31,isub)=field_up(ldn,ip10,ip20,ip31,isub)+&
     &               a10a20a31*up(ldn,iq)
          field_up(ldn,ip11,ip20,ip31,isub)=field_up(ldn,ip11,ip20,ip31,isub)+&
     &               a11a20a31*up(ldn,iq)
          field_up(ldn,ip10,ip21,ip31,isub)=field_up(ldn,ip10,ip21,ip31,isub)+&
     &               a10a21a31*up(ldn,iq)
          field_up(ldn,ip11,ip21,ip31,isub)=field_up(ldn,ip11,ip21,ip31,isub)+&
     &               a11a21a31*up(ldn,iq)
                 ENDDO    ! lda
              ENDDO        ! iq
           END IF          ! unrolled lda cases
#elif __MODE == __SCA
           isubl = topo%isublist(isub)
           DO ip = 1,store_info(isub)

              iq    = list_sub(isub,ip)

              x01 = (xp(1,iq)-min_sub(1,isubl))*dxxi
              x02 = (xp(2,iq)-min_sub(2,isubl))*dxyi
              x03 = (xp(3,iq)-min_sub(3,isubl))*dxzi

              ip10 = INT(x01) + 1
              ip20 = INT(x02) + 1
              ip30 = INT(x03) + 1

              ip11 = ip10 + 1
              ip21 = ip20 + 1
              ip31 = ip30 + 1

              xp1 = x01-REAL(ip10-1,mk)
              xp2 = x02-REAL(ip20-1,mk)
              xp3 = x03-REAL(ip30-1,mk)

              x10 = xp1
              x11 = x10 - 1.0_mk

              x20 = xp2
              x21 = x20 - 1.0_mk

              x30 = xp3
              x31 = x30 - 1.0_mk

              a10 = 1.0_mk - x10
              a20 = 1.0_mk - x20
              a30 = 1.0_mk - x30

              a11 = 1.0_mk + x11
              a21 = 1.0_mk + x21
              a31 = 1.0_mk + x31

              a10a20a30 = a10*a20*a30
              a10a20a31 = a10*a20*a31
              a10a21a30 = a10*a21*a30
              a10a21a31 = a10*a21*a31
              a11a20a30 = a11*a20*a30
              a11a20a31 = a11*a20*a31
              a11a21a30 = a11*a21*a30
              a11a21a31 = a11*a21*a31


              field_up(ip10,ip20,ip30,isub)=field_up(ip10,ip20,ip30,isub)+&
     &                   a10a20a30*up(iq)
              field_up(ip11,ip20,ip30,isub)=field_up(ip11,ip20,ip30,isub)+&
     &                   a11a20a30*up(iq)
              field_up(ip10,ip21,ip30,isub)=field_up(ip10,ip21,ip30,isub)+&
     &                   a10a21a30*up(iq)
              field_up(ip11,ip21,ip30,isub)=field_up(ip11,ip21,ip30,isub)+&
     &                   a11a21a30*up(iq)
              field_up(ip10,ip20,ip31,isub)=field_up(ip10,ip20,ip31,isub)+&
     &                   a10a20a31*up(iq)
              field_up(ip11,ip20,ip31,isub)=field_up(ip11,ip20,ip31,isub)+&
     &                   a11a20a31*up(iq)
              field_up(ip10,ip21,ip31,isub)=field_up(ip10,ip21,ip31,isub)+&
     &                   a10a21a31*up(iq)
              field_up(ip11,ip21,ip31,isub)=field_up(ip11,ip21,ip31,isub)+&
     &                   a11a21a31*up(iq)
           ENDDO        ! iq
#endif
        END DO              ! loop over subs
#elif __DIME == __2D
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,'p2m_interp_bsp2',    &
     &        'Currently bsp2 is not available. Use ppm_rmsh_remesh for other kernels.', &
     &           __LINE__,info)
#endif
      CALL substop('p2m_interp_bsp2',t0,info)
      RETURN

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_bsp2_dv_3d
#endif
#endif
#endif
