      !-------------------------------------------------------------------------
      !     Subroutine   :                   p2m_interp_mp4
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

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_mp4_ss_2d(Mesh,Field,field_up,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_mp4_ds_2d(Mesh,Field,field_up,xp,up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_mp4_sv_2d(Mesh,Field,field_up,lda,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_mp4_dv_2d(Mesh,Field,field_up,lda,xp,up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_mp4_ss_3d(Mesh,Field,field_up,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_mp4_ds_3d(Mesh,Field,field_up,xp,up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE p2m_interp_mp4_sv_3d(Mesh,Field,field_up,lda,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE p2m_interp_mp4_dv_3d(Mesh,Field,field_up,lda,xp,up,info)
#endif
#endif
#endif
     !!! Particle to mesh interpolation following the M-prime-4 scheme.
     !!!
     !!! The interpolation scheme is only implemented for 2D and 3D spaces. To
     !!! increase performance the inner loops over the number of properties to
     !!! be interpolated are unrolled for 2,3,4 and 5-vectors.
     !!! 
     !!! [NOTE]
     !!! This routine only performs the actual interpolation. It should not be
     !!! called directly by the user but instead the `ppm_interp_m2p`
     !!! routine should be used with the kernel argument set to 
     !!! `ppm_param_rmsh_kernel_mp4`.



      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_substart
      USE ppm_module_substop
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
      CLASS(ppm_t_equi_mesh_)                      :: Mesh
      !!! Mesh
      CLASS(ppm_t_field_)                          :: Field
      !!! Field
#if   __MODE == __SCA
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:    ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:  ) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
#elif __MODE == __VEC
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:  ) , POINTER      :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:) , POINTER      :: field_up
#endif
      !!! field onto which to interpolate
      INTEGER                       , INTENT(IN)   :: lda
      !!! leading dimension of up
#endif
      REAL(MK), DIMENSION(:,:), POINTER,INTENT(IN)   :: xp 
      !!! particle positions
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:),  POINTER,INTENT(IN)   :: up
#elif __MODE == __VEC
      REAL(MK) , DIMENSION(:,:),POINTER,INTENT(IN)   :: up 
#endif
      INTEGER                       , INTENT( OUT) :: info
      !!! Returns status, 0 upon success
      !------------------------------------------------------------------------!
      !  Local variables
      !------------------------------------------------------------------------!
      REAL(mk)                               :: dxxi,dxyi,dxzi
      REAL(mk)                               :: x1,x2,x3
      INTEGER                                :: i,j,k,l,ii,jj,kk
      INTEGER                                :: ip1,ip2,ip3
      INTEGER                                :: ip
      INTEGER, DIMENSION(6)                  :: bcdef
      INTEGER                                :: iq,nsubpatch,ipatch
      ! aliases
      REAL(mk)                               :: tim1s, tim1e
      REAL(mk)                               :: xp1,xp2,xp3
      REAL(mk)                               :: wx1,wx2,wx3
      REAL(mk), DIMENSION(ppm_dim)           :: x0
      REAL(mk)                               :: x01,x02,x03
      CHARACTER(len=256)                     :: msg
      CLASS(ppm_t_subpatch_),POINTER         :: p => NULL()
      REAL(mk)                               :: x01_l,x02_l,x03_l
      REAL(mk)                               :: x01_h,x02_h,x03_h
      REAL(mk)                               :: l10,l11,l12,l13
      REAL(mk)                               :: h10,h11,h12,h13
      REAL(mk)                               :: l20,l21,l22,l23
      REAL(mk)                               :: h20,h21,h22,h23
      REAL(mk)                               :: l30,l31,l32,l33
      REAL(mk)                               :: h30,h31,h32,h33
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


      start_subroutine("p2m_interp_mp4")

      dxxi = 1.0_mk/Mesh%h(1)
      dxyi = 1.0_mk/Mesh%h(2)
      IF(ppm_dim.EQ.3) dxzi = 1.0_mk/Mesh%h(3)


      !  loop over subpatches
      p => Mesh%subpatch%begin()
      ipatch = 1
      subpatch: DO WHILE (ASSOCIATED(p))
          CALL p%get_field(field_up,Field,info)
              or_fail("get_field failed for this subpatch")

#if __DIME == __3D

#if __MODE == __VEC

              !----------------------------------------------------------------!
              !  Unrolled versions for 3-vectors
              !----------------------------------------------------------------!
           IF (lda .EQ. 3) THEN
              DO ip = 1,store_info(ipatch)

                 iq    = list_sub(ipatch,ip)

                 x01 = xp(1,iq)*dxxi-p%istart(1) + 1
                 x02 = xp(2,iq)*dxyi-p%istart(2) + 1
                 x03 = xp(3,iq)*dxzi-p%istart(3) + 1

                 ip10 = FLOOR(x01)
                 ip20 = FLOOR(x02)
                 ip30 = FLOOR(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

#ifdef __NOMICROINSTRUCTIONS
      field_up(1,ip10,ip20,ip30)=field_up(1,ip10,ip20,ip30)+&
      &a10a20a30*up(1,iq)
      field_up(2,ip10,ip20,ip30)=field_up(2,ip10,ip20,ip30)+&
      &a10a20a30*up(2,iq)
      field_up(3,ip10,ip20,ip30)=field_up(3,ip10,ip20,ip30)+&
      &a10a20a30*up(3,iq)
      field_up(1,ip11,ip20,ip30)=field_up(1,ip11,ip20,ip30)+&
      &a11a20a30*up(1,iq)
      field_up(2,ip11,ip20,ip30)=field_up(2,ip11,ip20,ip30)+&
      &a11a20a30*up(2,iq)
      field_up(3,ip11,ip20,ip30)=field_up(3,ip11,ip20,ip30)+&
      &a11a20a30*up(3,iq)
      field_up(1,ip12,ip20,ip30)=field_up(1,ip12,ip20,ip30)+&
      &a12a20a30*up(1,iq)
      field_up(2,ip12,ip20,ip30)=field_up(2,ip12,ip20,ip30)+&
      &a12a20a30*up(2,iq)
      field_up(3,ip12,ip20,ip30)=field_up(3,ip12,ip20,ip30)+&
      &a12a20a30*up(3,iq)
      field_up(1,ip13,ip20,ip30)=field_up(1,ip13,ip20,ip30)+&
      &a13a20a30*up(1,iq)
      field_up(2,ip13,ip20,ip30)=field_up(2,ip13,ip20,ip30)+&
      &a13a20a30*up(2,iq)
      field_up(3,ip13,ip20,ip30)=field_up(3,ip13,ip20,ip30)+&
      &a13a20a30*up(3,iq)
      field_up(1,ip10,ip21,ip30)=field_up(1,ip10,ip21,ip30)+&
      &a10a21a30*up(1,iq)
      field_up(2,ip10,ip21,ip30)=field_up(2,ip10,ip21,ip30)+&
      &a10a21a30*up(2,iq)
      field_up(3,ip10,ip21,ip30)=field_up(3,ip10,ip21,ip30)+&
      &a10a21a30*up(3,iq)
      field_up(1,ip11,ip21,ip30)=field_up(1,ip11,ip21,ip30)+&
      &a11a21a30*up(1,iq)
      field_up(2,ip11,ip21,ip30)=field_up(2,ip11,ip21,ip30)+&
      &a11a21a30*up(2,iq)
      field_up(3,ip11,ip21,ip30)=field_up(3,ip11,ip21,ip30)+&
      &a11a21a30*up(3,iq)
      field_up(1,ip12,ip21,ip30)=field_up(1,ip12,ip21,ip30)+&
      &a12a21a30*up(1,iq)
      field_up(2,ip12,ip21,ip30)=field_up(2,ip12,ip21,ip30)+&
      &a12a21a30*up(2,iq)
      field_up(3,ip12,ip21,ip30)=field_up(3,ip12,ip21,ip30)+&
      &a12a21a30*up(3,iq)
      field_up(1,ip13,ip21,ip30)=field_up(1,ip13,ip21,ip30)+&
      &a13a21a30*up(1,iq)
      field_up(2,ip13,ip21,ip30)=field_up(2,ip13,ip21,ip30)+&
      &a13a21a30*up(2,iq)
      field_up(3,ip13,ip21,ip30)=field_up(3,ip13,ip21,ip30)+&
      &a13a21a30*up(3,iq)
      field_up(1,ip10,ip22,ip30)=field_up(1,ip10,ip22,ip30)+&
      &a10a22a30*up(1,iq)
      field_up(2,ip10,ip22,ip30)=field_up(2,ip10,ip22,ip30)+&
      &a10a22a30*up(2,iq)
      field_up(3,ip10,ip22,ip30)=field_up(3,ip10,ip22,ip30)+&
      &a10a22a30*up(3,iq)
      field_up(1,ip11,ip22,ip30)=field_up(1,ip11,ip22,ip30)+&
      &a11a22a30*up(1,iq)
      field_up(2,ip11,ip22,ip30)=field_up(2,ip11,ip22,ip30)+&
      &a11a22a30*up(2,iq)
      field_up(3,ip11,ip22,ip30)=field_up(3,ip11,ip22,ip30)+&
      &a11a22a30*up(3,iq)
      field_up(1,ip12,ip22,ip30)=field_up(1,ip12,ip22,ip30)+&
      &a12a22a30*up(1,iq)
      field_up(2,ip12,ip22,ip30)=field_up(2,ip12,ip22,ip30)+&
      &a12a22a30*up(2,iq)
      field_up(3,ip12,ip22,ip30)=field_up(3,ip12,ip22,ip30)+&
      &a12a22a30*up(3,iq)
      field_up(1,ip13,ip22,ip30)=field_up(1,ip13,ip22,ip30)+&
      &a13a22a30*up(1,iq)
      field_up(2,ip13,ip22,ip30)=field_up(2,ip13,ip22,ip30)+&
      &a13a22a30*up(2,iq)
      field_up(3,ip13,ip22,ip30)=field_up(3,ip13,ip22,ip30)+&
      &a13a22a30*up(3,iq)
      field_up(1,ip10,ip23,ip30)=field_up(1,ip10,ip23,ip30)+&
      &a10a23a30*up(1,iq)
      field_up(2,ip10,ip23,ip30)=field_up(2,ip10,ip23,ip30)+&
      &a10a23a30*up(2,iq)
      field_up(3,ip10,ip23,ip30)=field_up(3,ip10,ip23,ip30)+&
      &a10a23a30*up(3,iq)
      field_up(1,ip11,ip23,ip30)=field_up(1,ip11,ip23,ip30)+&
      &a11a23a30*up(1,iq)
      field_up(2,ip11,ip23,ip30)=field_up(2,ip11,ip23,ip30)+&
      &a11a23a30*up(2,iq)
      field_up(3,ip11,ip23,ip30)=field_up(3,ip11,ip23,ip30)+&
      &a11a23a30*up(3,iq)
      field_up(1,ip12,ip23,ip30)=field_up(1,ip12,ip23,ip30)+&
      &a12a23a30*up(1,iq)
      field_up(2,ip12,ip23,ip30)=field_up(2,ip12,ip23,ip30)+&
      &a12a23a30*up(2,iq)
      field_up(3,ip12,ip23,ip30)=field_up(3,ip12,ip23,ip30)+&
      &a12a23a30*up(3,iq)
      field_up(1,ip13,ip23,ip30)=field_up(1,ip13,ip23,ip30)+&
      &a13a23a30*up(1,iq)
      field_up(2,ip13,ip23,ip30)=field_up(2,ip13,ip23,ip30)+&
      &a13a23a30*up(2,iq)
      field_up(3,ip13,ip23,ip30)=field_up(3,ip13,ip23,ip30)+&
      &a13a23a30*up(3,iq)
      field_up(1,ip10,ip20,ip31)=field_up(1,ip10,ip20,ip31)+&
      &a10a20a31*up(1,iq)
      field_up(2,ip10,ip20,ip31)=field_up(2,ip10,ip20,ip31)+&
      &a10a20a31*up(2,iq)
      field_up(3,ip10,ip20,ip31)=field_up(3,ip10,ip20,ip31)+&
      &a10a20a31*up(3,iq)
      field_up(1,ip11,ip20,ip31)=field_up(1,ip11,ip20,ip31)+&
      &a11a20a31*up(1,iq)
      field_up(2,ip11,ip20,ip31)=field_up(2,ip11,ip20,ip31)+&
      &a11a20a31*up(2,iq)
      field_up(3,ip11,ip20,ip31)=field_up(3,ip11,ip20,ip31)+&
      &a11a20a31*up(3,iq)
      field_up(1,ip12,ip20,ip31)=field_up(1,ip12,ip20,ip31)+&
      &a12a20a31*up(1,iq)
      field_up(2,ip12,ip20,ip31)=field_up(2,ip12,ip20,ip31)+&
      &a12a20a31*up(2,iq)
      field_up(3,ip12,ip20,ip31)=field_up(3,ip12,ip20,ip31)+&
      &a12a20a31*up(3,iq)
      field_up(1,ip13,ip20,ip31)=field_up(1,ip13,ip20,ip31)+&
      &a13a20a31*up(1,iq)
      field_up(2,ip13,ip20,ip31)=field_up(2,ip13,ip20,ip31)+&
      &a13a20a31*up(2,iq)
      field_up(3,ip13,ip20,ip31)=field_up(3,ip13,ip20,ip31)+&
      &a13a20a31*up(3,iq)
      field_up(1,ip10,ip21,ip31)=field_up(1,ip10,ip21,ip31)+&
      &a10a21a31*up(1,iq)
      field_up(2,ip10,ip21,ip31)=field_up(2,ip10,ip21,ip31)+&
      &a10a21a31*up(2,iq)
      field_up(3,ip10,ip21,ip31)=field_up(3,ip10,ip21,ip31)+&
      &a10a21a31*up(3,iq)
      field_up(1,ip11,ip21,ip31)=field_up(1,ip11,ip21,ip31)+&
      &a11a21a31*up(1,iq)
      field_up(2,ip11,ip21,ip31)=field_up(2,ip11,ip21,ip31)+&
      &a11a21a31*up(2,iq)
      field_up(3,ip11,ip21,ip31)=field_up(3,ip11,ip21,ip31)+&
      &a11a21a31*up(3,iq)
      field_up(1,ip12,ip21,ip31)=field_up(1,ip12,ip21,ip31)+&
      &a12a21a31*up(1,iq)
      field_up(2,ip12,ip21,ip31)=field_up(2,ip12,ip21,ip31)+&
      &a12a21a31*up(2,iq)
      field_up(3,ip12,ip21,ip31)=field_up(3,ip12,ip21,ip31)+&
      &a12a21a31*up(3,iq)
      field_up(1,ip13,ip21,ip31)=field_up(1,ip13,ip21,ip31)+&
      &a13a21a31*up(1,iq)
      field_up(2,ip13,ip21,ip31)=field_up(2,ip13,ip21,ip31)+&
      &a13a21a31*up(2,iq)
      field_up(3,ip13,ip21,ip31)=field_up(3,ip13,ip21,ip31)+&
      &a13a21a31*up(3,iq)
      field_up(1,ip10,ip22,ip31)=field_up(1,ip10,ip22,ip31)+&
      &a10a22a31*up(1,iq)
      field_up(2,ip10,ip22,ip31)=field_up(2,ip10,ip22,ip31)+&
      &a10a22a31*up(2,iq)
      field_up(3,ip10,ip22,ip31)=field_up(3,ip10,ip22,ip31)+&
      &a10a22a31*up(3,iq)
      field_up(1,ip11,ip22,ip31)=field_up(1,ip11,ip22,ip31)+&
      &a11a22a31*up(1,iq)
      field_up(2,ip11,ip22,ip31)=field_up(2,ip11,ip22,ip31)+&
      &a11a22a31*up(2,iq)
      field_up(3,ip11,ip22,ip31)=field_up(3,ip11,ip22,ip31)+&
      &a11a22a31*up(3,iq)
      field_up(1,ip12,ip22,ip31)=field_up(1,ip12,ip22,ip31)+&
      &a12a22a31*up(1,iq)
      field_up(2,ip12,ip22,ip31)=field_up(2,ip12,ip22,ip31)+&
      &a12a22a31*up(2,iq)
      field_up(3,ip12,ip22,ip31)=field_up(3,ip12,ip22,ip31)+&
      &a12a22a31*up(3,iq)
      field_up(1,ip13,ip22,ip31)=field_up(1,ip13,ip22,ip31)+&
      &a13a22a31*up(1,iq)
      field_up(2,ip13,ip22,ip31)=field_up(2,ip13,ip22,ip31)+&
      &a13a22a31*up(2,iq)
      field_up(3,ip13,ip22,ip31)=field_up(3,ip13,ip22,ip31)+&
      &a13a22a31*up(3,iq)
      field_up(1,ip10,ip23,ip31)=field_up(1,ip10,ip23,ip31)+&
      &a10a23a31*up(1,iq)
      field_up(2,ip10,ip23,ip31)=field_up(2,ip10,ip23,ip31)+&
      &a10a23a31*up(2,iq)
      field_up(3,ip10,ip23,ip31)=field_up(3,ip10,ip23,ip31)+&
      &a10a23a31*up(3,iq)
      field_up(1,ip11,ip23,ip31)=field_up(1,ip11,ip23,ip31)+&
      &a11a23a31*up(1,iq)
      field_up(2,ip11,ip23,ip31)=field_up(2,ip11,ip23,ip31)+&
      &a11a23a31*up(2,iq)
      field_up(3,ip11,ip23,ip31)=field_up(3,ip11,ip23,ip31)+&
      &a11a23a31*up(3,iq)
      field_up(1,ip12,ip23,ip31)=field_up(1,ip12,ip23,ip31)+&
      &a12a23a31*up(1,iq)
      field_up(2,ip12,ip23,ip31)=field_up(2,ip12,ip23,ip31)+&
      &a12a23a31*up(2,iq)
      field_up(3,ip12,ip23,ip31)=field_up(3,ip12,ip23,ip31)+&
      &a12a23a31*up(3,iq)
      field_up(1,ip13,ip23,ip31)=field_up(1,ip13,ip23,ip31)+&
      &a13a23a31*up(1,iq)
      field_up(2,ip13,ip23,ip31)=field_up(2,ip13,ip23,ip31)+&
      &a13a23a31*up(2,iq)
      field_up(3,ip13,ip23,ip31)=field_up(3,ip13,ip23,ip31)+&
      &a13a23a31*up(3,iq)
      field_up(1,ip10,ip20,ip32)=field_up(1,ip10,ip20,ip32)+&
      &a10a20a32*up(1,iq)
      field_up(2,ip10,ip20,ip32)=field_up(2,ip10,ip20,ip32)+&
      &a10a20a32*up(2,iq)
      field_up(3,ip10,ip20,ip32)=field_up(3,ip10,ip20,ip32)+&
      &a10a20a32*up(3,iq)
      field_up(1,ip11,ip20,ip32)=field_up(1,ip11,ip20,ip32)+&
      &a11a20a32*up(1,iq)
      field_up(2,ip11,ip20,ip32)=field_up(2,ip11,ip20,ip32)+&
      &a11a20a32*up(2,iq)
      field_up(3,ip11,ip20,ip32)=field_up(3,ip11,ip20,ip32)+&
      &a11a20a32*up(3,iq)
      field_up(1,ip12,ip20,ip32)=field_up(1,ip12,ip20,ip32)+&
      &a12a20a32*up(1,iq)
      field_up(2,ip12,ip20,ip32)=field_up(2,ip12,ip20,ip32)+&
      &a12a20a32*up(2,iq)
      field_up(3,ip12,ip20,ip32)=field_up(3,ip12,ip20,ip32)+&
      &a12a20a32*up(3,iq)
      field_up(1,ip13,ip20,ip32)=field_up(1,ip13,ip20,ip32)+&
      &a13a20a32*up(1,iq)
      field_up(2,ip13,ip20,ip32)=field_up(2,ip13,ip20,ip32)+&
      &a13a20a32*up(2,iq)
      field_up(3,ip13,ip20,ip32)=field_up(3,ip13,ip20,ip32)+&
      &a13a20a32*up(3,iq)
      field_up(1,ip10,ip21,ip32)=field_up(1,ip10,ip21,ip32)+&
      &a10a21a32*up(1,iq)
      field_up(2,ip10,ip21,ip32)=field_up(2,ip10,ip21,ip32)+&
      &a10a21a32*up(2,iq)
      field_up(3,ip10,ip21,ip32)=field_up(3,ip10,ip21,ip32)+&
      &a10a21a32*up(3,iq)
      field_up(1,ip11,ip21,ip32)=field_up(1,ip11,ip21,ip32)+&
      &a11a21a32*up(1,iq)
      field_up(2,ip11,ip21,ip32)=field_up(2,ip11,ip21,ip32)+&
      &a11a21a32*up(2,iq)
      field_up(3,ip11,ip21,ip32)=field_up(3,ip11,ip21,ip32)+&
      &a11a21a32*up(3,iq)
      field_up(1,ip12,ip21,ip32)=field_up(1,ip12,ip21,ip32)+&
      &a12a21a32*up(1,iq)
      field_up(2,ip12,ip21,ip32)=field_up(2,ip12,ip21,ip32)+&
      &a12a21a32*up(2,iq)
      field_up(3,ip12,ip21,ip32)=field_up(3,ip12,ip21,ip32)+&
      &a12a21a32*up(3,iq)
      field_up(1,ip13,ip21,ip32)=field_up(1,ip13,ip21,ip32)+&
      &a13a21a32*up(1,iq)
      field_up(2,ip13,ip21,ip32)=field_up(2,ip13,ip21,ip32)+&
      &a13a21a32*up(2,iq)
      field_up(3,ip13,ip21,ip32)=field_up(3,ip13,ip21,ip32)+&
      &a13a21a32*up(3,iq)
      field_up(1,ip10,ip22,ip32)=field_up(1,ip10,ip22,ip32)+&
      &a10a22a32*up(1,iq)
      field_up(2,ip10,ip22,ip32)=field_up(2,ip10,ip22,ip32)+&
      &a10a22a32*up(2,iq)
      field_up(3,ip10,ip22,ip32)=field_up(3,ip10,ip22,ip32)+&
      &a10a22a32*up(3,iq)
      field_up(1,ip11,ip22,ip32)=field_up(1,ip11,ip22,ip32)+&
      &a11a22a32*up(1,iq)
      field_up(2,ip11,ip22,ip32)=field_up(2,ip11,ip22,ip32)+&
      &a11a22a32*up(2,iq)
      field_up(3,ip11,ip22,ip32)=field_up(3,ip11,ip22,ip32)+&
      &a11a22a32*up(3,iq)
      field_up(1,ip12,ip22,ip32)=field_up(1,ip12,ip22,ip32)+&
      &a12a22a32*up(1,iq)
      field_up(2,ip12,ip22,ip32)=field_up(2,ip12,ip22,ip32)+&
      &a12a22a32*up(2,iq)
      field_up(3,ip12,ip22,ip32)=field_up(3,ip12,ip22,ip32)+&
      &a12a22a32*up(3,iq)
      field_up(1,ip13,ip22,ip32)=field_up(1,ip13,ip22,ip32)+&
      &a13a22a32*up(1,iq)
      field_up(2,ip13,ip22,ip32)=field_up(2,ip13,ip22,ip32)+&
      &a13a22a32*up(2,iq)
      field_up(3,ip13,ip22,ip32)=field_up(3,ip13,ip22,ip32)+&
      &a13a22a32*up(3,iq)
      field_up(1,ip10,ip23,ip32)=field_up(1,ip10,ip23,ip32)+&
      &a10a23a32*up(1,iq)
      field_up(2,ip10,ip23,ip32)=field_up(2,ip10,ip23,ip32)+&
      &a10a23a32*up(2,iq)
      field_up(3,ip10,ip23,ip32)=field_up(3,ip10,ip23,ip32)+&
      &a10a23a32*up(3,iq)
      field_up(1,ip11,ip23,ip32)=field_up(1,ip11,ip23,ip32)+&
      &a11a23a32*up(1,iq)
      field_up(2,ip11,ip23,ip32)=field_up(2,ip11,ip23,ip32)+&
      &a11a23a32*up(2,iq)
      field_up(3,ip11,ip23,ip32)=field_up(3,ip11,ip23,ip32)+&
      &a11a23a32*up(3,iq)
      field_up(1,ip12,ip23,ip32)=field_up(1,ip12,ip23,ip32)+&
      &a12a23a32*up(1,iq)
      field_up(2,ip12,ip23,ip32)=field_up(2,ip12,ip23,ip32)+&
      &a12a23a32*up(2,iq)
      field_up(3,ip12,ip23,ip32)=field_up(3,ip12,ip23,ip32)+&
      &a12a23a32*up(3,iq)
      field_up(1,ip13,ip23,ip32)=field_up(1,ip13,ip23,ip32)+&
      &a13a23a32*up(1,iq)
      field_up(2,ip13,ip23,ip32)=field_up(2,ip13,ip23,ip32)+&
      &a13a23a32*up(2,iq)
      field_up(3,ip13,ip23,ip32)=field_up(3,ip13,ip23,ip32)+&
      &a13a23a32*up(3,iq)
      field_up(1,ip10,ip20,ip33)=field_up(1,ip10,ip20,ip33)+&
      &a10a20a33*up(1,iq)
      field_up(2,ip10,ip20,ip33)=field_up(2,ip10,ip20,ip33)+&
      &a10a20a33*up(2,iq)
      field_up(3,ip10,ip20,ip33)=field_up(3,ip10,ip20,ip33)+&
      &a10a20a33*up(3,iq)
      field_up(1,ip11,ip20,ip33)=field_up(1,ip11,ip20,ip33)+&
      &a11a20a33*up(1,iq)
      field_up(2,ip11,ip20,ip33)=field_up(2,ip11,ip20,ip33)+&
      &a11a20a33*up(2,iq)
      field_up(3,ip11,ip20,ip33)=field_up(3,ip11,ip20,ip33)+&
      &a11a20a33*up(3,iq)
      field_up(1,ip12,ip20,ip33)=field_up(1,ip12,ip20,ip33)+&
      &a12a20a33*up(1,iq)
      field_up(2,ip12,ip20,ip33)=field_up(2,ip12,ip20,ip33)+&
      &a12a20a33*up(2,iq)
      field_up(3,ip12,ip20,ip33)=field_up(3,ip12,ip20,ip33)+&
      &a12a20a33*up(3,iq)
      field_up(1,ip13,ip20,ip33)=field_up(1,ip13,ip20,ip33)+&
      &a13a20a33*up(1,iq)
      field_up(2,ip13,ip20,ip33)=field_up(2,ip13,ip20,ip33)+&
      &a13a20a33*up(2,iq)
      field_up(3,ip13,ip20,ip33)=field_up(3,ip13,ip20,ip33)+&
      &a13a20a33*up(3,iq)
      field_up(1,ip10,ip21,ip33)=field_up(1,ip10,ip21,ip33)+&
      &a10a21a33*up(1,iq)
      field_up(2,ip10,ip21,ip33)=field_up(2,ip10,ip21,ip33)+&
      &a10a21a33*up(2,iq)
      field_up(3,ip10,ip21,ip33)=field_up(3,ip10,ip21,ip33)+&
      &a10a21a33*up(3,iq)
      field_up(1,ip11,ip21,ip33)=field_up(1,ip11,ip21,ip33)+&
      &a11a21a33*up(1,iq)
      field_up(2,ip11,ip21,ip33)=field_up(2,ip11,ip21,ip33)+&
      &a11a21a33*up(2,iq)
      field_up(3,ip11,ip21,ip33)=field_up(3,ip11,ip21,ip33)+&
      &a11a21a33*up(3,iq)
      field_up(1,ip12,ip21,ip33)=field_up(1,ip12,ip21,ip33)+&
      &a12a21a33*up(1,iq)
      field_up(2,ip12,ip21,ip33)=field_up(2,ip12,ip21,ip33)+&
      &a12a21a33*up(2,iq)
      field_up(3,ip12,ip21,ip33)=field_up(3,ip12,ip21,ip33)+&
      &a12a21a33*up(3,iq)
      field_up(1,ip13,ip21,ip33)=field_up(1,ip13,ip21,ip33)+&
      &a13a21a33*up(1,iq)
      field_up(2,ip13,ip21,ip33)=field_up(2,ip13,ip21,ip33)+&
      &a13a21a33*up(2,iq)
      field_up(3,ip13,ip21,ip33)=field_up(3,ip13,ip21,ip33)+&
      &a13a21a33*up(3,iq)
      field_up(1,ip10,ip22,ip33)=field_up(1,ip10,ip22,ip33)+&
      &a10a22a33*up(1,iq)
      field_up(2,ip10,ip22,ip33)=field_up(2,ip10,ip22,ip33)+&
      &a10a22a33*up(2,iq)
      field_up(3,ip10,ip22,ip33)=field_up(3,ip10,ip22,ip33)+&
      &a10a22a33*up(3,iq)
      field_up(1,ip11,ip22,ip33)=field_up(1,ip11,ip22,ip33)+&
      &a11a22a33*up(1,iq)
      field_up(2,ip11,ip22,ip33)=field_up(2,ip11,ip22,ip33)+&
      &a11a22a33*up(2,iq)
      field_up(3,ip11,ip22,ip33)=field_up(3,ip11,ip22,ip33)+&
      &a11a22a33*up(3,iq)
      field_up(1,ip12,ip22,ip33)=field_up(1,ip12,ip22,ip33)+&
      &a12a22a33*up(1,iq)
      field_up(2,ip12,ip22,ip33)=field_up(2,ip12,ip22,ip33)+&
      &a12a22a33*up(2,iq)
      field_up(3,ip12,ip22,ip33)=field_up(3,ip12,ip22,ip33)+&
      &a12a22a33*up(3,iq)
      field_up(1,ip13,ip22,ip33)=field_up(1,ip13,ip22,ip33)+&
      &a13a22a33*up(1,iq)
      field_up(2,ip13,ip22,ip33)=field_up(2,ip13,ip22,ip33)+&
      &a13a22a33*up(2,iq)
      field_up(3,ip13,ip22,ip33)=field_up(3,ip13,ip22,ip33)+&
      &a13a22a33*up(3,iq)
      field_up(1,ip10,ip23,ip33)=field_up(1,ip10,ip23,ip33)+&
      &a10a23a33*up(1,iq)
      field_up(2,ip10,ip23,ip33)=field_up(2,ip10,ip23,ip33)+&
      &a10a23a33*up(2,iq)
      field_up(3,ip10,ip23,ip33)=field_up(3,ip10,ip23,ip33)+&
      &a10a23a33*up(3,iq)
      field_up(1,ip11,ip23,ip33)=field_up(1,ip11,ip23,ip33)+&
      &a11a23a33*up(1,iq)
      field_up(2,ip11,ip23,ip33)=field_up(2,ip11,ip23,ip33)+&
      &a11a23a33*up(2,iq)
      field_up(3,ip11,ip23,ip33)=field_up(3,ip11,ip23,ip33)+&
      &a11a23a33*up(3,iq)
      field_up(1,ip12,ip23,ip33)=field_up(1,ip12,ip23,ip33)+&
      &a12a23a33*up(1,iq)
      field_up(2,ip12,ip23,ip33)=field_up(2,ip12,ip23,ip33)+&
      &a12a23a33*up(2,iq)
      field_up(3,ip12,ip23,ip33)=field_up(3,ip12,ip23,ip33)+&
      &a12a23a33*up(3,iq)
      field_up(1,ip13,ip23,ip33)=field_up(1,ip13,ip23,ip33)+&
      &a13a23a33*up(1,iq)
      field_up(2,ip13,ip23,ip33)=field_up(2,ip13,ip23,ip33)+&
      &a13a23a33*up(2,iq)
      field_up(3,ip13,ip23,ip33)=field_up(3,ip13,ip23,ip33)+&
      &a13a23a33*up(3,iq)
#else
      field_up(1:3,ip10,ip20,ip30)=field_up(1:3,ip10,ip20,ip30)+&
      &a10a20a30*up(1:3,iq)
      field_up(1:3,ip11,ip20,ip30)=field_up(1:3,ip11,ip20,ip30)+&
      &a11a20a30*up(1:3,iq)
      field_up(1:3,ip12,ip20,ip30)=field_up(1:3,ip12,ip20,ip30)+&
      &a12a20a30*up(1:3,iq)
      field_up(1:3,ip13,ip20,ip30)=field_up(1:3,ip13,ip20,ip30)+&
      &a13a20a30*up(1:3,iq)
      field_up(1:3,ip10,ip21,ip30)=field_up(1:3,ip10,ip21,ip30)+&
      &a10a21a30*up(1:3,iq)
      field_up(1:3,ip11,ip21,ip30)=field_up(1:3,ip11,ip21,ip30)+&
      &a11a21a30*up(1:3,iq)
      field_up(1:3,ip12,ip21,ip30)=field_up(1:3,ip12,ip21,ip30)+&
      &a12a21a30*up(1:3,iq)
      field_up(1:3,ip13,ip21,ip30)=field_up(1:3,ip13,ip21,ip30)+&
      &a13a21a30*up(1:3,iq)
      field_up(1:3,ip10,ip22,ip30)=field_up(1:3,ip10,ip22,ip30)+&
      &a10a22a30*up(1:3,iq)
      field_up(1:3,ip11,ip22,ip30)=field_up(1:3,ip11,ip22,ip30)+&
      &a11a22a30*up(1:3,iq)
      field_up(1:3,ip12,ip22,ip30)=field_up(1:3,ip12,ip22,ip30)+&
      &a12a22a30*up(1:3,iq)
      field_up(1:3,ip13,ip22,ip30)=field_up(1:3,ip13,ip22,ip30)+&
      &a13a22a30*up(1:3,iq)
      field_up(1:3,ip10,ip23,ip30)=field_up(1:3,ip10,ip23,ip30)+&
      &a10a23a30*up(1:3,iq)
      field_up(1:3,ip11,ip23,ip30)=field_up(1:3,ip11,ip23,ip30)+&
      &a11a23a30*up(1:3,iq)
      field_up(1:3,ip12,ip23,ip30)=field_up(1:3,ip12,ip23,ip30)+&
      &a12a23a30*up(1:3,iq)
      field_up(1:3,ip13,ip23,ip30)=field_up(1:3,ip13,ip23,ip30)+&
      &a13a23a30*up(1:3,iq)
      field_up(1:3,ip10,ip20,ip31)=field_up(1:3,ip10,ip20,ip31)+&
      &a10a20a31*up(1:3,iq)
      field_up(1:3,ip11,ip20,ip31)=field_up(1:3,ip11,ip20,ip31)+&
      &a11a20a31*up(1:3,iq)
      field_up(1:3,ip12,ip20,ip31)=field_up(1:3,ip12,ip20,ip31)+&
      &a12a20a31*up(1:3,iq)
      field_up(1:3,ip13,ip20,ip31)=field_up(1:3,ip13,ip20,ip31)+&
      &a13a20a31*up(1:3,iq)
      field_up(1:3,ip10,ip21,ip31)=field_up(1:3,ip10,ip21,ip31)+&
      &a10a21a31*up(1:3,iq)
      field_up(1:3,ip11,ip21,ip31)=field_up(1:3,ip11,ip21,ip31)+&
      &a11a21a31*up(1:3,iq)
      field_up(1:3,ip12,ip21,ip31)=field_up(1:3,ip12,ip21,ip31)+&
      &a12a21a31*up(1:3,iq)
      field_up(1:3,ip13,ip21,ip31)=field_up(1:3,ip13,ip21,ip31)+&
      &a13a21a31*up(1:3,iq)
      field_up(1:3,ip10,ip22,ip31)=field_up(1:3,ip10,ip22,ip31)+&
      &a10a22a31*up(1:3,iq)
      field_up(1:3,ip11,ip22,ip31)=field_up(1:3,ip11,ip22,ip31)+&
      &a11a22a31*up(1:3,iq)
      field_up(1:3,ip12,ip22,ip31)=field_up(1:3,ip12,ip22,ip31)+&
      &a12a22a31*up(1:3,iq)
      field_up(1:3,ip13,ip22,ip31)=field_up(1:3,ip13,ip22,ip31)+&
      &a13a22a31*up(1:3,iq)
      field_up(1:3,ip10,ip23,ip31)=field_up(1:3,ip10,ip23,ip31)+&
      &a10a23a31*up(1:3,iq)
      field_up(1:3,ip11,ip23,ip31)=field_up(1:3,ip11,ip23,ip31)+&
      &a11a23a31*up(1:3,iq)
      field_up(1:3,ip12,ip23,ip31)=field_up(1:3,ip12,ip23,ip31)+&
      &a12a23a31*up(1:3,iq)
      field_up(1:3,ip13,ip23,ip31)=field_up(1:3,ip13,ip23,ip31)+&
      &a13a23a31*up(1:3,iq)
      field_up(1:3,ip10,ip20,ip32)=field_up(1:3,ip10,ip20,ip32)+&
      &a10a20a32*up(1:3,iq)
      field_up(1:3,ip11,ip20,ip32)=field_up(1:3,ip11,ip20,ip32)+&
      &a11a20a32*up(1:3,iq)
      field_up(1:3,ip12,ip20,ip32)=field_up(1:3,ip12,ip20,ip32)+&
      &a12a20a32*up(1:3,iq)
      field_up(1:3,ip13,ip20,ip32)=field_up(1:3,ip13,ip20,ip32)+&
      &a13a20a32*up(1:3,iq)
      field_up(1:3,ip10,ip21,ip32)=field_up(1:3,ip10,ip21,ip32)+&
      &a10a21a32*up(1:3,iq)
      field_up(1:3,ip11,ip21,ip32)=field_up(1:3,ip11,ip21,ip32)+&
      &a11a21a32*up(1:3,iq)
      field_up(1:3,ip12,ip21,ip32)=field_up(1:3,ip12,ip21,ip32)+&
      &a12a21a32*up(1:3,iq)
      field_up(1:3,ip13,ip21,ip32)=field_up(1:3,ip13,ip21,ip32)+&
      &a13a21a32*up(1:3,iq)
      field_up(1:3,ip10,ip22,ip32)=field_up(1:3,ip10,ip22,ip32)+&
      &a10a22a32*up(1:3,iq)
      field_up(1:3,ip11,ip22,ip32)=field_up(1:3,ip11,ip22,ip32)+&
      &a11a22a32*up(1:3,iq)
      field_up(1:3,ip12,ip22,ip32)=field_up(1:3,ip12,ip22,ip32)+&
      &a12a22a32*up(1:3,iq)
      field_up(1:3,ip13,ip22,ip32)=field_up(1:3,ip13,ip22,ip32)+&
      &a13a22a32*up(1:3,iq)
      field_up(1:3,ip10,ip23,ip32)=field_up(1:3,ip10,ip23,ip32)+&
      &a10a23a32*up(1:3,iq)
      field_up(1:3,ip11,ip23,ip32)=field_up(1:3,ip11,ip23,ip32)+&
      &a11a23a32*up(1:3,iq)
      field_up(1:3,ip12,ip23,ip32)=field_up(1:3,ip12,ip23,ip32)+&
      &a12a23a32*up(1:3,iq)
      field_up(1:3,ip13,ip23,ip32)=field_up(1:3,ip13,ip23,ip32)+&
      &a13a23a32*up(1:3,iq)
      field_up(1:3,ip10,ip20,ip33)=field_up(1:3,ip10,ip20,ip33)+&
      &a10a20a33*up(1:3,iq)
      field_up(1:3,ip11,ip20,ip33)=field_up(1:3,ip11,ip20,ip33)+&
      &a11a20a33*up(1:3,iq)
      field_up(1:3,ip12,ip20,ip33)=field_up(1:3,ip12,ip20,ip33)+&
      &a12a20a33*up(1:3,iq)
      field_up(1:3,ip13,ip20,ip33)=field_up(1:3,ip13,ip20,ip33)+&
      &a13a20a33*up(1:3,iq)
      field_up(1:3,ip10,ip21,ip33)=field_up(1:3,ip10,ip21,ip33)+&
      &a10a21a33*up(1:3,iq)
      field_up(1:3,ip11,ip21,ip33)=field_up(1:3,ip11,ip21,ip33)+&
      &a11a21a33*up(1:3,iq)
      field_up(1:3,ip12,ip21,ip33)=field_up(1:3,ip12,ip21,ip33)+&
      &a12a21a33*up(1:3,iq)
      field_up(1:3,ip13,ip21,ip33)=field_up(1:3,ip13,ip21,ip33)+&
      &a13a21a33*up(1:3,iq)
      field_up(1:3,ip10,ip22,ip33)=field_up(1:3,ip10,ip22,ip33)+&
      &a10a22a33*up(1:3,iq)
      field_up(1:3,ip11,ip22,ip33)=field_up(1:3,ip11,ip22,ip33)+&
      &a11a22a33*up(1:3,iq)
      field_up(1:3,ip12,ip22,ip33)=field_up(1:3,ip12,ip22,ip33)+&
      &a12a22a33*up(1:3,iq)
      field_up(1:3,ip13,ip22,ip33)=field_up(1:3,ip13,ip22,ip33)+&
      &a13a22a33*up(1:3,iq)
      field_up(1:3,ip10,ip23,ip33)=field_up(1:3,ip10,ip23,ip33)+&
      &a10a23a33*up(1:3,iq)
      field_up(1:3,ip11,ip23,ip33)=field_up(1:3,ip11,ip23,ip33)+&
      &a11a23a33*up(1:3,iq)
      field_up(1:3,ip12,ip23,ip33)=field_up(1:3,ip12,ip23,ip33)+&
      &a12a23a33*up(1:3,iq)
      field_up(1:3,ip13,ip23,ip33)=field_up(1:3,ip13,ip23,ip33)+&
      &a13a23a33*up(1:3,iq)

#endif
              END DO
              !----------------------------------------------------------------!
              !  Unrolled versions for 2-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 2) THEN
               DO ip = 1,store_info(ipatch)
                   iq    = list_sub(ipatch,ip)

                   x01 = xp(1,iq)*dxxi-p%istart(1) + 1
                   x02 = xp(2,iq)*dxyi-p%istart(2) + 1
                   x03 = xp(3,iq)*dxzi-p%istart(3) + 1

                   ip10 = FLOOR(x01)
                   ip20 = FLOOR(x02)
                   ip30 = FLOOR(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

#ifdef __NOMICROINSTRUCTIONS
      field_up(1,ip10,ip20,ip30)=field_up(1,ip10,ip20,ip30)+&
      &a10a20a30*up(1,iq)
      field_up(2,ip10,ip20,ip30)=field_up(2,ip10,ip20,ip30)+&
      &a10a20a30*up(2,iq)
      field_up(1,ip11,ip20,ip30)=field_up(1,ip11,ip20,ip30)+&
      &a11a20a30*up(1,iq)
      field_up(2,ip11,ip20,ip30)=field_up(2,ip11,ip20,ip30)+&
      &a11a20a30*up(2,iq)
      field_up(1,ip12,ip20,ip30)=field_up(1,ip12,ip20,ip30)+&
      &a12a20a30*up(1,iq)
      field_up(2,ip12,ip20,ip30)=field_up(2,ip12,ip20,ip30)+&
      &a12a20a30*up(2,iq)
      field_up(1,ip13,ip20,ip30)=field_up(1,ip13,ip20,ip30)+&
      &a13a20a30*up(1,iq)
      field_up(2,ip13,ip20,ip30)=field_up(2,ip13,ip20,ip30)+&
      &a13a20a30*up(2,iq)
      field_up(1,ip10,ip21,ip30)=field_up(1,ip10,ip21,ip30)+&
      &a10a21a30*up(1,iq)
      field_up(2,ip10,ip21,ip30)=field_up(2,ip10,ip21,ip30)+&
      &a10a21a30*up(2,iq)
      field_up(1,ip11,ip21,ip30)=field_up(1,ip11,ip21,ip30)+&
      &a11a21a30*up(1,iq)
      field_up(2,ip11,ip21,ip30)=field_up(2,ip11,ip21,ip30)+&
      &a11a21a30*up(2,iq)
      field_up(1,ip12,ip21,ip30)=field_up(1,ip12,ip21,ip30)+&
      &a12a21a30*up(1,iq)
      field_up(2,ip12,ip21,ip30)=field_up(2,ip12,ip21,ip30)+&
      &a12a21a30*up(2,iq)
      field_up(1,ip13,ip21,ip30)=field_up(1,ip13,ip21,ip30)+&
      &a13a21a30*up(1,iq)
      field_up(2,ip13,ip21,ip30)=field_up(2,ip13,ip21,ip30)+&
      &a13a21a30*up(2,iq)
      field_up(1,ip10,ip22,ip30)=field_up(1,ip10,ip22,ip30)+&
      &a10a22a30*up(1,iq)
      field_up(2,ip10,ip22,ip30)=field_up(2,ip10,ip22,ip30)+&
      &a10a22a30*up(2,iq)
      field_up(1,ip11,ip22,ip30)=field_up(1,ip11,ip22,ip30)+&
      &a11a22a30*up(1,iq)
      field_up(2,ip11,ip22,ip30)=field_up(2,ip11,ip22,ip30)+&
      &a11a22a30*up(2,iq)
      field_up(1,ip12,ip22,ip30)=field_up(1,ip12,ip22,ip30)+&
      &a12a22a30*up(1,iq)
      field_up(2,ip12,ip22,ip30)=field_up(2,ip12,ip22,ip30)+&
      &a12a22a30*up(2,iq)
      field_up(1,ip13,ip22,ip30)=field_up(1,ip13,ip22,ip30)+&
      &a13a22a30*up(1,iq)
      field_up(2,ip13,ip22,ip30)=field_up(2,ip13,ip22,ip30)+&
      &a13a22a30*up(2,iq)
      field_up(1,ip10,ip23,ip30)=field_up(1,ip10,ip23,ip30)+&
      &a10a23a30*up(1,iq)
      field_up(2,ip10,ip23,ip30)=field_up(2,ip10,ip23,ip30)+&
      &a10a23a30*up(2,iq)
      field_up(1,ip11,ip23,ip30)=field_up(1,ip11,ip23,ip30)+&
      &a11a23a30*up(1,iq)
      field_up(2,ip11,ip23,ip30)=field_up(2,ip11,ip23,ip30)+&
      &a11a23a30*up(2,iq)
      field_up(1,ip12,ip23,ip30)=field_up(1,ip12,ip23,ip30)+&
      &a12a23a30*up(1,iq)
      field_up(2,ip12,ip23,ip30)=field_up(2,ip12,ip23,ip30)+&
      &a12a23a30*up(2,iq)
      field_up(1,ip13,ip23,ip30)=field_up(1,ip13,ip23,ip30)+&
      &a13a23a30*up(1,iq)
      field_up(2,ip13,ip23,ip30)=field_up(2,ip13,ip23,ip30)+&
      &a13a23a30*up(2,iq)
      field_up(1,ip10,ip20,ip31)=field_up(1,ip10,ip20,ip31)+&
      &a10a20a31*up(1,iq)
      field_up(2,ip10,ip20,ip31)=field_up(2,ip10,ip20,ip31)+&
      &a10a20a31*up(2,iq)
      field_up(1,ip11,ip20,ip31)=field_up(1,ip11,ip20,ip31)+&
      &a11a20a31*up(1,iq)
      field_up(2,ip11,ip20,ip31)=field_up(2,ip11,ip20,ip31)+&
      &a11a20a31*up(2,iq)
      field_up(1,ip12,ip20,ip31)=field_up(1,ip12,ip20,ip31)+&
      &a12a20a31*up(1,iq)
      field_up(2,ip12,ip20,ip31)=field_up(2,ip12,ip20,ip31)+&
      &a12a20a31*up(2,iq)
      field_up(1,ip13,ip20,ip31)=field_up(1,ip13,ip20,ip31)+&
      &a13a20a31*up(1,iq)
      field_up(2,ip13,ip20,ip31)=field_up(2,ip13,ip20,ip31)+&
      &a13a20a31*up(2,iq)
      field_up(1,ip10,ip21,ip31)=field_up(1,ip10,ip21,ip31)+&
      &a10a21a31*up(1,iq)
      field_up(2,ip10,ip21,ip31)=field_up(2,ip10,ip21,ip31)+&
      &a10a21a31*up(2,iq)
      field_up(1,ip11,ip21,ip31)=field_up(1,ip11,ip21,ip31)+&
      &a11a21a31*up(1,iq)
      field_up(2,ip11,ip21,ip31)=field_up(2,ip11,ip21,ip31)+&
      &a11a21a31*up(2,iq)
      field_up(1,ip12,ip21,ip31)=field_up(1,ip12,ip21,ip31)+&
      &a12a21a31*up(1,iq)
      field_up(2,ip12,ip21,ip31)=field_up(2,ip12,ip21,ip31)+&
      &a12a21a31*up(2,iq)
      field_up(1,ip13,ip21,ip31)=field_up(1,ip13,ip21,ip31)+&
      &a13a21a31*up(1,iq)
      field_up(2,ip13,ip21,ip31)=field_up(2,ip13,ip21,ip31)+&
      &a13a21a31*up(2,iq)
      field_up(1,ip10,ip22,ip31)=field_up(1,ip10,ip22,ip31)+&
      &a10a22a31*up(1,iq)
      field_up(2,ip10,ip22,ip31)=field_up(2,ip10,ip22,ip31)+&
      &a10a22a31*up(2,iq)
      field_up(1,ip11,ip22,ip31)=field_up(1,ip11,ip22,ip31)+&
      &a11a22a31*up(1,iq)
      field_up(2,ip11,ip22,ip31)=field_up(2,ip11,ip22,ip31)+&
      &a11a22a31*up(2,iq)
      field_up(1,ip12,ip22,ip31)=field_up(1,ip12,ip22,ip31)+&
      &a12a22a31*up(1,iq)
      field_up(2,ip12,ip22,ip31)=field_up(2,ip12,ip22,ip31)+&
      &a12a22a31*up(2,iq)
      field_up(1,ip13,ip22,ip31)=field_up(1,ip13,ip22,ip31)+&
      &a13a22a31*up(1,iq)
      field_up(2,ip13,ip22,ip31)=field_up(2,ip13,ip22,ip31)+&
      &a13a22a31*up(2,iq)
      field_up(1,ip10,ip23,ip31)=field_up(1,ip10,ip23,ip31)+&
      &a10a23a31*up(1,iq)
      field_up(2,ip10,ip23,ip31)=field_up(2,ip10,ip23,ip31)+&
      &a10a23a31*up(2,iq)
      field_up(1,ip11,ip23,ip31)=field_up(1,ip11,ip23,ip31)+&
      &a11a23a31*up(1,iq)
      field_up(2,ip11,ip23,ip31)=field_up(2,ip11,ip23,ip31)+&
      &a11a23a31*up(2,iq)
      field_up(1,ip12,ip23,ip31)=field_up(1,ip12,ip23,ip31)+&
      &a12a23a31*up(1,iq)
      field_up(2,ip12,ip23,ip31)=field_up(2,ip12,ip23,ip31)+&
      &a12a23a31*up(2,iq)
      field_up(1,ip13,ip23,ip31)=field_up(1,ip13,ip23,ip31)+&
      &a13a23a31*up(1,iq)
      field_up(2,ip13,ip23,ip31)=field_up(2,ip13,ip23,ip31)+&
      &a13a23a31*up(2,iq)
      field_up(1,ip10,ip20,ip32)=field_up(1,ip10,ip20,ip32)+&
      &a10a20a32*up(1,iq)
      field_up(2,ip10,ip20,ip32)=field_up(2,ip10,ip20,ip32)+&
      &a10a20a32*up(2,iq)
      field_up(1,ip11,ip20,ip32)=field_up(1,ip11,ip20,ip32)+&
      &a11a20a32*up(1,iq)
      field_up(2,ip11,ip20,ip32)=field_up(2,ip11,ip20,ip32)+&
      &a11a20a32*up(2,iq)
      field_up(1,ip12,ip20,ip32)=field_up(1,ip12,ip20,ip32)+&
      &a12a20a32*up(1,iq)
      field_up(2,ip12,ip20,ip32)=field_up(2,ip12,ip20,ip32)+&
      &a12a20a32*up(2,iq)
      field_up(1,ip13,ip20,ip32)=field_up(1,ip13,ip20,ip32)+&
      &a13a20a32*up(1,iq)
      field_up(2,ip13,ip20,ip32)=field_up(2,ip13,ip20,ip32)+&
      &a13a20a32*up(2,iq)
      field_up(1,ip10,ip21,ip32)=field_up(1,ip10,ip21,ip32)+&
      &a10a21a32*up(1,iq)
      field_up(2,ip10,ip21,ip32)=field_up(2,ip10,ip21,ip32)+&
      &a10a21a32*up(2,iq)
      field_up(1,ip11,ip21,ip32)=field_up(1,ip11,ip21,ip32)+&
      &a11a21a32*up(1,iq)
      field_up(2,ip11,ip21,ip32)=field_up(2,ip11,ip21,ip32)+&
      &a11a21a32*up(2,iq)
      field_up(1,ip12,ip21,ip32)=field_up(1,ip12,ip21,ip32)+&
      &a12a21a32*up(1,iq)
      field_up(2,ip12,ip21,ip32)=field_up(2,ip12,ip21,ip32)+&
      &a12a21a32*up(2,iq)
      field_up(1,ip13,ip21,ip32)=field_up(1,ip13,ip21,ip32)+&
      &a13a21a32*up(1,iq)
      field_up(2,ip13,ip21,ip32)=field_up(2,ip13,ip21,ip32)+&
      &a13a21a32*up(2,iq)
      field_up(1,ip10,ip22,ip32)=field_up(1,ip10,ip22,ip32)+&
      &a10a22a32*up(1,iq)
      field_up(2,ip10,ip22,ip32)=field_up(2,ip10,ip22,ip32)+&
      &a10a22a32*up(2,iq)
      field_up(1,ip11,ip22,ip32)=field_up(1,ip11,ip22,ip32)+&
      &a11a22a32*up(1,iq)
      field_up(2,ip11,ip22,ip32)=field_up(2,ip11,ip22,ip32)+&
      &a11a22a32*up(2,iq)
      field_up(1,ip12,ip22,ip32)=field_up(1,ip12,ip22,ip32)+&
      &a12a22a32*up(1,iq)
      field_up(2,ip12,ip22,ip32)=field_up(2,ip12,ip22,ip32)+&
      &a12a22a32*up(2,iq)
      field_up(1,ip13,ip22,ip32)=field_up(1,ip13,ip22,ip32)+&
      &a13a22a32*up(1,iq)
      field_up(2,ip13,ip22,ip32)=field_up(2,ip13,ip22,ip32)+&
      &a13a22a32*up(2,iq)
      field_up(1,ip10,ip23,ip32)=field_up(1,ip10,ip23,ip32)+&
      &a10a23a32*up(1,iq)
      field_up(2,ip10,ip23,ip32)=field_up(2,ip10,ip23,ip32)+&
      &a10a23a32*up(2,iq)
      field_up(1,ip11,ip23,ip32)=field_up(1,ip11,ip23,ip32)+&
      &a11a23a32*up(1,iq)
      field_up(2,ip11,ip23,ip32)=field_up(2,ip11,ip23,ip32)+&
      &a11a23a32*up(2,iq)
      field_up(1,ip12,ip23,ip32)=field_up(1,ip12,ip23,ip32)+&
      &a12a23a32*up(1,iq)
      field_up(2,ip12,ip23,ip32)=field_up(2,ip12,ip23,ip32)+&
      &a12a23a32*up(2,iq)
      field_up(1,ip13,ip23,ip32)=field_up(1,ip13,ip23,ip32)+&
      &a13a23a32*up(1,iq)
      field_up(2,ip13,ip23,ip32)=field_up(2,ip13,ip23,ip32)+&
      &a13a23a32*up(2,iq)
      field_up(1,ip10,ip20,ip33)=field_up(1,ip10,ip20,ip33)+&
      &a10a20a33*up(1,iq)
      field_up(2,ip10,ip20,ip33)=field_up(2,ip10,ip20,ip33)+&
      &a10a20a33*up(2,iq)
      field_up(1,ip11,ip20,ip33)=field_up(1,ip11,ip20,ip33)+&
      &a11a20a33*up(1,iq)
      field_up(2,ip11,ip20,ip33)=field_up(2,ip11,ip20,ip33)+&
      &a11a20a33*up(2,iq)
      field_up(1,ip12,ip20,ip33)=field_up(1,ip12,ip20,ip33)+&
      &a12a20a33*up(1,iq)
      field_up(2,ip12,ip20,ip33)=field_up(2,ip12,ip20,ip33)+&
      &a12a20a33*up(2,iq)
      field_up(1,ip13,ip20,ip33)=field_up(1,ip13,ip20,ip33)+&
      &a13a20a33*up(1,iq)
      field_up(2,ip13,ip20,ip33)=field_up(2,ip13,ip20,ip33)+&
      &a13a20a33*up(2,iq)
      field_up(1,ip10,ip21,ip33)=field_up(1,ip10,ip21,ip33)+&
      &a10a21a33*up(1,iq)
      field_up(2,ip10,ip21,ip33)=field_up(2,ip10,ip21,ip33)+&
      &a10a21a33*up(2,iq)
      field_up(1,ip11,ip21,ip33)=field_up(1,ip11,ip21,ip33)+&
      &a11a21a33*up(1,iq)
      field_up(2,ip11,ip21,ip33)=field_up(2,ip11,ip21,ip33)+&
      &a11a21a33*up(2,iq)
      field_up(1,ip12,ip21,ip33)=field_up(1,ip12,ip21,ip33)+&
      &a12a21a33*up(1,iq)
      field_up(2,ip12,ip21,ip33)=field_up(2,ip12,ip21,ip33)+&
      &a12a21a33*up(2,iq)
      field_up(1,ip13,ip21,ip33)=field_up(1,ip13,ip21,ip33)+&
      &a13a21a33*up(1,iq)
      field_up(2,ip13,ip21,ip33)=field_up(2,ip13,ip21,ip33)+&
      &a13a21a33*up(2,iq)
      field_up(1,ip10,ip22,ip33)=field_up(1,ip10,ip22,ip33)+&
      &a10a22a33*up(1,iq)
      field_up(2,ip10,ip22,ip33)=field_up(2,ip10,ip22,ip33)+&
      &a10a22a33*up(2,iq)
      field_up(1,ip11,ip22,ip33)=field_up(1,ip11,ip22,ip33)+&
      &a11a22a33*up(1,iq)
      field_up(2,ip11,ip22,ip33)=field_up(2,ip11,ip22,ip33)+&
      &a11a22a33*up(2,iq)
      field_up(1,ip12,ip22,ip33)=field_up(1,ip12,ip22,ip33)+&
      &a12a22a33*up(1,iq)
      field_up(2,ip12,ip22,ip33)=field_up(2,ip12,ip22,ip33)+&
      &a12a22a33*up(2,iq)
      field_up(1,ip13,ip22,ip33)=field_up(1,ip13,ip22,ip33)+&
      &a13a22a33*up(1,iq)
      field_up(2,ip13,ip22,ip33)=field_up(2,ip13,ip22,ip33)+&
      &a13a22a33*up(2,iq)
      field_up(1,ip10,ip23,ip33)=field_up(1,ip10,ip23,ip33)+&
      &a10a23a33*up(1,iq)
      field_up(2,ip10,ip23,ip33)=field_up(2,ip10,ip23,ip33)+&
      &a10a23a33*up(2,iq)
      field_up(1,ip11,ip23,ip33)=field_up(1,ip11,ip23,ip33)+&
      &a11a23a33*up(1,iq)
      field_up(2,ip11,ip23,ip33)=field_up(2,ip11,ip23,ip33)+&
      &a11a23a33*up(2,iq)
      field_up(1,ip12,ip23,ip33)=field_up(1,ip12,ip23,ip33)+&
      &a12a23a33*up(1,iq)
      field_up(2,ip12,ip23,ip33)=field_up(2,ip12,ip23,ip33)+&
      &a12a23a33*up(2,iq)
      field_up(1,ip13,ip23,ip33)=field_up(1,ip13,ip23,ip33)+&
      &a13a23a33*up(1,iq)
      field_up(2,ip13,ip23,ip33)=field_up(2,ip13,ip23,ip33)+&
      &a13a23a33*up(2,iq)
#else
      field_up(1:2,ip10,ip20,ip30)=field_up(1:2,ip10,ip20,ip30)+&
      &a10a20a30*up(1:2,iq)
      field_up(1:2,ip11,ip20,ip30)=field_up(1:2,ip11,ip20,ip30)+&
      &a11a20a30*up(1:2,iq)
      field_up(1:2,ip12,ip20,ip30)=field_up(1:2,ip12,ip20,ip30)+&
      &a12a20a30*up(1:2,iq)
      field_up(1:2,ip13,ip20,ip30)=field_up(1:2,ip13,ip20,ip30)+&
      &a13a20a30*up(1:2,iq)
      field_up(1:2,ip10,ip21,ip30)=field_up(1:2,ip10,ip21,ip30)+&
      &a10a21a30*up(1:2,iq)
      field_up(1:2,ip11,ip21,ip30)=field_up(1:2,ip11,ip21,ip30)+&
      &a11a21a30*up(1:2,iq)
      field_up(1:2,ip12,ip21,ip30)=field_up(1:2,ip12,ip21,ip30)+&
      &a12a21a30*up(1:2,iq)
      field_up(1:2,ip13,ip21,ip30)=field_up(1:2,ip13,ip21,ip30)+&
      &a13a21a30*up(1:2,iq)
      field_up(1:2,ip10,ip22,ip30)=field_up(1:2,ip10,ip22,ip30)+&
      &a10a22a30*up(1:2,iq)
      field_up(1:2,ip11,ip22,ip30)=field_up(1:2,ip11,ip22,ip30)+&
      &a11a22a30*up(1:2,iq)
      field_up(1:2,ip12,ip22,ip30)=field_up(1:2,ip12,ip22,ip30)+&
      &a12a22a30*up(1:2,iq)
      field_up(1:2,ip13,ip22,ip30)=field_up(1:2,ip13,ip22,ip30)+&
      &a13a22a30*up(1:2,iq)
      field_up(1:2,ip10,ip23,ip30)=field_up(1:2,ip10,ip23,ip30)+&
      &a10a23a30*up(1:2,iq)
      field_up(1:2,ip11,ip23,ip30)=field_up(1:2,ip11,ip23,ip30)+&
      &a11a23a30*up(1:2,iq)
      field_up(1:2,ip12,ip23,ip30)=field_up(1:2,ip12,ip23,ip30)+&
      &a12a23a30*up(1:2,iq)
      field_up(1:2,ip13,ip23,ip30)=field_up(1:2,ip13,ip23,ip30)+&
      &a13a23a30*up(1:2,iq)
      field_up(1:2,ip10,ip20,ip31)=field_up(1:2,ip10,ip20,ip31)+&
      &a10a20a31*up(1:2,iq)
      field_up(1:2,ip11,ip20,ip31)=field_up(1:2,ip11,ip20,ip31)+&
      &a11a20a31*up(1:2,iq)
      field_up(1:2,ip12,ip20,ip31)=field_up(1:2,ip12,ip20,ip31)+&
      &a12a20a31*up(1:2,iq)
      field_up(1:2,ip13,ip20,ip31)=field_up(1:2,ip13,ip20,ip31)+&
      &a13a20a31*up(1:2,iq)
      field_up(1:2,ip10,ip21,ip31)=field_up(1:2,ip10,ip21,ip31)+&
      &a10a21a31*up(1:2,iq)
      field_up(1:2,ip11,ip21,ip31)=field_up(1:2,ip11,ip21,ip31)+&
      &a11a21a31*up(1:2,iq)
      field_up(1:2,ip12,ip21,ip31)=field_up(1:2,ip12,ip21,ip31)+&
      &a12a21a31*up(1:2,iq)
      field_up(1:2,ip13,ip21,ip31)=field_up(1:2,ip13,ip21,ip31)+&
      &a13a21a31*up(1:2,iq)
      field_up(1:2,ip10,ip22,ip31)=field_up(1:2,ip10,ip22,ip31)+&
      &a10a22a31*up(1:2,iq)
      field_up(1:2,ip11,ip22,ip31)=field_up(1:2,ip11,ip22,ip31)+&
      &a11a22a31*up(1:2,iq)
      field_up(1:2,ip12,ip22,ip31)=field_up(1:2,ip12,ip22,ip31)+&
      &a12a22a31*up(1:2,iq)
      field_up(1:2,ip13,ip22,ip31)=field_up(1:2,ip13,ip22,ip31)+&
      &a13a22a31*up(1:2,iq)
      field_up(1:2,ip10,ip23,ip31)=field_up(1:2,ip10,ip23,ip31)+&
      &a10a23a31*up(1:2,iq)
      field_up(1:2,ip11,ip23,ip31)=field_up(1:2,ip11,ip23,ip31)+&
      &a11a23a31*up(1:2,iq)
      field_up(1:2,ip12,ip23,ip31)=field_up(1:2,ip12,ip23,ip31)+&
      &a12a23a31*up(1:2,iq)
      field_up(1:2,ip13,ip23,ip31)=field_up(1:2,ip13,ip23,ip31)+&
      &a13a23a31*up(1:2,iq)
      field_up(1:2,ip10,ip20,ip32)=field_up(1:2,ip10,ip20,ip32)+&
      &a10a20a32*up(1:2,iq)
      field_up(1:2,ip11,ip20,ip32)=field_up(1:2,ip11,ip20,ip32)+&
      &a11a20a32*up(1:2,iq)
      field_up(1:2,ip12,ip20,ip32)=field_up(1:2,ip12,ip20,ip32)+&
      &a12a20a32*up(1:2,iq)
      field_up(1:2,ip13,ip20,ip32)=field_up(1:2,ip13,ip20,ip32)+&
      &a13a20a32*up(1:2,iq)
      field_up(1:2,ip10,ip21,ip32)=field_up(1:2,ip10,ip21,ip32)+&
      &a10a21a32*up(1:2,iq)
      field_up(1:2,ip11,ip21,ip32)=field_up(1:2,ip11,ip21,ip32)+&
      &a11a21a32*up(1:2,iq)
      field_up(1:2,ip12,ip21,ip32)=field_up(1:2,ip12,ip21,ip32)+&
      &a12a21a32*up(1:2,iq)
      field_up(1:2,ip13,ip21,ip32)=field_up(1:2,ip13,ip21,ip32)+&
      &a13a21a32*up(1:2,iq)
      field_up(1:2,ip10,ip22,ip32)=field_up(1:2,ip10,ip22,ip32)+&
      &a10a22a32*up(1:2,iq)
      field_up(1:2,ip11,ip22,ip32)=field_up(1:2,ip11,ip22,ip32)+&
      &a11a22a32*up(1:2,iq)
      field_up(1:2,ip12,ip22,ip32)=field_up(1:2,ip12,ip22,ip32)+&
      &a12a22a32*up(1:2,iq)
      field_up(1:2,ip13,ip22,ip32)=field_up(1:2,ip13,ip22,ip32)+&
      &a13a22a32*up(1:2,iq)
      field_up(1:2,ip10,ip23,ip32)=field_up(1:2,ip10,ip23,ip32)+&
      &a10a23a32*up(1:2,iq)
      field_up(1:2,ip11,ip23,ip32)=field_up(1:2,ip11,ip23,ip32)+&
      &a11a23a32*up(1:2,iq)
      field_up(1:2,ip12,ip23,ip32)=field_up(1:2,ip12,ip23,ip32)+&
      &a12a23a32*up(1:2,iq)
      field_up(1:2,ip13,ip23,ip32)=field_up(1:2,ip13,ip23,ip32)+&
      &a13a23a32*up(1:2,iq)
      field_up(1:2,ip10,ip20,ip33)=field_up(1:2,ip10,ip20,ip33)+&
      &a10a20a33*up(1:2,iq)
      field_up(1:2,ip11,ip20,ip33)=field_up(1:2,ip11,ip20,ip33)+&
      &a11a20a33*up(1:2,iq)
      field_up(1:2,ip12,ip20,ip33)=field_up(1:2,ip12,ip20,ip33)+&
      &a12a20a33*up(1:2,iq)
      field_up(1:2,ip13,ip20,ip33)=field_up(1:2,ip13,ip20,ip33)+&
      &a13a20a33*up(1:2,iq)
      field_up(1:2,ip10,ip21,ip33)=field_up(1:2,ip10,ip21,ip33)+&
      &a10a21a33*up(1:2,iq)
      field_up(1:2,ip11,ip21,ip33)=field_up(1:2,ip11,ip21,ip33)+&
      &a11a21a33*up(1:2,iq)
      field_up(1:2,ip12,ip21,ip33)=field_up(1:2,ip12,ip21,ip33)+&
      &a12a21a33*up(1:2,iq)
      field_up(1:2,ip13,ip21,ip33)=field_up(1:2,ip13,ip21,ip33)+&
      &a13a21a33*up(1:2,iq)
      field_up(1:2,ip10,ip22,ip33)=field_up(1:2,ip10,ip22,ip33)+&
      &a10a22a33*up(1:2,iq)
      field_up(1:2,ip11,ip22,ip33)=field_up(1:2,ip11,ip22,ip33)+&
      &a11a22a33*up(1:2,iq)
      field_up(1:2,ip12,ip22,ip33)=field_up(1:2,ip12,ip22,ip33)+&
      &a12a22a33*up(1:2,iq)
      field_up(1:2,ip13,ip22,ip33)=field_up(1:2,ip13,ip22,ip33)+&
      &a13a22a33*up(1:2,iq)
      field_up(1:2,ip10,ip23,ip33)=field_up(1:2,ip10,ip23,ip33)+&
      &a10a23a33*up(1:2,iq)
      field_up(1:2,ip11,ip23,ip33)=field_up(1:2,ip11,ip23,ip33)+&
      &a11a23a33*up(1:2,iq)
      field_up(1:2,ip12,ip23,ip33)=field_up(1:2,ip12,ip23,ip33)+&
      &a12a23a33*up(1:2,iq)
      field_up(1:2,ip13,ip23,ip33)=field_up(1:2,ip13,ip23,ip33)+&
      &a13a23a33*up(1:2,iq)

#endif
              ENDDO
              !----------------------------------------------------------------!
              !  Unrolled versions for 1-vectors
              !----------------------------------------------------------------!
           ELSEIF (lda .EQ. 1) THEN
               DO ip = 1,store_info(ipatch)
                   iq    = list_sub(ipatch,ip)

                 x01 = xp(1,iq)*dxxi-p%istart(1) + 1
                 x02 = xp(2,iq)*dxyi-p%istart(2) + 1
                 x03 = xp(3,iq)*dxzi-p%istart(3) + 1

                 ip10 = FLOOR(x01)
                 ip20 = FLOOR(x02)
                 ip30 = FLOOR(x03)

                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33


                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

             field_up(1,ip10,ip20,ip30)=field_up(1,ip10,ip20,ip30)+&
     &                  a10a20a30*up(1,iq)
             field_up(1,ip11,ip20,ip30)=field_up(1,ip11,ip20,ip30)+&
     &                  a11a20a30*up(1,iq)
             field_up(1,ip12,ip20,ip30)=field_up(1,ip12,ip20,ip30)+&
     &                  a12a20a30*up(1,iq)
             field_up(1,ip13,ip20,ip30)=field_up(1,ip13,ip20,ip30)+&
     &                  a13a20a30*up(1,iq)
             field_up(1,ip10,ip21,ip30)=field_up(1,ip10,ip21,ip30)+&
     &                  a10a21a30*up(1,iq)
             field_up(1,ip11,ip21,ip30)=field_up(1,ip11,ip21,ip30)+&
     &                  a11a21a30*up(1,iq)
             field_up(1,ip12,ip21,ip30)=field_up(1,ip12,ip21,ip30)+&
     &                  a12a21a30*up(1,iq)
             field_up(1,ip13,ip21,ip30)=field_up(1,ip13,ip21,ip30)+&
     &                  a13a21a30*up(1,iq)
             field_up(1,ip10,ip22,ip30)=field_up(1,ip10,ip22,ip30)+&
     &                  a10a22a30*up(1,iq)
             field_up(1,ip11,ip22,ip30)=field_up(1,ip11,ip22,ip30)+&
     &                  a11a22a30*up(1,iq)
             field_up(1,ip12,ip22,ip30)=field_up(1,ip12,ip22,ip30)+&
     &                  a12a22a30*up(1,iq)
             field_up(1,ip13,ip22,ip30)=field_up(1,ip13,ip22,ip30)+&
     &                  a13a22a30*up(1,iq)
             field_up(1,ip10,ip23,ip30)=field_up(1,ip10,ip23,ip30)+&
     &                  a10a23a30*up(1,iq)
             field_up(1,ip11,ip23,ip30)=field_up(1,ip11,ip23,ip30)+&
     &                  a11a23a30*up(1,iq)
             field_up(1,ip12,ip23,ip30)=field_up(1,ip12,ip23,ip30)+&
     &                  a12a23a30*up(1,iq)
             field_up(1,ip13,ip23,ip30)=field_up(1,ip13,ip23,ip30)+&
     &                  a13a23a30*up(1,iq)
             field_up(1,ip10,ip20,ip31)=field_up(1,ip10,ip20,ip31)+&
     &                  a10a20a31*up(1,iq)
             field_up(1,ip11,ip20,ip31)=field_up(1,ip11,ip20,ip31)+&
     &                  a11a20a31*up(1,iq)
             field_up(1,ip12,ip20,ip31)=field_up(1,ip12,ip20,ip31)+&
     &                  a12a20a31*up(1,iq)
             field_up(1,ip13,ip20,ip31)=field_up(1,ip13,ip20,ip31)+&
     &                  a13a20a31*up(1,iq)
             field_up(1,ip10,ip21,ip31)=field_up(1,ip10,ip21,ip31)+&
     &                  a10a21a31*up(1,iq)
             field_up(1,ip11,ip21,ip31)=field_up(1,ip11,ip21,ip31)+&
     &                  a11a21a31*up(1,iq)
             field_up(1,ip12,ip21,ip31)=field_up(1,ip12,ip21,ip31)+&
     &                  a12a21a31*up(1,iq)
             field_up(1,ip13,ip21,ip31)=field_up(1,ip13,ip21,ip31)+&
     &                  a13a21a31*up(1,iq)
             field_up(1,ip10,ip22,ip31)=field_up(1,ip10,ip22,ip31)+&
     &                  a10a22a31*up(1,iq)
             field_up(1,ip11,ip22,ip31)=field_up(1,ip11,ip22,ip31)+&
     &                  a11a22a31*up(1,iq)
             field_up(1,ip12,ip22,ip31)=field_up(1,ip12,ip22,ip31)+&
     &                  a12a22a31*up(1,iq)
             field_up(1,ip13,ip22,ip31)=field_up(1,ip13,ip22,ip31)+&
     &                  a13a22a31*up(1,iq)
             field_up(1,ip10,ip23,ip31)=field_up(1,ip10,ip23,ip31)+&
     &                  a10a23a31*up(1,iq)
             field_up(1,ip11,ip23,ip31)=field_up(1,ip11,ip23,ip31)+&
     &                  a11a23a31*up(1,iq)
             field_up(1,ip12,ip23,ip31)=field_up(1,ip12,ip23,ip31)+&
     &                  a12a23a31*up(1,iq)
             field_up(1,ip13,ip23,ip31)=field_up(1,ip13,ip23,ip31)+&
     &                  a13a23a31*up(1,iq)
             field_up(1,ip10,ip20,ip32)=field_up(1,ip10,ip20,ip32)+&
     &                  a10a20a32*up(1,iq)
             field_up(1,ip11,ip20,ip32)=field_up(1,ip11,ip20,ip32)+&
     &                  a11a20a32*up(1,iq)
             field_up(1,ip12,ip20,ip32)=field_up(1,ip12,ip20,ip32)+&
     &                  a12a20a32*up(1,iq)
             field_up(1,ip13,ip20,ip32)=field_up(1,ip13,ip20,ip32)+&
     &                  a13a20a32*up(1,iq)
             field_up(1,ip10,ip21,ip32)=field_up(1,ip10,ip21,ip32)+&
     &                  a10a21a32*up(1,iq)
             field_up(1,ip11,ip21,ip32)=field_up(1,ip11,ip21,ip32)+&
     &                  a11a21a32*up(1,iq)
             field_up(1,ip12,ip21,ip32)=field_up(1,ip12,ip21,ip32)+&
     &                  a12a21a32*up(1,iq)
             field_up(1,ip13,ip21,ip32)=field_up(1,ip13,ip21,ip32)+&
     &                  a13a21a32*up(1,iq)
             field_up(1,ip10,ip22,ip32)=field_up(1,ip10,ip22,ip32)+&
     &                  a10a22a32*up(1,iq)
             field_up(1,ip11,ip22,ip32)=field_up(1,ip11,ip22,ip32)+&
     &                  a11a22a32*up(1,iq)
             field_up(1,ip12,ip22,ip32)=field_up(1,ip12,ip22,ip32)+&
     &                  a12a22a32*up(1,iq)
             field_up(1,ip13,ip22,ip32)=field_up(1,ip13,ip22,ip32)+&
     &                  a13a22a32*up(1,iq)
             field_up(1,ip10,ip23,ip32)=field_up(1,ip10,ip23,ip32)+&
     &                  a10a23a32*up(1,iq)
             field_up(1,ip11,ip23,ip32)=field_up(1,ip11,ip23,ip32)+&
     &                  a11a23a32*up(1,iq)
             field_up(1,ip12,ip23,ip32)=field_up(1,ip12,ip23,ip32)+&
     &                  a12a23a32*up(1,iq)
             field_up(1,ip13,ip23,ip32)=field_up(1,ip13,ip23,ip32)+&
     &                  a13a23a32*up(1,iq)
             field_up(1,ip10,ip20,ip33)=field_up(1,ip10,ip20,ip33)+&
     &                  a10a20a33*up(1,iq)
             field_up(1,ip11,ip20,ip33)=field_up(1,ip11,ip20,ip33)+&
     &                  a11a20a33*up(1,iq)
             field_up(1,ip12,ip20,ip33)=field_up(1,ip12,ip20,ip33)+&
     &                  a12a20a33*up(1,iq)
             field_up(1,ip13,ip20,ip33)=field_up(1,ip13,ip20,ip33)+&
     &                  a13a20a33*up(1,iq)
             field_up(1,ip10,ip21,ip33)=field_up(1,ip10,ip21,ip33)+&
     &                  a10a21a33*up(1,iq)
             field_up(1,ip11,ip21,ip33)=field_up(1,ip11,ip21,ip33)+&
     &                  a11a21a33*up(1,iq)
             field_up(1,ip12,ip21,ip33)=field_up(1,ip12,ip21,ip33)+&
     &                  a12a21a33*up(1,iq)
             field_up(1,ip13,ip21,ip33)=field_up(1,ip13,ip21,ip33)+&
     &                  a13a21a33*up(1,iq)
             field_up(1,ip10,ip22,ip33)=field_up(1,ip10,ip22,ip33)+&
     &                  a10a22a33*up(1,iq)
             field_up(1,ip11,ip22,ip33)=field_up(1,ip11,ip22,ip33)+&
     &                  a11a22a33*up(1,iq)
             field_up(1,ip12,ip22,ip33)=field_up(1,ip12,ip22,ip33)+&
     &                  a12a22a33*up(1,iq)
             field_up(1,ip13,ip22,ip33)=field_up(1,ip13,ip22,ip33)+&
     &                  a13a22a33*up(1,iq)
             field_up(1,ip10,ip23,ip33)=field_up(1,ip10,ip23,ip33)+&
     &                  a10a23a33*up(1,iq)
             field_up(1,ip11,ip23,ip33)=field_up(1,ip11,ip23,ip33)+&
     &                  a11a23a33*up(1,iq)
             field_up(1,ip12,ip23,ip33)=field_up(1,ip12,ip23,ip33)+&
     &                  a12a23a33*up(1,iq)
             field_up(1,ip13,ip23,ip33)=field_up(1,ip13,ip23,ip33)+&
     &                  a13a23a33*up(1,iq)

              ENDDO
              !----------------------------------------------------------------!
              !  All other lda are NOT UNROLLED. Vectorization will be over lda!
              !----------------------------------------------------------------!
           ELSE
               DO ip = 1,store_info(ipatch)
                   iq    = list_sub(ipatch,ip)


                 x01 = xp(1,iq)*dxxi-p%istart(1) + 1
                 x02 = xp(2,iq)*dxyi-p%istart(2) + 1
                 x03 = xp(3,iq)*dxzi-p%istart(3) + 1

                 ip10 = FLOOR(x01)
                 ip20 = FLOOR(x02)
                 ip30 = FLOOR(x03)


                 ip11 = ip10 + 1
                 ip21 = ip20 + 1
                 ip31 = ip30 + 1

                 ip12 = ip11 + 1
                 ip22 = ip21 + 1
                 ip32 = ip31 + 1

                 ip13 = ip11 + 2
                 ip23 = ip21 + 2
                 ip33 = ip31 + 2

                 xp1 = x01-REAL(ip10,mk)
                 xp2 = x02-REAL(ip20,mk)
                 xp3 = x03-REAL(ip30,mk)

                 x10 = xp1 + 1.0_mk
                 x11 = x10 - 1.0_mk
                 x12 = x10 - 2.0_mk
                 x13 = x10 - 3.0_mk

                 x20 = xp2 + 1.0_mk
                 x21 = x20 - 1.0_mk
                 x22 = x20 - 2.0_mk
                 x23 = x20 - 3.0_mk

                 x30 = xp3 + 1.0_mk
                 x31 = x30 - 1.0_mk
                 x32 = x30 - 2.0_mk
                 x33 = x30 - 3.0_mk

                 a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
                 a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
                 a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

                 a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
                 a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
                 a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

                 a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
                 a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
                 a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

                 a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
                 a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
                 a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

                 a10a20a30 = a10*a20*a30
                 a10a20a31 = a10*a20*a31
                 a10a20a32 = a10*a20*a32
                 a10a20a33 = a10*a20*a33
                 a10a21a30 = a10*a21*a30
                 a10a21a31 = a10*a21*a31
                 a10a21a32 = a10*a21*a32
                 a10a21a33 = a10*a21*a33
                 a10a22a30 = a10*a22*a30
                 a10a22a31 = a10*a22*a31
                 a10a22a32 = a10*a22*a32
                 a10a22a33 = a10*a22*a33
                 a10a23a30 = a10*a23*a30
                 a10a23a31 = a10*a23*a31
                 a10a23a32 = a10*a23*a32
                 a10a23a33 = a10*a23*a33
                 a11a20a30 = a11*a20*a30
                 a11a20a31 = a11*a20*a31
                 a11a20a32 = a11*a20*a32
                 a11a20a33 = a11*a20*a33
                 a11a21a30 = a11*a21*a30
                 a11a21a31 = a11*a21*a31
                 a11a21a32 = a11*a21*a32
                 a11a21a33 = a11*a21*a33
                 a11a22a30 = a11*a22*a30
                 a11a22a31 = a11*a22*a31
                 a11a22a32 = a11*a22*a32
                 a11a22a33 = a11*a22*a33
                 a11a23a30 = a11*a23*a30
                 a11a23a31 = a11*a23*a31
                 a11a23a32 = a11*a23*a32
                 a11a23a33 = a11*a23*a33
                 a12a20a30 = a12*a20*a30
                 a12a20a31 = a12*a20*a31
                 a12a20a32 = a12*a20*a32
                 a12a20a33 = a12*a20*a33
                 a12a21a30 = a12*a21*a30
                 a12a21a31 = a12*a21*a31
                 a12a21a32 = a12*a21*a32
                 a12a21a33 = a12*a21*a33
                 a12a22a30 = a12*a22*a30
                 a12a22a31 = a12*a22*a31
                 a12a22a32 = a12*a22*a32
                 a12a22a33 = a12*a22*a33
                 a12a23a30 = a12*a23*a30
                 a12a23a31 = a12*a23*a31
                 a12a23a32 = a12*a23*a32
                 a12a23a33 = a12*a23*a33
                 a13a20a30 = a13*a20*a30
                 a13a20a31 = a13*a20*a31
                 a13a20a32 = a13*a20*a32
                 a13a20a33 = a13*a20*a33
                 a13a21a30 = a13*a21*a30
                 a13a21a31 = a13*a21*a31
                 a13a21a32 = a13*a21*a32
                 a13a21a33 = a13*a21*a33
                 a13a22a30 = a13*a22*a30
                 a13a22a31 = a13*a22*a31
                 a13a22a32 = a13*a22*a32
                 a13a22a33 = a13*a22*a33
                 a13a23a30 = a13*a23*a30
                 a13a23a31 = a13*a23*a31
                 a13a23a32 = a13*a23*a32
                 a13a23a33 = a13*a23*a33

                 DO ldn=1,lda

        field_up(ldn,ip10,ip20,ip30)=field_up(ldn,ip10,ip20,ip30)+&
     &             a10a20a30*up(ldn,iq)
        field_up(ldn,ip11,ip20,ip30)=field_up(ldn,ip11,ip20,ip30)+&
     &             a11a20a30*up(ldn,iq)
        field_up(ldn,ip12,ip20,ip30)=field_up(ldn,ip12,ip20,ip30)+&
     &             a12a20a30*up(ldn,iq)
        field_up(ldn,ip13,ip20,ip30)=field_up(ldn,ip13,ip20,ip30)+&
     &             a13a20a30*up(ldn,iq)
        field_up(ldn,ip10,ip21,ip30)=field_up(ldn,ip10,ip21,ip30)+&
     &             a10a21a30*up(ldn,iq)
        field_up(ldn,ip11,ip21,ip30)=field_up(ldn,ip11,ip21,ip30)+&
     &             a11a21a30*up(ldn,iq)
        field_up(ldn,ip12,ip21,ip30)=field_up(ldn,ip12,ip21,ip30)+&
     &             a12a21a30*up(ldn,iq)
        field_up(ldn,ip13,ip21,ip30)=field_up(ldn,ip13,ip21,ip30)+&
     &             a13a21a30*up(ldn,iq)
        field_up(ldn,ip10,ip22,ip30)=field_up(ldn,ip10,ip22,ip30)+&
     &             a10a22a30*up(ldn,iq)
        field_up(ldn,ip11,ip22,ip30)=field_up(ldn,ip11,ip22,ip30)+&
     &             a11a22a30*up(ldn,iq)
        field_up(ldn,ip12,ip22,ip30)=field_up(ldn,ip12,ip22,ip30)+&
     &             a12a22a30*up(ldn,iq)
        field_up(ldn,ip13,ip22,ip30)=field_up(ldn,ip13,ip22,ip30)+&
     &             a13a22a30*up(ldn,iq)
        field_up(ldn,ip10,ip23,ip30)=field_up(ldn,ip10,ip23,ip30)+&
     &             a10a23a30*up(ldn,iq)
        field_up(ldn,ip11,ip23,ip30)=field_up(ldn,ip11,ip23,ip30)+&
     &             a11a23a30*up(ldn,iq)
        field_up(ldn,ip12,ip23,ip30)=field_up(ldn,ip12,ip23,ip30)+&
     &             a12a23a30*up(ldn,iq)
        field_up(ldn,ip13,ip23,ip30)=field_up(ldn,ip13,ip23,ip30)+&
     &             a13a23a30*up(ldn,iq)
        field_up(ldn,ip10,ip20,ip31)=field_up(ldn,ip10,ip20,ip31)+&
     &             a10a20a31*up(ldn,iq)
        field_up(ldn,ip11,ip20,ip31)=field_up(ldn,ip11,ip20,ip31)+&
     &             a11a20a31*up(ldn,iq)
        field_up(ldn,ip12,ip20,ip31)=field_up(ldn,ip12,ip20,ip31)+&
     &             a12a20a31*up(ldn,iq)
        field_up(ldn,ip13,ip20,ip31)=field_up(ldn,ip13,ip20,ip31)+&
     &             a13a20a31*up(ldn,iq)
        field_up(ldn,ip10,ip21,ip31)=field_up(ldn,ip10,ip21,ip31)+&
     &             a10a21a31*up(ldn,iq)
        field_up(ldn,ip11,ip21,ip31)=field_up(ldn,ip11,ip21,ip31)+&
     &             a11a21a31*up(ldn,iq)
        field_up(ldn,ip12,ip21,ip31)=field_up(ldn,ip12,ip21,ip31)+&
     &             a12a21a31*up(ldn,iq)
        field_up(ldn,ip13,ip21,ip31)=field_up(ldn,ip13,ip21,ip31)+&
     &             a13a21a31*up(ldn,iq)
        field_up(ldn,ip10,ip22,ip31)=field_up(ldn,ip10,ip22,ip31)+&
     &             a10a22a31*up(ldn,iq)
        field_up(ldn,ip11,ip22,ip31)=field_up(ldn,ip11,ip22,ip31)+&
     &             a11a22a31*up(ldn,iq)
        field_up(ldn,ip12,ip22,ip31)=field_up(ldn,ip12,ip22,ip31)+&
     &             a12a22a31*up(ldn,iq)
        field_up(ldn,ip13,ip22,ip31)=field_up(ldn,ip13,ip22,ip31)+&
     &             a13a22a31*up(ldn,iq)
        field_up(ldn,ip10,ip23,ip31)=field_up(ldn,ip10,ip23,ip31)+&
     &             a10a23a31*up(ldn,iq)
        field_up(ldn,ip11,ip23,ip31)=field_up(ldn,ip11,ip23,ip31)+&
     &             a11a23a31*up(ldn,iq)
        field_up(ldn,ip12,ip23,ip31)=field_up(ldn,ip12,ip23,ip31)+&
     &             a12a23a31*up(ldn,iq)
        field_up(ldn,ip13,ip23,ip31)=field_up(ldn,ip13,ip23,ip31)+&
     &             a13a23a31*up(ldn,iq)
        field_up(ldn,ip10,ip20,ip32)=field_up(ldn,ip10,ip20,ip32)+&
     &             a10a20a32*up(ldn,iq)
        field_up(ldn,ip11,ip20,ip32)=field_up(ldn,ip11,ip20,ip32)+&
     &             a11a20a32*up(ldn,iq)
        field_up(ldn,ip12,ip20,ip32)=field_up(ldn,ip12,ip20,ip32)+&
     &             a12a20a32*up(ldn,iq)
        field_up(ldn,ip13,ip20,ip32)=field_up(ldn,ip13,ip20,ip32)+&
     &             a13a20a32*up(ldn,iq)
        field_up(ldn,ip10,ip21,ip32)=field_up(ldn,ip10,ip21,ip32)+&
     &             a10a21a32*up(ldn,iq)
        field_up(ldn,ip11,ip21,ip32)=field_up(ldn,ip11,ip21,ip32)+&
     &             a11a21a32*up(ldn,iq)
        field_up(ldn,ip12,ip21,ip32)=field_up(ldn,ip12,ip21,ip32)+&
     &             a12a21a32*up(ldn,iq)
        field_up(ldn,ip13,ip21,ip32)=field_up(ldn,ip13,ip21,ip32)+&
     &             a13a21a32*up(ldn,iq)
        field_up(ldn,ip10,ip22,ip32)=field_up(ldn,ip10,ip22,ip32)+&
     &             a10a22a32*up(ldn,iq)
        field_up(ldn,ip11,ip22,ip32)=field_up(ldn,ip11,ip22,ip32)+&
     &             a11a22a32*up(ldn,iq)
        field_up(ldn,ip12,ip22,ip32)=field_up(ldn,ip12,ip22,ip32)+&
     &             a12a22a32*up(ldn,iq)
        field_up(ldn,ip13,ip22,ip32)=field_up(ldn,ip13,ip22,ip32)+&
     &             a13a22a32*up(ldn,iq)
        field_up(ldn,ip10,ip23,ip32)=field_up(ldn,ip10,ip23,ip32)+&
     &             a10a23a32*up(ldn,iq)
        field_up(ldn,ip11,ip23,ip32)=field_up(ldn,ip11,ip23,ip32)+&
     &             a11a23a32*up(ldn,iq)
        field_up(ldn,ip12,ip23,ip32)=field_up(ldn,ip12,ip23,ip32)+&
     &             a12a23a32*up(ldn,iq)
        field_up(ldn,ip13,ip23,ip32)=field_up(ldn,ip13,ip23,ip32)+&
     &             a13a23a32*up(ldn,iq)
        field_up(ldn,ip10,ip20,ip33)=field_up(ldn,ip10,ip20,ip33)+&
     &             a10a20a33*up(ldn,iq)
        field_up(ldn,ip11,ip20,ip33)=field_up(ldn,ip11,ip20,ip33)+&
     &             a11a20a33*up(ldn,iq)
        field_up(ldn,ip12,ip20,ip33)=field_up(ldn,ip12,ip20,ip33)+&
     &             a12a20a33*up(ldn,iq)
        field_up(ldn,ip13,ip20,ip33)=field_up(ldn,ip13,ip20,ip33)+&
     &             a13a20a33*up(ldn,iq)
        field_up(ldn,ip10,ip21,ip33)=field_up(ldn,ip10,ip21,ip33)+&
     &             a10a21a33*up(ldn,iq)
        field_up(ldn,ip11,ip21,ip33)=field_up(ldn,ip11,ip21,ip33)+&
     &             a11a21a33*up(ldn,iq)
        field_up(ldn,ip12,ip21,ip33)=field_up(ldn,ip12,ip21,ip33)+&
     &             a12a21a33*up(ldn,iq)
        field_up(ldn,ip13,ip21,ip33)=field_up(ldn,ip13,ip21,ip33)+&
     &             a13a21a33*up(ldn,iq)
        field_up(ldn,ip10,ip22,ip33)=field_up(ldn,ip10,ip22,ip33)+&
     &             a10a22a33*up(ldn,iq)
        field_up(ldn,ip11,ip22,ip33)=field_up(ldn,ip11,ip22,ip33)+&
     &             a11a22a33*up(ldn,iq)
        field_up(ldn,ip12,ip22,ip33)=field_up(ldn,ip12,ip22,ip33)+&
     &             a12a22a33*up(ldn,iq)
        field_up(ldn,ip13,ip22,ip33)=field_up(ldn,ip13,ip22,ip33)+&
     &             a13a22a33*up(ldn,iq)
        field_up(ldn,ip10,ip23,ip33)=field_up(ldn,ip10,ip23,ip33)+&
     &             a10a23a33*up(ldn,iq)
        field_up(ldn,ip11,ip23,ip33)=field_up(ldn,ip11,ip23,ip33)+&
     &             a11a23a33*up(ldn,iq)
        field_up(ldn,ip12,ip23,ip33)=field_up(ldn,ip12,ip23,ip33)+&
     &             a12a23a33*up(ldn,iq)
        field_up(ldn,ip13,ip23,ip33)=field_up(ldn,ip13,ip23,ip33)+&
     &             a13a23a33*up(ldn,iq)

                 ENDDO    ! lda
              ENDDO        ! iq
           END IF          ! unrolled lda cases
#elif __MODE == __SCA
         DO ip = 1,store_info(ipatch)
              iq    = list_sub(ipatch,ip)

              x01 = xp(1,iq)*dxxi-p%istart(1) + 1

              ip10 = FLOOR(x01)
              ip11 = ip10 + 1
              ip12 = ip11 + 1
              ip13 = ip11 + 2

              x01_l = x01 - p%lo_a(1)
              x01_h = p%hi_a(1) - x01

              ! some kind of Heaviside function
              !    for particles near the left boundary
              l10 = 1+FLOOR(x01_l      )-INT(x01_l)
              l11 = 1+FLOOR(x01_l+1._mk)-INT(x01_l+1._mk)
              l12 = 1+FLOOR(x01_l+2._mk)-INT(x01_l+2._mk)

              !    for particles near the right boundary
              h11 = INT(x01_h- 1._mk)-FLOOR(x01_h - 1._mk)
              h12 = INT(x01_h- 2._mk)-FLOOR(x01_h - 2._mk)
              h13 = INT(x01_h- 3._mk)-FLOOR(x01_h - 3._mk)

              x10 = (x01-REAL(ip10,mk)) * l10       + 1._mk
              x11 = (x01-REAL(ip11,mk)) * l11 * h11 + 1._mk
              x12 = (x01-REAL(ip12,mk)) * l12 * h12 + 1._mk
              x13 = (x01-REAL(ip13,mk))       * h13 + 1._mk

              ip10 = (ip10-p%lo_a(1) + ABS(ip10-p%lo_a(1)))/2 + p%lo_a(1)
              ip11 = (ip11-p%lo_a(1) + ABS(ip11-p%lo_a(1)))/2 + p%lo_a(1)
              ip12 = (ip12-p%lo_a(1) + ABS(ip12-p%lo_a(1)))/2 + p%lo_a(1)

              ip11 = p%hi_a(1) - (p%hi_a(1)-ip11 + ABS(p%hi_a(1)-ip11))/2
              ip12 = p%hi_a(1) - (p%hi_a(1)-ip12 + ABS(p%hi_a(1)-ip12))/2
              ip13 = p%hi_a(1) - (p%hi_a(1)-ip13 + ABS(p%hi_a(1)-ip13))/2

              !stdout("p%lo_a",'p%lo_a')
              !stdout("p%hi_a",'p%hi_a')
              !stdout("ip10",ip10,'floor(x01)')
              !stdout("ip11",ip11,'floor(x01)+1')
              !stdout("ip12",ip12,'floor(x01)+2')
              !stdout("ip13",ip13,'floor(x01)+3')
              !stdout("x10",x10)
              !stdout("x11",x11)
              !stdout("x12",x12)
              !stdout("x13",x13)
              !stdout("l1j",l10,l11,l12)
              !stdout("h1j",h11,h12,h13)

              !check_true("ip10.GE.p%lo_a(1)")
              !check_true("ip10.LE.p%hi_a(1)")
              !check_true("ip11.GE.p%lo_a(1)")
              !check_true("ip11.LE.p%hi_a(1)")
              !check_true("ip12.GE.p%lo_a(1)")
              !check_true("ip12.LE.p%hi_a(1)")
              !check_true("ip13.GE.p%lo_a(1)")
              !check_true("ip13.LE.p%hi_a(1)")


              x02 = xp(2,iq)*dxyi-p%istart(2) + 1

              ip20 = FLOOR(x02)
              ip21 = ip20 + 1
              ip22 = ip21 + 1
              ip23 = ip21 + 2

              x02_l = x02 - p%lo_a(2)
              x02_h = p%hi_a(2) - x02

              ! some kind of Heaviside function
              !    for particles near the left boundary
              l20 = 1+FLOOR(x02_l      )-INT(x02_l)
              l21 = 1+FLOOR(x02_l+1._mk)-INT(x02_l+1._mk)
              l22 = 1+FLOOR(x02_l+2._mk)-INT(x02_l+2._mk)

              !    for particles near the right boundary
              h21 = INT(x02_h - 1._mk)-FLOOR(x02_h - 1._mk)
              h22 = INT(x02_h - 2._mk)-FLOOR(x02_h - 2._mk)
              h23 = INT(x02_h - 3._mk)-FLOOR(x02_h - 3._mk)

              x20 = (x02-REAL(ip20,mk)) * l20       + 1._mk
              x21 = (x02-REAL(ip21,mk)) * l21 * h21 + 1._mk
              x22 = (x02-REAL(ip22,mk)) * l22 * h22 + 1._mk
              x23 = (x02-REAL(ip23,mk))       * h23 + 1._mk

              ip20 = (ip20-p%lo_a(2) + ABS(ip20-p%lo_a(2)))/2 + p%lo_a(2)
              ip21 = (ip21-p%lo_a(2) + ABS(ip21-p%lo_a(2)))/2 + p%lo_a(2)
              ip22 = (ip22-p%lo_a(2) + ABS(ip22-p%lo_a(2)))/2 + p%lo_a(2)

              ip21 = p%hi_a(2) - (p%hi_a(2)-ip21 + ABS(p%hi_a(2)-ip21))/2
              ip22 = p%hi_a(2) - (p%hi_a(2)-ip22 + ABS(p%hi_a(2)-ip22))/2
              ip23 = p%hi_a(2) - (p%hi_a(2)-ip23 + ABS(p%hi_a(2)-ip23))/2



              x03 = xp(3,iq)*dxzi-p%istart(3) + 1

              ip30 = FLOOR(x03)
              ip31 = ip30 + 1
              ip32 = ip31 + 1
              ip33 = ip31 + 2

              x03_l = x03 - p%lo_a(3)
              x03_h = p%hi_a(3) - x03

              ! some kind of Heaviside function
              !    for particles near the left boundary
              l30 = 1+FLOOR(x03_l      )-INT(x03_l)
              l31 = 1+FLOOR(x03_l+1._mk)-INT(x03_l+1._mk)
              l32 = 1+FLOOR(x03_l+2._mk)-INT(x03_l+2._mk)

              !    for particles near the right boundary
              h31 = INT(x03_h - 1._mk)-FLOOR(x03_h - 1._mk)
              h32 = INT(x03_h - 2._mk)-FLOOR(x03_h - 2._mk)
              h33 = INT(x03_h - 3._mk)-FLOOR(x03_h - 3._mk)

              x30 = (x01-REAL(ip30,mk)) * l30       + 1._mk
              x31 = (x01-REAL(ip31,mk)) * l31 * h31 + 1._mk
              x32 = (x01-REAL(ip32,mk)) * l32 * h32 + 1._mk
              x33 = (x01-REAL(ip33,mk))       * h33 + 1._mk

              ip30 = (ip30-p%lo_a(3) + ABS(ip30-p%lo_a(3)))/2 + p%lo_a(3)
              ip31 = (ip31-p%lo_a(3) + ABS(ip31-p%lo_a(3)))/2 + p%lo_a(3)
              ip32 = (ip32-p%lo_a(3) + ABS(ip32-p%lo_a(3)))/2 + p%lo_a(3)

              ip31 = p%hi_a(3) - (p%hi_a(3)-ip31 + ABS(p%hi_a(3)-ip31))/2
              ip32 = p%hi_a(3) - (p%hi_a(3)-ip32 + ABS(p%hi_a(3)-ip32))/2
              ip33 = p%hi_a(3) - (p%hi_a(3)-ip33 + ABS(p%hi_a(3)-ip33))/2




              a10 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x10)*x10)*x10
              a20 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x20)*x20)*x20
              a30 = 2.0_mk + (-4.0_mk+(2.5_mk-0.5_mk*x30)*x30)*x30

              a11 = 1.0_mk + (-2.5_mk+1.5_mk*x11)*x11**2
              a21 = 1.0_mk + (-2.5_mk+1.5_mk*x21)*x21**2
              a31 = 1.0_mk + (-2.5_mk+1.5_mk*x31)*x31**2

              a12 = 1.0_mk + (-2.5_mk-1.5_mk*x12)*x12**2
              a22 = 1.0_mk + (-2.5_mk-1.5_mk*x22)*x22**2
              a32 = 1.0_mk + (-2.5_mk-1.5_mk*x32)*x32**2

              a13 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x13)*x13)*x13
              a23 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x23)*x23)*x23
              a33 = 2.0_mk + (4.0_mk + (2.5_mk+0.5_mk*x33)*x33)*x33

              a10a20a30 = a10*a20*a30
              a10a20a31 = a10*a20*a31
              a10a20a32 = a10*a20*a32
              a10a20a33 = a10*a20*a33
              a10a21a30 = a10*a21*a30
              a10a21a31 = a10*a21*a31
              a10a21a32 = a10*a21*a32
              a10a21a33 = a10*a21*a33
              a10a22a30 = a10*a22*a30
              a10a22a31 = a10*a22*a31
              a10a22a32 = a10*a22*a32
              a10a22a33 = a10*a22*a33
              a10a23a30 = a10*a23*a30
              a10a23a31 = a10*a23*a31
              a10a23a32 = a10*a23*a32
              a10a23a33 = a10*a23*a33
              a11a20a30 = a11*a20*a30
              a11a20a31 = a11*a20*a31
              a11a20a32 = a11*a20*a32
              a11a20a33 = a11*a20*a33
              a11a21a30 = a11*a21*a30
              a11a21a31 = a11*a21*a31
              a11a21a32 = a11*a21*a32
              a11a21a33 = a11*a21*a33
              a11a22a30 = a11*a22*a30
              a11a22a31 = a11*a22*a31
              a11a22a32 = a11*a22*a32
              a11a22a33 = a11*a22*a33
              a11a23a30 = a11*a23*a30
              a11a23a31 = a11*a23*a31
              a11a23a32 = a11*a23*a32
              a11a23a33 = a11*a23*a33
              a12a20a30 = a12*a20*a30
              a12a20a31 = a12*a20*a31
              a12a20a32 = a12*a20*a32
              a12a20a33 = a12*a20*a33
              a12a21a30 = a12*a21*a30
              a12a21a31 = a12*a21*a31
              a12a21a32 = a12*a21*a32
              a12a21a33 = a12*a21*a33
              a12a22a30 = a12*a22*a30
              a12a22a31 = a12*a22*a31
              a12a22a32 = a12*a22*a32
              a12a22a33 = a12*a22*a33
              a12a23a30 = a12*a23*a30
              a12a23a31 = a12*a23*a31
              a12a23a32 = a12*a23*a32
              a12a23a33 = a12*a23*a33
              a13a20a30 = a13*a20*a30
              a13a20a31 = a13*a20*a31
              a13a20a32 = a13*a20*a32
              a13a20a33 = a13*a20*a33
              a13a21a30 = a13*a21*a30
              a13a21a31 = a13*a21*a31
              a13a21a32 = a13*a21*a32
              a13a21a33 = a13*a21*a33
              a13a22a30 = a13*a22*a30
              a13a22a31 = a13*a22*a31
              a13a22a32 = a13*a22*a32
              a13a22a33 = a13*a22*a33
              a13a23a30 = a13*a23*a30
              a13a23a31 = a13*a23*a31
              a13a23a32 = a13*a23*a32
              a13a23a33 = a13*a23*a33

              field_up(ip10,ip20,ip30)=field_up(ip10,ip20,ip30)+&
     &                   a10a20a30*up(iq)
              field_up(ip11,ip20,ip30)=field_up(ip11,ip20,ip30)+&
     &                   a11a20a30*up(iq)
              field_up(ip12,ip20,ip30)=field_up(ip12,ip20,ip30)+&
     &                   a12a20a30*up(iq)
              field_up(ip13,ip20,ip30)=field_up(ip13,ip20,ip30)+&
     &                   a13a20a30*up(iq)
              field_up(ip10,ip21,ip30)=field_up(ip10,ip21,ip30)+&
     &                   a10a21a30*up(iq)
              field_up(ip11,ip21,ip30)=field_up(ip11,ip21,ip30)+&
     &                   a11a21a30*up(iq)
              field_up(ip12,ip21,ip30)=field_up(ip12,ip21,ip30)+&
     &                   a12a21a30*up(iq)
              field_up(ip13,ip21,ip30)=field_up(ip13,ip21,ip30)+&
     &                   a13a21a30*up(iq)
              field_up(ip10,ip22,ip30)=field_up(ip10,ip22,ip30)+&
     &                   a10a22a30*up(iq)
              field_up(ip11,ip22,ip30)=field_up(ip11,ip22,ip30)+&
     &                   a11a22a30*up(iq)
              field_up(ip12,ip22,ip30)=field_up(ip12,ip22,ip30)+&
     &                   a12a22a30*up(iq)
              field_up(ip13,ip22,ip30)=field_up(ip13,ip22,ip30)+&
     &                   a13a22a30*up(iq)
              field_up(ip10,ip23,ip30)=field_up(ip10,ip23,ip30)+&
     &                   a10a23a30*up(iq)
              field_up(ip11,ip23,ip30)=field_up(ip11,ip23,ip30)+&
     &                   a11a23a30*up(iq)
              field_up(ip12,ip23,ip30)=field_up(ip12,ip23,ip30)+&
     &                   a12a23a30*up(iq)
              field_up(ip13,ip23,ip30)=field_up(ip13,ip23,ip30)+&
     &                   a13a23a30*up(iq)
              field_up(ip10,ip20,ip31)=field_up(ip10,ip20,ip31)+&
     &                   a10a20a31*up(iq)
              field_up(ip11,ip20,ip31)=field_up(ip11,ip20,ip31)+&
     &                   a11a20a31*up(iq)
              field_up(ip12,ip20,ip31)=field_up(ip12,ip20,ip31)+&
     &                   a12a20a31*up(iq)
              field_up(ip13,ip20,ip31)=field_up(ip13,ip20,ip31)+&
     &                   a13a20a31*up(iq)
              field_up(ip10,ip21,ip31)=field_up(ip10,ip21,ip31)+&
     &                   a10a21a31*up(iq)
              field_up(ip11,ip21,ip31)=field_up(ip11,ip21,ip31)+&
     &                   a11a21a31*up(iq)
              field_up(ip12,ip21,ip31)=field_up(ip12,ip21,ip31)+&
     &                   a12a21a31*up(iq)
              field_up(ip13,ip21,ip31)=field_up(ip13,ip21,ip31)+&
     &                   a13a21a31*up(iq)
              field_up(ip10,ip22,ip31)=field_up(ip10,ip22,ip31)+&
     &                   a10a22a31*up(iq)
              field_up(ip11,ip22,ip31)=field_up(ip11,ip22,ip31)+&
     &                   a11a22a31*up(iq)
              field_up(ip12,ip22,ip31)=field_up(ip12,ip22,ip31)+&
     &                   a12a22a31*up(iq)
              field_up(ip13,ip22,ip31)=field_up(ip13,ip22,ip31)+&
     &                   a13a22a31*up(iq)
              field_up(ip10,ip23,ip31)=field_up(ip10,ip23,ip31)+&
     &                   a10a23a31*up(iq)
              field_up(ip11,ip23,ip31)=field_up(ip11,ip23,ip31)+&
     &                   a11a23a31*up(iq)
              field_up(ip12,ip23,ip31)=field_up(ip12,ip23,ip31)+&
     &                   a12a23a31*up(iq)
              field_up(ip13,ip23,ip31)=field_up(ip13,ip23,ip31)+&
     &                   a13a23a31*up(iq)
              field_up(ip10,ip20,ip32)=field_up(ip10,ip20,ip32)+&
     &                   a10a20a32*up(iq)
              field_up(ip11,ip20,ip32)=field_up(ip11,ip20,ip32)+&
     &                   a11a20a32*up(iq)
              field_up(ip12,ip20,ip32)=field_up(ip12,ip20,ip32)+&
     &                   a12a20a32*up(iq)
              field_up(ip13,ip20,ip32)=field_up(ip13,ip20,ip32)+&
     &                   a13a20a32*up(iq)
              field_up(ip10,ip21,ip32)=field_up(ip10,ip21,ip32)+&
     &                   a10a21a32*up(iq)
              field_up(ip11,ip21,ip32)=field_up(ip11,ip21,ip32)+&
     &                   a11a21a32*up(iq)
              field_up(ip12,ip21,ip32)=field_up(ip12,ip21,ip32)+&
     &                   a12a21a32*up(iq)
              field_up(ip13,ip21,ip32)=field_up(ip13,ip21,ip32)+&
     &                   a13a21a32*up(iq)
              field_up(ip10,ip22,ip32)=field_up(ip10,ip22,ip32)+&
     &                   a10a22a32*up(iq)
              field_up(ip11,ip22,ip32)=field_up(ip11,ip22,ip32)+&
     &                   a11a22a32*up(iq)
              field_up(ip12,ip22,ip32)=field_up(ip12,ip22,ip32)+&
     &                   a12a22a32*up(iq)
              field_up(ip13,ip22,ip32)=field_up(ip13,ip22,ip32)+&
     &                   a13a22a32*up(iq)
              field_up(ip10,ip23,ip32)=field_up(ip10,ip23,ip32)+&
     &                   a10a23a32*up(iq)
              field_up(ip11,ip23,ip32)=field_up(ip11,ip23,ip32)+&
     &                   a11a23a32*up(iq)
              field_up(ip12,ip23,ip32)=field_up(ip12,ip23,ip32)+&
     &                   a12a23a32*up(iq)
              field_up(ip13,ip23,ip32)=field_up(ip13,ip23,ip32)+&
     &                   a13a23a32*up(iq)
              field_up(ip10,ip20,ip33)=field_up(ip10,ip20,ip33)+&
     &                   a10a20a33*up(iq)
              field_up(ip11,ip20,ip33)=field_up(ip11,ip20,ip33)+&
     &                   a11a20a33*up(iq)
              field_up(ip12,ip20,ip33)=field_up(ip12,ip20,ip33)+&
     &                   a12a20a33*up(iq)
              field_up(ip13,ip20,ip33)=field_up(ip13,ip20,ip33)+&
     &                   a13a20a33*up(iq)
              field_up(ip10,ip21,ip33)=field_up(ip10,ip21,ip33)+&
     &                   a10a21a33*up(iq)
              field_up(ip11,ip21,ip33)=field_up(ip11,ip21,ip33)+&
     &                   a11a21a33*up(iq)
              field_up(ip12,ip21,ip33)=field_up(ip12,ip21,ip33)+&
     &                   a12a21a33*up(iq)
              field_up(ip13,ip21,ip33)=field_up(ip13,ip21,ip33)+&
     &                   a13a21a33*up(iq)
              field_up(ip10,ip22,ip33)=field_up(ip10,ip22,ip33)+&
     &                   a10a22a33*up(iq)
              field_up(ip11,ip22,ip33)=field_up(ip11,ip22,ip33)+&
     &                   a11a22a33*up(iq)
              field_up(ip12,ip22,ip33)=field_up(ip12,ip22,ip33)+&
     &                   a12a22a33*up(iq)
              field_up(ip13,ip22,ip33)=field_up(ip13,ip22,ip33)+&
     &                   a13a22a33*up(iq)
              field_up(ip10,ip23,ip33)=field_up(ip10,ip23,ip33)+&
     &                   a10a23a33*up(iq)
              field_up(ip11,ip23,ip33)=field_up(ip11,ip23,ip33)+&
     &                   a11a23a33*up(iq)
              field_up(ip12,ip23,ip33)=field_up(ip12,ip23,ip33)+&
     &                   a12a23a33*up(iq)
              field_up(ip13,ip23,ip33)=field_up(ip13,ip23,ip33)+&
     &                   a13a23a33*up(iq)
           ENDDO        ! iq
#endif

#elif __DIME == __2D
         DO ip = 1,store_info(ipatch)
              iq    = list_sub(ipatch,ip)

              x01 = xp(1,iq)*dxxi-p%istart(1) + 1
              x02 = xp(2,iq)*dxyi-p%istart(2) + 1

              ip1 = FLOOR(x01)+1
              ip2 = FLOOR(x02)+1
!
              xp1 = x01-FLOOR(x01)
              xp2 = x02-FLOOR(x02)

              DO jj = -1,2
                  IF (jj+ip2.LT.p%lo_a(2)) CYCLE
                  IF (jj+ip2.GT.p%hi_a(2)) CYCLE

                 x2 = ABS(xp2 - REAL(jj,mk))

                 IF(x2.LT.1.0_mk) THEN
                    wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                 ELSE
                    wx2 = 2.0_mk + (-4.0_mk + &
     &                         (2.5_mk - 0.5_mk * x2)*x2)*x2
                 END IF

                 DO ii    = - 1,2
                  IF (ii+ip1.LT.p%lo_a(1)) CYCLE
                  IF (ii+ip1.GT.p%hi_a(1)) CYCLE

                    x1 = ABS(xp1 - REAL(ii,mk))

                    IF(x1.LT.1.0_MK) THEN
                       wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                             1.5_mk*x1)
                    ELSE
                       wx1 =  2.0_mk + (-4.0_mk + &
     &                             (2.5_mk - 0.5_mk*x1)*x1)*x1
                    END IF

#if __MODE == __SCA
                    field_up(ii+ip1,jj+ip2) &
     &                         = field_up(ii+ip1,jj+ip2) &
     &                                          + wx1*wx2*up(iq)
#else
                    DO ldn=1,lda
                       field_up(ldn,ii+ip1,jj+ip2) &
     &                            = field_up(ldn,ii+ip1,jj+ip2) &
     &                                             + wx1*wx2*up(ldn,iq)
                    ENDDO
#endif
                 END DO !ii
              END DO !jj
           END DO !ip
#endif

          p => Mesh%subpatch%next()
          ipatch = ipatch + 1
        ENDDO subpatch       ! loop over subpatches

      end_subroutine()
      RETURN

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE p2m_interp_mp4_dv_3d
#endif
#endif
#endif
