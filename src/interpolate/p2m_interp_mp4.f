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
      INTEGER,PARAMETER                      :: one=1
      INTEGER                                :: ip10l,ip20l,ip03l
      INTEGER                                :: ip10h,ip20h,ip03h
      INTEGER                                :: l10,l11,l12,l13
      INTEGER                                :: h10,h11,h12,h13
      INTEGER                                :: l20,l21,l22,l23
      INTEGER                                :: h20,h21,h22,h23
      INTEGER                                :: l30,l31,l32,l33
      INTEGER                                :: h30,h31,h32,h33

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

#include "comp_weights_mp4.inc"

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

#include "comp_weights_mp4.inc"

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

#include "comp_weights_mp4.inc"

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

#include "comp_weights_mp4.inc"

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

#include "comp_weights_mp4.inc"

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
