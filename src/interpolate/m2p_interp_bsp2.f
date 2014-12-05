      !-------------------------------------------------------------------------
      !     Subroutine   :                   m2p_interp_bsp2
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
      SUBROUTINE m2p_interp_bsp2_ss_2d(Mesh,Field,field_up,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_ds_2d(Mesh,Field,field_up,xp,up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_sv_2d(Mesh,Field,field_up,lda,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_dv_2d(Mesh,Field,field_up,lda,xp,up,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_ss_3d(Mesh,Field,field_up,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_ds_3d(Mesh,Field,field_up,xp,up,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_sv_3d(Mesh,Field,field_up,lda,xp,up,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_bsp2_dv_3d(Mesh,Field,field_up,lda,xp,up,info)
#endif
#endif
#endif
     !!! Mesh to particle interpolation following the 2nd order B-spline scheme.
     !!!
     !!! The interpolation scheme is only implemented for 2D and 3D spaces. To
     !!! increase performance the inner loops over the number of properties to
     !!! be interpolated are unrolled for 2,3,4 and 5-vectors.
     !!!
     !!! [NOTE]
     !!! This routine only performs the actual interpolation. It should not be
     !!! called directly by the user but instead the `ppm_interp_m2p`
     !!! routine should be used with the kernel argument set to
     !!! `ppm_param_rmsh_kernel_bsp2`.
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      CLASS(ppm_t_equi_mesh_)                          :: Mesh
      !!! Mesh
      CLASS(ppm_t_field_)                              :: Field
      !!! Field
#if   __MODE == __SCA
      REAL(MK) , DIMENSION(:)         , POINTER        :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:  ) ,     POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:) ,     POINTER        :: field_up
#endif
      !!! field from which to interpolate
#elif __MODE == __VEC
      INTEGER                         , INTENT(IN   )  :: lda
      !!! leading dimension of up
      REAL(MK) , DIMENSION(:,:)       , POINTER        :: up
      !!! particle weights onto which to interpolate
#if   __DIME == __2D
      REAL(MK) , DIMENSION(:,:,:  ) ,   POINTER        :: field_up
#elif __DIME == __3D
      REAL(MK) , DIMENSION(:,:,:,:) ,   POINTER        :: field_up
#endif
      !!! field from which to interpolate
#endif
      REAL(MK), DIMENSION(:,:)       ,   INTENT(IN   ) :: xp
      !!! particle positions
      INTEGER                        ,   INTENT(  OUT) :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      REAL(MK),  DIMENSION(ppm_dim) :: dxi
      REAL(MK)                      :: x1,x2,x3

      INTEGER :: i,j,k,ii,jj,kk
      INTEGER :: ip,ip1,ip2,ip3
      INTEGER :: iq,ipatch

      ! aliases
      CLASS(ppm_t_subpatch_), POINTER :: p

      REAL(MK), DIMENSION(ppm_dim) :: x0
      REAL(MK)                     :: xp1,xp2,xp3

      INTEGER :: ldn

      !-------------------------------------------------------------------------
      !  Variables for unrolled versions
      !-------------------------------------------------------------------------

#if   __DIME == __2D
      REAL(MK) :: x10,x11,x20,x21
      REAL(MK) :: a10,a11,a20,a21
      INTEGER  :: ip10,ip11,ip20,ip21
      REAL(MK) :: a10a20
      REAL(MK) :: a10a21
      REAL(MK) :: a11a20
      REAL(MK) :: a11a21
#elif __DIME == __3D
      REAL(MK) :: x10,x11,x12,x13,x20,x21,x22,x23,x30,x31,x32,x33
      REAL(MK) :: a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33
      INTEGER  :: ip10,ip11,ip12,ip13,ip20,ip21,ip22,ip23,ip30,ip31,ip32,ip33
      REAL(MK) :: a10a20a30
      REAL(MK) :: a10a20a31
      REAL(MK) :: a10a20a32
      REAL(MK) :: a10a20a33
      REAL(MK) :: a10a21a30
      REAL(MK) :: a10a21a31
      REAL(MK) :: a10a21a32
      REAL(MK) :: a10a21a33
      REAL(MK) :: a10a22a30
      REAL(MK) :: a10a22a31
      REAL(MK) :: a10a22a32
      REAL(MK) :: a10a22a33
      REAL(MK) :: a10a23a30
      REAL(MK) :: a10a23a31
      REAL(MK) :: a10a23a32
      REAL(MK) :: a10a23a33
      REAL(MK) :: a11a20a30
      REAL(MK) :: a11a20a31
      REAL(MK) :: a11a20a32
      REAL(MK) :: a11a20a33
      REAL(MK) :: a11a21a30
      REAL(MK) :: a11a21a31
      REAL(MK) :: a11a21a32
      REAL(MK) :: a11a21a33
      REAL(MK) :: a11a22a30
      REAL(MK) :: a11a22a31
      REAL(MK) :: a11a22a32
      REAL(MK) :: a11a22a33
      REAL(MK) :: a11a23a30
      REAL(MK) :: a11a23a31
      REAL(MK) :: a11a23a32
      REAL(MK) :: a11a23a33
      REAL(MK) :: a12a20a30
      REAL(MK) :: a12a20a31
      REAL(MK) :: a12a20a32
      REAL(MK) :: a12a20a33
      REAL(MK) :: a12a21a30
      REAL(MK) :: a12a21a31
      REAL(MK) :: a12a21a32
      REAL(MK) :: a12a21a33
      REAL(MK) :: a12a22a30
      REAL(MK) :: a12a22a31
      REAL(MK) :: a12a22a32
      REAL(MK) :: a12a22a33
      REAL(MK) :: a12a23a30
      REAL(MK) :: a12a23a31
      REAL(MK) :: a12a23a32
      REAL(MK) :: a12a23a33
      REAL(MK) :: a13a20a30
      REAL(MK) :: a13a20a31
      REAL(MK) :: a13a20a32
      REAL(MK) :: a13a20a33
      REAL(MK) :: a13a21a30
      REAL(MK) :: a13a21a31
      REAL(MK) :: a13a21a32
      REAL(MK) :: a13a21a33
      REAL(MK) :: a13a22a30
      REAL(MK) :: a13a22a31
      REAL(MK) :: a13a22a32
      REAL(MK) :: a13a22a33
      REAL(MK) :: a13a23a30
      REAL(MK) :: a13a23a31
      REAL(MK) :: a13a23a32
      REAL(MK) :: a13a23a33
#endif


      start_subroutine("m2p_interp_mp4")

      dxi = 1.0_MK/Mesh%h

      !  loop over subpatches
      p => Mesh%subpatch%begin()
      ipatch = 1
      DO WHILE (ASSOCIATED(p))
         CALL p%get_field(Field,field_up,info)
         or_fail("get_field failed for this subpatch")

#if  __DIME == __2D
            !-------------------------------------------------------------------
            !  --- 2D ---
            !-------------------------------------------------------------------
#if  __MODE == __SCA
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(ipatch)
               iq    = list_sub(ipatch,ip)

               x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
               x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1

               ip10 = FLOOR(x0(1)) + 1
               ip20 = FLOOR(x0(2)) + 1

               ip11 = ip10 + 1
               ip21 = ip20 + 1

               xp1 = x0(1)-REAL(ip10-1,MK)
               xp2 = x0(2)-REAL(ip20-1,MK)

               x10 = xp1
               x11 = x10 - 1.0_MK

               x20 = xp2
               x21 = x20 - 1.0_MK

               a10 = 1.0_MK - x10
               a20 = 1.0_MK - x20

               a11 = 1.0_MK + x11
               a21 = 1.0_MK + x21

               a10a20 = a10*a20
               a10a21 = a10*a21

               a11a20 = a11*a20
               a11a21 = a11*a21

               up(iq) = up(iq) + a10a20*field_up(ip10,ip20)
               up(iq) = up(iq) + a10a21*field_up(ip10,ip21)
               up(iq) = up(iq) + a11a20*field_up(ip11,ip20)
               up(iq) = up(iq) + a11a21*field_up(ip11,ip21)
            END DO  ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)

                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + a10a20*field_up(1,ip10,ip20)
                  up(1,iq) = up(1,iq) + a10a21*field_up(1,ip10,ip21)
                  up(1,iq) = up(1,iq) + a11a20*field_up(1,ip11,ip20)
                  up(1,iq) = up(1,iq) + a11a21*field_up(1,ip11,ip21)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)

                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + a10a20*field_up(1,ip10,ip20)
                  up(1,iq) = up(1,iq) + a10a21*field_up(1,ip10,ip21)
                  up(1,iq) = up(1,iq) + a11a20*field_up(1,ip11,ip20)
                  up(1,iq) = up(1,iq) + a11a21*field_up(1,ip11,ip21)

                  up(2,iq) = up(2,iq) + a10a20*field_up(2,ip10,ip20)
                  up(2,iq) = up(2,iq) + a10a21*field_up(2,ip10,ip21)
                  up(2,iq) = up(2,iq) + a11a20*field_up(2,ip11,ip20)
                  up(2,iq) = up(2,iq) + a11a21*field_up(2,ip11,ip21)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)

                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21

                  up(1,iq) = up(1,iq) + a10a20*field_up(1,ip10,ip20)
                  up(1,iq) = up(1,iq) + a10a21*field_up(1,ip10,ip21)
                  up(1,iq) = up(1,iq) + a11a20*field_up(1,ip11,ip20)
                  up(1,iq) = up(1,iq) + a11a21*field_up(1,ip11,ip21)

                  up(2,iq) = up(2,iq) + a10a20*field_up(2,ip10,ip20)
                  up(2,iq) = up(2,iq) + a10a21*field_up(2,ip10,ip21)
                  up(2,iq) = up(2,iq) + a11a20*field_up(2,ip11,ip20)
                  up(2,iq) = up(2,iq) + a11a21*field_up(2,ip11,ip21)

                  up(3,iq) = up(3,iq) + a10a20*field_up(3,ip10,ip20)
                  up(3,iq) = up(3,iq) + a10a21*field_up(3,ip10,ip21)
                  up(3,iq) = up(3,iq) + a11a20*field_up(3,ip11,ip20)
                  up(3,iq) = up(3,iq) + a11a21*field_up(3,ip11,ip21)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)

                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21

                  a10a20 = a10*a20
                  a10a21 = a10*a21

                  a11a20 = a11*a20
                  a11a21 = a11*a21
                  DO ldn=1,lda
                     up(ldn,iq) = up(ldn,iq) + a10a20*field_up(ldn,ip10,ip20)
                     up(ldn,iq) = up(ldn,iq) + a10a21*field_up(ldn,ip10,ip21)
                     up(ldn,iq) = up(ldn,iq) + a11a20*field_up(ldn,ip11,ip20)
                     up(ldn,iq) = up(ldn,iq) + a11a21*field_up(ldn,ip11,ip21)
                  END DO   ! ldn
               END DO ! end loop over particles in the current subdomain
            END IF
#endif
#elif __DIME == __3D
            !-------------------------------------------------------------------
            !  --- 3D ---
            !-------------------------------------------------------------------
#if   __MODE == __SCA
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(ipatch)
               iq    = list_sub(ipatch,ip)
               x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
               x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1
               x0(3) = xp(3,iq)*dxi(3)-p%istart(3) + 1

               ip10 = FLOOR(x0(1)) + 1
               ip20 = FLOOR(x0(2)) + 1
               ip30 = FLOOR(x0(3)) + 1

               ip11 = ip10 + 1
               ip21 = ip20 + 1
               ip31 = ip30 + 1

               xp1 = x0(1)-REAL(ip10-1,MK)
               xp2 = x0(2)-REAL(ip20-1,MK)
               xp3 = x0(3)-REAL(ip30-1,MK)

               x10 = xp1
               x11 = x10 - 1.0_MK

               x20 = xp2
               x21 = x20 - 1.0_MK

               x30 = xp3
               x31 = x30 - 1.0_MK

               a10 = 1.0_MK - x10
               a20 = 1.0_MK - x20
               a30 = 1.0_MK - x30

               a11 = 1.0_MK + x11
               a21 = 1.0_MK + x21
               a31 = 1.0_MK + x31

               a10a20a30 = a10*a20*a30
               a10a20a31 = a10*a20*a31

               a10a21a30 = a10*a21*a30
               a10a21a31 = a10*a21*a31

               a11a20a30 = a11*a20*a30
               a11a20a31 = a11*a20*a31

               a11a21a30 = a11*a21*a30
               a11a21a31 = a11*a21*a31

               up(iq) = up(iq) + &
     &                     a10a20a30*field_up(ip10,ip20,ip30)
               up(iq) = up(iq) + &
     &                     a10a20a31*field_up(ip10,ip20,ip31)

               up(iq) = up(iq) + &
     &                     a10a21a30*field_up(ip10,ip21,ip30)
               up(iq) = up(iq) + &
     &                     a10a21a31*field_up(ip10,ip21,ip31)

               up(iq) = up(iq) + &
     &                     a11a20a30*field_up(ip11,ip20,ip30)
               up(iq) = up(iq) + &
     &                     a11a20a31*field_up(ip11,ip20,ip31)

               up(iq) = up(iq) + &
     &                     a11a21a30*field_up(ip11,ip21,ip30)
               up(iq) = up(iq) + &
     &                     a11a21a31*field_up(ip11,ip21,ip31)
            END DO  ! iq
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)
                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1
                  x0(3) = xp(3,iq)*dxi(3)-p%istart(3) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1
                  ip30 = FLOOR(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)
                  xp3 = x0(3)-REAL(ip30-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  x30 = xp3
                  x31 = x30 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20
                  a30 = 1.0_MK - x30

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21
                  a31 = 1.0_MK + x31

                  a10a20a30 = a10*a20*a30
                  a10a20a31 = a10*a20*a31

                  a10a21a30 = a10*a21*a30
                  a10a21a31 = a10*a21*a31

                  a11a20a30 = a11*a20*a30
                  a11a20a31 = a11*a20*a31

                  a11a21a30 = a11*a21*a30
                  a11a21a31 = a11*a21*a31

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)
                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1
                  x0(3) = xp(3,iq)*dxi(3)-p%istart(3) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1
                  ip30 = FLOOR(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)
                  xp3 = x0(3)-REAL(ip30-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  x30 = xp3
                  x31 = x30 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20
                  a30 = 1.0_MK - x30

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21
                  a31 = 1.0_MK + x31

                  a10a20a30 = a10*a20*a30
                  a10a20a31 = a10*a20*a31

                  a10a21a30 = a10*a21*a30
                  a10a21a31 = a10*a21*a31

                  a11a20a30 = a11*a20*a30
                  a11a20a31 = a11*a20*a31

                  a11a21a30 = a11*a21*a30
                  a11a21a31 = a11*a21*a31

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31)
               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)
                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1
                  x0(3) = xp(3,iq)*dxi(3)-p%istart(3) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1
                  ip30 = FLOOR(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)
                  xp3 = x0(3)-REAL(ip30-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  x30 = xp3
                  x31 = x30 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20
                  a30 = 1.0_MK - x30

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21
                  a31 = 1.0_MK + x31

                  a10a20a30 = a10*a20*a30
                  a10a20a31 = a10*a20*a31

                  a10a21a30 = a10*a21*a30
                  a10a21a31 = a10*a21*a31

                  a11a20a30 = a11*a20*a30
                  a11a20a31 = a11*a20*a31

                  a11a21a30 = a11*a21*a30
                  a11a21a31 = a11*a21*a31

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31)

                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31)

                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31)

                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31)

                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31)

                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31)
               END DO ! end loop over particles in the current subdomain

               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(ipatch)
                  iq    = list_sub(ipatch,ip)
                  x0(1) = xp(1,iq)*dxi(1)-p%istart(1) + 1
                  x0(2) = xp(2,iq)*dxi(2)-p%istart(2) + 1
                  x0(3) = xp(3,iq)*dxi(3)-p%istart(3) + 1

                  ip10 = FLOOR(x0(1)) + 1
                  ip20 = FLOOR(x0(2)) + 1
                  ip30 = FLOOR(x0(3)) + 1

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  xp1 = x0(1)-REAL(ip10-1,MK)
                  xp2 = x0(2)-REAL(ip20-1,MK)
                  xp3 = x0(3)-REAL(ip30-1,MK)

                  x10 = xp1
                  x11 = x10 - 1.0_MK

                  x20 = xp2
                  x21 = x20 - 1.0_MK

                  x30 = xp3
                  x31 = x30 - 1.0_MK

                  a10 = 1.0_MK - x10
                  a20 = 1.0_MK - x20
                  a30 = 1.0_MK - x30

                  a11 = 1.0_MK + x11
                  a21 = 1.0_MK + x21
                  a31 = 1.0_MK + x31

                  a10a20a30 = a10*a20*a30
                  a10a20a31 = a10*a20*a31

                  a10a21a30 = a10*a21*a30
                  a10a21a31 = a10*a21*a31

                  a11a20a30 = a11*a20*a30
                  a11a20a31 = a11*a20*a31

                  a11a21a30 = a11*a21*a30
                  a11a21a31 = a11*a21*a31

                  DO ldn=1,lda
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a30*field_up(ldn,ip10,ip20,ip30)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a31*field_up(ldn,ip10,ip20,ip31)

                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a30*field_up(ldn,ip10,ip21,ip30)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a31*field_up(ldn,ip10,ip21,ip31)

                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a30*field_up(ldn,ip11,ip20,ip30)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a31*field_up(ldn,ip11,ip20,ip31)

                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a30*field_up(ldn,ip11,ip21,ip30)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a31*field_up(ldn,ip11,ip21,ip31)
                  END DO   ! ldn
               END DO ! end loop over particles in the current subdomain
            END IF ! lda unroll
#endif
#endif
          p => Mesh%subpatch%next()
          ipatch = ipatch + 1
      ENDDO !subpatch


      end_subroutine()

#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_bsp2_dv_3d
#endif
#endif
#endif
