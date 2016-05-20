      !-------------------------------------------------------------------------
      !     Subroutine   :                   m2p_interp_mp4
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
      SUBROUTINE m2p_interp_mp4_ss_2d(topoid,meshid,field_up,xp,up,dx,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_mp4_ds_2d(topoid,meshid,field_up,xp,up,dx,ghostsize,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_mp4_sv_2d(topoid,meshid,field_up,lda,xp,up,dx,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_mp4_dv_2d(topoid,meshid,field_up,lda,xp,up,dx,ghostsize,info)
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_mp4_ss_3d(topoid,meshid,field_up,xp,up,dx,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_mp4_ds_3d(topoid,meshid,field_up,xp,up,dx,ghostsize,info)
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE m2p_interp_mp4_sv_3d(topoid,meshid,field_up,lda,xp,up,dx,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE m2p_interp_mp4_dv_3d(topoid,meshid,field_up,lda,xp,up,dx,ghostsize,info)
#endif
#endif
#endif
     !!! Mesh to particle interpolation following the MP4 scheme.
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
      USE ppm_module_data
      USE ppm_module_data_rmsh
      USE ppm_module_data_mesh
      USE ppm_module_typedef
      USE ppm_module_write
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      ! Arguments
      !-------------------------------------------------------------------------
      INTEGER                        , INTENT(IN   ) :: topoid
      !!! topology identifier of target
      INTEGER                        , INTENT(IN   ) :: meshid
      !!! id of the mesh (user)
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
      INTEGER                        , INTENT(  OUT) :: info
      !!! Returns 0 upon success
      !-------------------------------------------------------------------------
      ! Local variables
      !-------------------------------------------------------------------------
      INTEGER,  DIMENSION(:,:)     , POINTER :: istart
      REAL(mk), DIMENSION(:)       , POINTER :: min_phys
      REAL(mk), DIMENSION(:)       , POINTER :: max_phys
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
      INTEGER                                :: isub,ifrom,ito,ip,iopt,isubl
      INTEGER                                :: max_partnumber,idom,nlist2,idoml
      INTEGER, DIMENSION(ppm_dim)            :: Nm
      INTEGER                                :: nsubs
      INTEGER, DIMENSION(6)                  :: bcdef
      INTEGER                                :: iq
      LOGICAL                                :: internal_weights,lok
      ! aliases
      REAL(mk), DIMENSION(:,:),      POINTER :: min_sub
      REAL(mk), DIMENSION(:,:),      POINTER :: max_sub
      REAL(mk)                               :: myeps
      REAL(mk)                               :: tim1s, tim1e
      REAL(mk)                               :: xp1,xp2,xp3
      REAL(mk)                               :: wx1,wx2,wx3
      INTEGER                                :: ldn
      REAL(mk), DIMENSION(ppm_dim)           :: x0
      TYPE(ppm_t_equi_mesh), POINTER         :: mesh
      TYPE(ppm_t_topo)     , POINTER         :: topo
      !-------------------------------------------------------------------------
      !  Variables for unrolled versions
      !-------------------------------------------------------------------------


#if   __DIME == __2D
      REAL(mk) :: a10a20
      REAL(mk) :: a10a21
      REAL(mk) :: a11a20
      REAL(mk) :: a11a21
#elif __DIME == __3D
      REAL(mk) :: x10,x11,x12,x13,x20,x21,x22,x23,x30,x31,x32,x33
      REAL(mk) :: a10,a11,a12,a13,a20,a21,a22,a23,a30,a31,a32,a33
      INTEGER  :: ip10,ip11,ip12,ip13,ip20,ip21,ip22,ip23,ip30,ip31,ip32,ip33
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
#endif

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('m2p_interp_mp4',t0,info)

      topo => ppm_topo(topoid)%t
      mesh => topo%mesh(meshid)

      istart => mesh%istart

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


      dxi = 1.0_mk/dx

         !  loop over subs
         DO isub = 1,topo%nsublist
#if __DIME == __2D
            !-------------------------------------------------------------------
            !  --- 2D ---
            !-------------------------------------------------------------------
#if __MODE == __SCA
            isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               x0(1) = (xp(1,iq)-min_phys(1))*dxi(1)
               x0(2) = (xp(2,iq)-min_phys(2))*dxi(2)

               ip1 = INT(x0(1))+2-istart(1,isubl)
               ip2 = INT(x0(2))+2-istart(2,isubl)

               xp1 = x0(1)-AINT(x0(1))
               xp2 = x0(2)-AINT(x0(2))

               DO jj = -1,2
                  x2 = ABS(xp2 - REAL(jj,mk))
                  IF(x2.LT.1.0_mk) THEN
                     wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                  ELSE
                     wx2 = 2.0_mk + (-4.0_mk + &
     &                          (2.5_mk - 0.5_mk * x2)*x2)*x2
                  END IF

                  DO ii    = - 1,2
                     x1 = ABS(xp1 - REAL(ii,mk))
                     IF(x1.LT.1.0_MK) THEN
                        wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                              1.5_mk*x1)
                     ELSE
                        wx1 =  2.0_mk + (-4.0_mk + &
     &                              (2.5_mk - 0.5_mk*x1)*x1)*x1
                     END IF
                     up(iq) = up(iq) + wx1*wx2*field_up(ii+ip1,jj+ip2,isub)
                  END DO
               END DO
            END DO ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            ! This will only vectorize over lda
            isubl = topo%isublist(isub)
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)
               x0(1) = (xp(1,iq)-min_phys(1))*dxi(1)
               x0(2) = (xp(2,iq)-min_phys(2))*dxi(2)

               ip1 = INT(x0(1))+2-istart(1,isubl)
               ip2 = INT(x0(2))+2-istart(2,isubl)

               xp1 = x0(1)-AINT(x0(1))
               xp2 = x0(2)-AINT(x0(2))

               DO jj = -1,2
                  x2 = ABS(xp2 - REAL(jj,mk))
                  IF(x2.LT.1.0_mk) THEN
                     wx2 = 1.0_mk - x2**2*(2.5_mk-1.5_mk*x2)
                  ELSE
                     wx2 = 2.0_mk + (-4.0_mk + &
     &                          (2.5_mk - 0.5_mk * x2)*x2)*x2
                  END IF

                  DO ii    = - 1,2
                     x1 = ABS(xp1 - REAL(ii,mk))
                     IF(x1.LT.1.0_MK) THEN
                        wx1 =  1.0_mk - x1**2*(2.5_mk - &
     &                              1.5_mk*x1)
                     ELSE
                        wx1 =  2.0_mk + (-4.0_mk + &
     &                              (2.5_mk - 0.5_mk*x1)*x1)*x1
                     END IF
                     DO ldn=1,lda
                        up(ldn,iq) = up(ldn,iq) + wx1*wx2*     &
     &                                                   field_up(ldn,ii+ip1,jj+ip2,isub)
                     END DO
                  END DO
               END DO
            END DO ! end loop over particles in the current subdomain
#endif
#elif __DIME == __3D
            !-------------------------------------------------------------------
            !  --- 3D ---
            !-------------------------------------------------------------------
#if   __MODE == __SCA
            isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
            DO ip = 1,store_info(isub)
               iq    = list_sub(isub,ip)

               x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
               x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
               x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

               ip10 = INT(x0(1))
               ip20 = INT(x0(2))
               ip30 = INT(x0(3))

               ip11 = ip10 + 1
               ip21 = ip20 + 1
               ip31 = ip30 + 1

               ip12 = ip11 + 1
               ip22 = ip21 + 1
               ip32 = ip31 + 1

               ip13 = ip11 + 2
               ip23 = ip21 + 2
               ip33 = ip31 + 2

               xp1 = x0(1)-REAL(ip10,mk)
               xp2 = x0(2)-REAL(ip20,mk)
               xp3 = x0(3)-REAL(ip30,mk)

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

               up(iq) = up(iq) + &
     &                     a10a20a30*field_up(ip10,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a20a31*field_up(ip10,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a20a32*field_up(ip10,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a20a33*field_up(ip10,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a21a30*field_up(ip10,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a21a31*field_up(ip10,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a21a32*field_up(ip10,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a21a33*field_up(ip10,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a22a30*field_up(ip10,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a22a31*field_up(ip10,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a22a32*field_up(ip10,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a22a33*field_up(ip10,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a10a23a30*field_up(ip10,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a10a23a31*field_up(ip10,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a10a23a32*field_up(ip10,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a10a23a33*field_up(ip10,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a20a30*field_up(ip11,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a20a31*field_up(ip11,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a20a32*field_up(ip11,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a20a33*field_up(ip11,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a21a30*field_up(ip11,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a21a31*field_up(ip11,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a21a32*field_up(ip11,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a21a33*field_up(ip11,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a22a30*field_up(ip11,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a22a31*field_up(ip11,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a22a32*field_up(ip11,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a22a33*field_up(ip11,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a11a23a30*field_up(ip11,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a11a23a31*field_up(ip11,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a11a23a32*field_up(ip11,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a11a23a33*field_up(ip11,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a20a30*field_up(ip12,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a20a31*field_up(ip12,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a20a32*field_up(ip12,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a20a33*field_up(ip12,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a21a30*field_up(ip12,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a21a31*field_up(ip12,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a21a32*field_up(ip12,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a21a33*field_up(ip12,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a22a30*field_up(ip12,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a22a31*field_up(ip12,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a22a32*field_up(ip12,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a22a33*field_up(ip12,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a12a23a30*field_up(ip12,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a12a23a31*field_up(ip12,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a12a23a32*field_up(ip12,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a12a23a33*field_up(ip12,ip23,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a20a30*field_up(ip13,ip20,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a20a31*field_up(ip13,ip20,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a20a32*field_up(ip13,ip20,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a20a33*field_up(ip13,ip20,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a21a30*field_up(ip13,ip21,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a21a31*field_up(ip13,ip21,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a21a32*field_up(ip13,ip21,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a21a33*field_up(ip13,ip21,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a22a30*field_up(ip13,ip22,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a22a31*field_up(ip13,ip22,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a22a32*field_up(ip13,ip22,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a22a33*field_up(ip13,ip22,ip33,isub)
               up(iq) = up(iq) + &
     &                     a13a23a30*field_up(ip13,ip23,ip30,isub)
               up(iq) = up(iq) + &
     &                     a13a23a31*field_up(ip13,ip23,ip31,isub)
               up(iq) = up(iq) + &
     &                     a13a23a32*field_up(ip13,ip23,ip32,isub)
               up(iq) = up(iq) + &
     &                     a13a23a33*field_up(ip13,ip23,ip33,isub)
            END DO ! end loop over particles in the current subdomain
#elif __MODE == __VEC
            !-------------------------------------------------------------------
            !  Unrolled version for 1-vectors
            !-------------------------------------------------------------------
            IF(lda.EQ.1) THEN
               isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 2-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.2) THEN
               isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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

                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a32*field_up(2,ip10,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a33*field_up(2,ip10,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a32*field_up(2,ip10,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a33*field_up(2,ip10,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a30*field_up(2,ip10,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a31*field_up(2,ip10,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a32*field_up(2,ip10,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a33*field_up(2,ip10,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a30*field_up(2,ip10,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a31*field_up(2,ip10,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a32*field_up(2,ip10,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a33*field_up(2,ip10,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a32*field_up(2,ip11,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a33*field_up(2,ip11,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a32*field_up(2,ip11,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a33*field_up(2,ip11,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a30*field_up(2,ip11,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a31*field_up(2,ip11,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a32*field_up(2,ip11,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a33*field_up(2,ip11,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a30*field_up(2,ip11,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a31*field_up(2,ip11,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a32*field_up(2,ip11,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a33*field_up(2,ip11,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a30*field_up(2,ip12,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a31*field_up(2,ip12,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a32*field_up(2,ip12,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a33*field_up(2,ip12,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a30*field_up(2,ip12,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a31*field_up(2,ip12,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a32*field_up(2,ip12,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a33*field_up(2,ip12,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a30*field_up(2,ip12,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a31*field_up(2,ip12,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a32*field_up(2,ip12,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a33*field_up(2,ip12,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a30*field_up(2,ip12,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a31*field_up(2,ip12,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a32*field_up(2,ip12,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a33*field_up(2,ip12,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a30*field_up(2,ip13,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a31*field_up(2,ip13,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a32*field_up(2,ip13,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a33*field_up(2,ip13,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a30*field_up(2,ip13,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a31*field_up(2,ip13,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a32*field_up(2,ip13,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a33*field_up(2,ip13,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a30*field_up(2,ip13,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a31*field_up(2,ip13,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a32*field_up(2,ip13,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a33*field_up(2,ip13,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a30*field_up(2,ip13,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a31*field_up(2,ip13,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a32*field_up(2,ip13,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a33*field_up(2,ip13,ip23,ip33,isub)


               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  Unrolled version for 3-vectors
               !----------------------------------------------------------------
            ELSEIF(lda.EQ.3) THEN
               isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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


                  up(1,iq) = up(1,iq) + &
     &                        a10a20a30*field_up(1,ip10,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a31*field_up(1,ip10,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a32*field_up(1,ip10,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a20a33*field_up(1,ip10,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a30*field_up(1,ip10,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a31*field_up(1,ip10,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a32*field_up(1,ip10,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a21a33*field_up(1,ip10,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a30*field_up(1,ip10,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a31*field_up(1,ip10,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a32*field_up(1,ip10,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a22a33*field_up(1,ip10,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a30*field_up(1,ip10,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a31*field_up(1,ip10,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a32*field_up(1,ip10,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a10a23a33*field_up(1,ip10,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a30*field_up(1,ip11,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a31*field_up(1,ip11,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a32*field_up(1,ip11,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a20a33*field_up(1,ip11,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a30*field_up(1,ip11,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a31*field_up(1,ip11,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a32*field_up(1,ip11,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a21a33*field_up(1,ip11,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a30*field_up(1,ip11,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a31*field_up(1,ip11,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a32*field_up(1,ip11,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a22a33*field_up(1,ip11,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a30*field_up(1,ip11,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a31*field_up(1,ip11,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a32*field_up(1,ip11,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a11a23a33*field_up(1,ip11,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a30*field_up(1,ip12,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a31*field_up(1,ip12,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a32*field_up(1,ip12,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a20a33*field_up(1,ip12,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a30*field_up(1,ip12,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a31*field_up(1,ip12,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a32*field_up(1,ip12,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a21a33*field_up(1,ip12,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a30*field_up(1,ip12,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a31*field_up(1,ip12,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a32*field_up(1,ip12,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a22a33*field_up(1,ip12,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a30*field_up(1,ip12,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a31*field_up(1,ip12,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a32*field_up(1,ip12,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a12a23a33*field_up(1,ip12,ip23,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a30*field_up(1,ip13,ip20,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a31*field_up(1,ip13,ip20,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a32*field_up(1,ip13,ip20,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a20a33*field_up(1,ip13,ip20,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a30*field_up(1,ip13,ip21,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a31*field_up(1,ip13,ip21,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a32*field_up(1,ip13,ip21,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a21a33*field_up(1,ip13,ip21,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a30*field_up(1,ip13,ip22,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a31*field_up(1,ip13,ip22,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a32*field_up(1,ip13,ip22,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a22a33*field_up(1,ip13,ip22,ip33,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a30*field_up(1,ip13,ip23,ip30,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a31*field_up(1,ip13,ip23,ip31,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a32*field_up(1,ip13,ip23,ip32,isub)
                  up(1,iq) = up(1,iq) + &
     &                        a13a23a33*field_up(1,ip13,ip23,ip33,isub)

                  up(2,iq) = up(2,iq) + &
     &                        a10a20a30*field_up(2,ip10,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a31*field_up(2,ip10,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a32*field_up(2,ip10,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a20a33*field_up(2,ip10,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a30*field_up(2,ip10,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a31*field_up(2,ip10,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a32*field_up(2,ip10,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a21a33*field_up(2,ip10,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a30*field_up(2,ip10,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a31*field_up(2,ip10,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a32*field_up(2,ip10,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a22a33*field_up(2,ip10,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a30*field_up(2,ip10,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a31*field_up(2,ip10,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a32*field_up(2,ip10,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a10a23a33*field_up(2,ip10,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a30*field_up(2,ip11,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a31*field_up(2,ip11,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a32*field_up(2,ip11,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a20a33*field_up(2,ip11,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a30*field_up(2,ip11,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a31*field_up(2,ip11,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a32*field_up(2,ip11,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a21a33*field_up(2,ip11,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a30*field_up(2,ip11,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a31*field_up(2,ip11,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a32*field_up(2,ip11,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a22a33*field_up(2,ip11,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a30*field_up(2,ip11,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a31*field_up(2,ip11,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a32*field_up(2,ip11,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a11a23a33*field_up(2,ip11,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a30*field_up(2,ip12,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a31*field_up(2,ip12,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a32*field_up(2,ip12,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a20a33*field_up(2,ip12,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a30*field_up(2,ip12,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a31*field_up(2,ip12,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a32*field_up(2,ip12,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a21a33*field_up(2,ip12,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a30*field_up(2,ip12,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a31*field_up(2,ip12,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a32*field_up(2,ip12,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a22a33*field_up(2,ip12,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a30*field_up(2,ip12,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a31*field_up(2,ip12,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a32*field_up(2,ip12,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a12a23a33*field_up(2,ip12,ip23,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a30*field_up(2,ip13,ip20,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a31*field_up(2,ip13,ip20,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a32*field_up(2,ip13,ip20,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a20a33*field_up(2,ip13,ip20,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a30*field_up(2,ip13,ip21,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a31*field_up(2,ip13,ip21,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a32*field_up(2,ip13,ip21,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a21a33*field_up(2,ip13,ip21,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a30*field_up(2,ip13,ip22,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a31*field_up(2,ip13,ip22,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a32*field_up(2,ip13,ip22,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a22a33*field_up(2,ip13,ip22,ip33,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a30*field_up(2,ip13,ip23,ip30,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a31*field_up(2,ip13,ip23,ip31,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a32*field_up(2,ip13,ip23,ip32,isub)
                  up(2,iq) = up(2,iq) + &
     &                        a13a23a33*field_up(2,ip13,ip23,ip33,isub)

                  up(3,iq) = up(3,iq) + &
     &                        a10a20a30*field_up(3,ip10,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a31*field_up(3,ip10,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a32*field_up(3,ip10,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a20a33*field_up(3,ip10,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a30*field_up(3,ip10,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a31*field_up(3,ip10,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a32*field_up(3,ip10,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a21a33*field_up(3,ip10,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a30*field_up(3,ip10,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a31*field_up(3,ip10,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a32*field_up(3,ip10,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a22a33*field_up(3,ip10,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a30*field_up(3,ip10,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a31*field_up(3,ip10,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a32*field_up(3,ip10,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a10a23a33*field_up(3,ip10,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a30*field_up(3,ip11,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a31*field_up(3,ip11,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a32*field_up(3,ip11,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a20a33*field_up(3,ip11,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a30*field_up(3,ip11,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a31*field_up(3,ip11,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a32*field_up(3,ip11,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a21a33*field_up(3,ip11,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a30*field_up(3,ip11,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a31*field_up(3,ip11,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a32*field_up(3,ip11,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a22a33*field_up(3,ip11,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a30*field_up(3,ip11,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a31*field_up(3,ip11,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a32*field_up(3,ip11,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a11a23a33*field_up(3,ip11,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a30*field_up(3,ip12,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a31*field_up(3,ip12,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a32*field_up(3,ip12,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a20a33*field_up(3,ip12,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a30*field_up(3,ip12,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a31*field_up(3,ip12,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a32*field_up(3,ip12,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a21a33*field_up(3,ip12,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a30*field_up(3,ip12,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a31*field_up(3,ip12,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a32*field_up(3,ip12,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a22a33*field_up(3,ip12,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a30*field_up(3,ip12,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a31*field_up(3,ip12,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a32*field_up(3,ip12,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a12a23a33*field_up(3,ip12,ip23,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a30*field_up(3,ip13,ip20,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a31*field_up(3,ip13,ip20,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a32*field_up(3,ip13,ip20,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a20a33*field_up(3,ip13,ip20,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a30*field_up(3,ip13,ip21,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a31*field_up(3,ip13,ip21,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a32*field_up(3,ip13,ip21,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a21a33*field_up(3,ip13,ip21,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a30*field_up(3,ip13,ip22,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a31*field_up(3,ip13,ip22,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a32*field_up(3,ip13,ip22,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a22a33*field_up(3,ip13,ip22,ip33,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a30*field_up(3,ip13,ip23,ip30,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a31*field_up(3,ip13,ip23,ip31,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a32*field_up(3,ip13,ip23,ip32,isub)
                  up(3,iq) = up(3,iq) + &
     &                        a13a23a33*field_up(3,ip13,ip23,ip33,isub)

               END DO ! end loop over particles in the current subdomain
               !----------------------------------------------------------------
               !  All other lda are not unrolled. This will vectorize over lda!
               !----------------------------------------------------------------
            ELSE
               isubl = topo%isublist(isub)
#ifdef __SXF90
!CDIR NODEP
#endif
               DO ip = 1,store_info(isub)
                  iq    = list_sub(isub,ip)

                  x0(1) = (xp(1,iq)-min_sub(1,isubl))*dxi(1)
                  x0(2) = (xp(2,iq)-min_sub(2,isubl))*dxi(2)
                  x0(3) = (xp(3,iq)-min_sub(3,isubl))*dxi(3)

                  ip10 = INT(x0(1))
                  ip20 = INT(x0(2))
                  ip30 = INT(x0(3))

                  ip11 = ip10 + 1
                  ip21 = ip20 + 1
                  ip31 = ip30 + 1

                  ip12 = ip11 + 1
                  ip22 = ip21 + 1
                  ip32 = ip31 + 1

                  ip13 = ip11 + 2
                  ip23 = ip21 + 2
                  ip33 = ip31 + 2

                  xp1 = x0(1)-REAL(ip10,mk)
                  xp2 = x0(2)-REAL(ip20,mk)
                  xp3 = x0(3)-REAL(ip30,mk)

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
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a30*field_up(ldn,ip10,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a31*field_up(ldn,ip10,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a32*field_up(ldn,ip10,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a20a33*field_up(ldn,ip10,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a30*field_up(ldn,ip10,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a31*field_up(ldn,ip10,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a32*field_up(ldn,ip10,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a21a33*field_up(ldn,ip10,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a30*field_up(ldn,ip10,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a31*field_up(ldn,ip10,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a32*field_up(ldn,ip10,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a22a33*field_up(ldn,ip10,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a30*field_up(ldn,ip10,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a31*field_up(ldn,ip10,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a32*field_up(ldn,ip10,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a10a23a33*field_up(ldn,ip10,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a30*field_up(ldn,ip11,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a31*field_up(ldn,ip11,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a32*field_up(ldn,ip11,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a20a33*field_up(ldn,ip11,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a30*field_up(ldn,ip11,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a31*field_up(ldn,ip11,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a32*field_up(ldn,ip11,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a21a33*field_up(ldn,ip11,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a30*field_up(ldn,ip11,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a31*field_up(ldn,ip11,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a32*field_up(ldn,ip11,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a22a33*field_up(ldn,ip11,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a30*field_up(ldn,ip11,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a31*field_up(ldn,ip11,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a32*field_up(ldn,ip11,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a11a23a33*field_up(ldn,ip11,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a30*field_up(ldn,ip12,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a31*field_up(ldn,ip12,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a32*field_up(ldn,ip12,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a20a33*field_up(ldn,ip12,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a30*field_up(ldn,ip12,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a31*field_up(ldn,ip12,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a32*field_up(ldn,ip12,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a21a33*field_up(ldn,ip12,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a30*field_up(ldn,ip12,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a31*field_up(ldn,ip12,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a32*field_up(ldn,ip12,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a22a33*field_up(ldn,ip12,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a30*field_up(ldn,ip12,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a31*field_up(ldn,ip12,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a32*field_up(ldn,ip12,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a12a23a33*field_up(ldn,ip12,ip23,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a30*field_up(ldn,ip13,ip20,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a31*field_up(ldn,ip13,ip20,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a32*field_up(ldn,ip13,ip20,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a20a33*field_up(ldn,ip13,ip20,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a30*field_up(ldn,ip13,ip21,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a31*field_up(ldn,ip13,ip21,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a32*field_up(ldn,ip13,ip21,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a21a33*field_up(ldn,ip13,ip21,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a30*field_up(ldn,ip13,ip22,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a31*field_up(ldn,ip13,ip22,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a32*field_up(ldn,ip13,ip22,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a22a33*field_up(ldn,ip13,ip22,ip33,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a30*field_up(ldn,ip13,ip23,ip30,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a31*field_up(ldn,ip13,ip23,ip31,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a32*field_up(ldn,ip13,ip23,ip32,isub)
                     up(ldn,iq) = up(ldn,iq) + &
     &                           a13a23a33*field_up(ldn,ip13,ip23,ip33,isub)
                  END DO ! lda
               END DO ! end loop over particles in the current subdomain
            END IF ! unrolled lda cases
#endif
#endif
         END DO ! isub


#if   __DIME == __2D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_ss_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_ds_2d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_sv_2d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_dv_2d
#endif
#endif
#elif __DIME == __3D
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_ss_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_ds_3d
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_sv_3d
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE m2p_interp_mp4_dv_3d
#endif
#endif
#endif
