      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_interp_m2p
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
program ppm_test_interp_m2p

    USE ppm_module_core
    USE ppm_module_typedef
    USE ppm_module_mktopo
    USE ppm_module_topo_get
    USE ppm_module_interp_m2p
    USE ppm_module_interp_p2m
    USE ppm_module_init
    USE ppm_module_finalize
    USE ppm_module_map
    USE ppm_module_data_rmsh
    USE ppm_module_data_mesh
    USE ppm_module_mpi
    IMPLICIT NONE

#include "../../ppm_define.h"

    INTEGER, PARAMETER              :: debug = 0
    INTEGER, PARAMETER              :: MK = ppm_kind_double
    REAL(MK),parameter              :: pi = 3.1415926535897931_MK
    INTEGER                         :: ndim,nspec
    INTEGER                         :: decomp
    INTEGER                         :: assig
    INTEGER                         :: tolexp
    REAL(MK)                        :: tol
    INTEGER                         :: info
    INTEGER                         :: topoid,meshid
    INTEGER,  parameter             :: ngrid = 2383
    INTEGER,  parameter             :: npgrid = 3029
    REAL(MK),DIMENSION(:,:),POINTER :: xp,wp
    REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,len_phys,h,p_h
    INTEGER, DIMENSION(:  ),POINTER :: ghostsize
    INTEGER                         :: i,j,ai,aj,p_i
    INTEGER, DIMENSION(6)           :: bcdef
    REAL(MK),DIMENSION(:  ),POINTER :: cost
    INTEGER, DIMENSION(:,:),POINTER :: istart,ndata
    INTEGER, DIMENSION(:  ),POINTER :: nm
    REAL(MK),DIMENSION(:,:),POINTER :: minsub,maxsub
    INTEGER, DIMENSION(:  ),POINTER :: sub2proc
    INTEGER                         :: np,mp
    INTEGER, PARAMETER              :: kernel = ppm_param_rmsh_kernel_mp4
    REAL(MK),DIMENSION(:,:,:,:  ), POINTER :: field_wp ! 2d  field_up(ldn,i,j,isub)
    !REAL(MK),DIMENSION(:,:,:,:,:), POINTER :: field_up ! 3d  field_up(ldn,i,j,k,isub)
    REAL(MK),DIMENSION(:  ),POINTER :: field_x
    REAL(MK),DIMENSION(2)           :: x
    TYPE(ppm_t_topo), POINTER       :: topo
    REAL(MK)                        :: maxm3
    INTEGER, PARAMETER              :: nmom = 10
    INTEGER, DIMENSION(2,nmom)      :: alpha
    REAL(MK),DIMENSION(nmom)        :: f_moments, p_moments
    REAL(MK),DIMENSION(2)           :: sigma,mu
    !----------------
    ! setup
    !----------------
    tol = 10.0_MK**-6 !100000.0_MK*EPSILON(1.0_MK)
    tolexp = INT(LOG10(tol))
    ndim = 2
    nspec = 1
    data ((alpha(ai,aj), ai=1,2), aj=1,nmom) /0,0, 1,0, 0,1, 2,0, 0,2, &
   &                                          1,1, 3,0, 0,3, 2,1, 1,2/

    ALLOCATE(len_phys(ndim),min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),p_h(ndim),field_x(ndim),STAT=info)

    DO i=1,ndim
        min_phys(i) = 0.0_MK
        max_phys(i) = 10.0_MK
        ghostsize(i) = 2
        len_phys(i) = max_phys(i) - min_phys(i)
    ENDDO
    bcdef(1:6) = ppm_param_bcdef_periodic

    NULLIFY(xp)
    NULLIFY(wp)

#ifdef __MPI
    CALL MPI_Init(info)
#endif

    CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = npgrid**2
    mp = 0

    ALLOCATE(xp(ndim,np),wp(nspec,np),STAT=info)
    xp = 0.0_MK
    wp = 0.0_MK


    ALLOCATE(nm(ndim),STAT=info)
    DO i=1,ndim
        nm(i) = ngrid
    ENDDO

    DO i=1,ndim
        p_h(i) = (max_phys(i) - min_phys(i)) / REAL(npgrid,mk)
    ENDDO

    DO j=1,npgrid
        DO i=1,npgrid
            p_i = i + (j-1)*npgrid
            xp(1,p_i) = REAL(i-1,mk)*p_h(1)
            xp(2,p_i) = REAL(j-1,mk)*p_h(2)
        ENDDO
    ENDDO

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_user_defined
    assig  = ppm_param_assign_user_defined

    topoid = 0
    meshid = -1

    ALLOCATE(minsub(ndim,1),maxsub(ndim,1),sub2proc(1),cost(1),STAT=info)

    cost(1) = 1.0_MK

    DO i=1,ndim
        minsub(i,1)=min_phys(i)
        maxsub(i,1)=max_phys(i)
    ENDDO

    sub2proc(1) = 0

    CALL ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=1,user_sub2proc=sub2proc)

    ALLOCATE(field_wp(nspec,(1-ghostsize(1)):(ndata(1,1)+ghostsize(1)),         &
    &        (1-ghostsize(2)):(ndata(2,1)+ghostsize(2)),1),STAT=info) ! 2d
    !ALLOCATE(field_wp(nspec,ndata(1,1),ndata(2,1),ndata(3,1),1),STAT=info) ! 3d

    field_wp = 0.0_MK

    DO i=1,ndim
        h(i) = (max_phys(i) - min_phys(i)) / REAL(ndata(i,1)-1,mk)
    ENDDO
    !----------------
    ! setup mesh data
    !----------------

    DO i=1,ndim
        mu(i) = 0.0_MK !(max_phys(i) - min_phys(i))/2.0_MK
    ENDDO
    sigma = 2.5_MK

    DO j=1,ndata(2,1)
        DO i=1,ndata(1,1)
            field_x(1) = min_phys(1) + h(1)*REAL(i-1,mk)
            field_x(2) = min_phys(2) + h(2)*REAL(j-1,mk)
            field_wp(1,i,j,1) = 0.0001_MK/(2.0_MK*pi*sigma(1)*sigma(2))*&
   &                 exp(-0.5_MK*(((field_x(1)-mu(1))**2/sigma(1)**2)+  &
   &                              ((field_x(2)-mu(2))**2/sigma(2)**2)))
        ENDDO
    ENDDO

    print *, 'finished initialization'

   CALL ppm_map_field_ghost_get(topoid,meshid,ghostsize,info)
   CALL ppm_map_field_push(topoid,meshid,field_wp,1,info)
   CALL ppm_map_field_send(info)
   CALL ppm_map_field_pop(topoid,meshid,field_wp,1,ghostsize,info)

    !----------------
    ! m --> p
    !----------------
    CALL ppm_interp_m2p(topoid,meshid,xp,np,wp,1,kernel,ghostsize,field_wp,info)

    !----------------
    ! test m --> p
    !----------------
    f_moments = 0.0_MK
    p_moments = 0.0_MK
    do p_i = 1,np
        x = xp(:,p_i)
        IF (xp(1,p_i).GE.len_phys(1)) x(1)  = xp(1,p_i) - len_phys(1)
        IF (xp(2,p_i).GE.len_phys(2)) x(2)  = xp(2,p_i) - len_phys(2)
        do aj = 2,nmom
            p_moments(aj) = p_moments(aj) + &
 &                      p_h(1)*p_h(2)*wp(1,p_i)*x(1)**alpha(1,aj)* &
 &                      x(2)**alpha(2,aj)
        ENDDO
        p_moments(1) = p_moments(1) + wp(1,p_i)*p_h(1)*p_h(2)
    ENDDO
    DO j=1,ndata(2,1)-1
        DO i=1,ndata(1,1)-1
               field_x(1) = min_phys(1) + h(1)*REAL(i-1,mk)
               field_x(2) = min_phys(2) + h(2)*REAL(j-1,mk)
               x = field_x
               IF (field_x(1).GE.len_phys(1)) x(1)  = field_x(1) - len_phys(1)
               IF (field_x(2).GE.len_phys(2)) x(2)  = field_x(2) - len_phys(2)
               do aj = 2,nmom
                   f_moments(aj) = f_moments(aj) + &
 &                      h(1)*h(2)*field_wp(1,i,j,1)*x(1)**alpha(1,aj)* &
 &                      x(2)**alpha(2,aj)
            ENDDO
            f_moments(1) = f_moments(1) + field_wp(1,i,j,1)*h(1)*h(2)
        ENDDO
    ENDDO

    do aj = 1,6
        IF (ABS(p_moments(aj) - f_moments(aj)) .GT. tol) then
            print *, 'failed at moment: ', aj
            print *, 'particle moments: ',p_moments
            print *, 'field moments: ',   f_moments
            stop 'ERROR: m2p interpolation: moments not conserved.'
        ENDIF
    ENDDO


    !----------------
    ! teardown...
    !----------------
    print *, 'cleaning up'
8000 CONTINUE

    CALL ppm_finalize(info)

#ifdef __MPI
    CALL MPI_Finalize(info)
#endif

    DEALLOCATE(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

END PROGRAM ppm_test_interp_m2p
