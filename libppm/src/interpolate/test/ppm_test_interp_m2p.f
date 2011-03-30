      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_interp_m2p
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
program ppm_test_interp_m2p

    use ppm_module_user
    use ppm_module_typedef
    use ppm_module_mktopo
    use ppm_module_topo_get
    use ppm_module_interp_m2p
    use ppm_module_interp_p2m
    use ppm_module_init
    use ppm_module_finalize
    use ppm_module_map
    USE ppm_module_data_rmsh
    USE ppm_module_data_mesh

    implicit none

#include "../../ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif

    integer, parameter              :: debug = 0
    integer, parameter              :: MK = ppm_kind_double
    real(mk),parameter              :: pi = 3.1415926535897931_mk
    integer                         :: ndim,nspec
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    real(MK)                        :: tol
    integer                         :: info
    integer                         :: topoid,meshid
    integer,  parameter             :: ngrid = 2383
    integer,  parameter             :: npgrid = 3029
    real(MK),dimension(:,:),pointer :: xp,wp
    real(MK),dimension(:  ),pointer :: min_phys,max_phys,len_phys,h,p_h
    integer, dimension(:  ),pointer :: ghostsize
    integer                         :: i,j,ai,aj,p_i
    integer, dimension(6)           :: bcdef
    real(MK),dimension(:  ),pointer :: cost
    integer, dimension(:,:),pointer :: istart,ndata
    integer, dimension(:  ),pointer :: nm
    real(MK),dimension(:,:),pointer :: minsub,maxsub
    integer, dimension(:  ),pointer :: sub2proc
    integer                         :: np,mp
    integer, parameter              :: kernel = ppm_param_rmsh_kernel_mp4
    real(MK),dimension(:,:,:,:  ), pointer :: field_wp ! 2d  field_up(ldn,i,j,isub)
    !real(MK),dimension(:,:,:,:,:), pointer :: field_up ! 3d  field_up(ldn,i,j,k,isub)
    real(MK),dimension(:  ),pointer :: field_x
    real(mk),dimension(2)           :: x
    type(ppm_t_topo), pointer       :: topo
    real(MK)                        :: maxm3
    integer, parameter              :: nmom = 10
    integer, dimension(2,nmom)      :: alpha
    real(MK),dimension(nmom)        :: f_moments, p_moments
    real(MK),dimension(2)           :: sigma,mu
    !----------------
    ! setup
    !----------------
    tol = 10.0_mk**-6 !100000.0_mk*epsilon(1.0_mk)
    tolexp = int(log10(tol))
    ndim = 2
    nspec = 1
    data ((alpha(ai,aj), ai=1,2), aj=1,nmom) /0,0, 1,0, 0,1, 2,0, 0,2, &
   &                                          1,1, 3,0, 0,3, 2,1, 1,2/

    allocate(len_phys(ndim),min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),p_h(ndim),field_x(ndim),stat=info)

    do i=1,ndim
        min_phys(i) = 0.0_mk
        max_phys(i) = 10.0_mk
        ghostsize(i) = 2
        len_phys(i) = max_phys(i) - min_phys(i)
    enddo
    bcdef(1:6) = ppm_param_bcdef_periodic

    nullify(xp)
    nullify(wp)

#ifdef __MPI
    call MPI_Init(info)
#endif

    call ppm_init(ndim,MK,tolexp,0,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = npgrid**2
    mp = 0

    allocate(xp(ndim,np),wp(nspec,np),stat=info)
    xp = 0.0_mk
    wp = 0.0_mk


    allocate(nm(ndim),stat=info)
    do i=1,ndim
        nm(i) = ngrid
    enddo

    do i=1,ndim
        p_h(i) = (max_phys(i) - min_phys(i)) / real(npgrid,mk)
    enddo

    do j=1,npgrid
        do i=1,npgrid
            p_i = i + (j-1)*npgrid
            xp(1,p_i) = real(i-1,mk)*p_h(1)
            xp(2,p_i) = real(j-1,mk)*p_h(2)
        enddo
    enddo

    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_user_defined
    assig  = ppm_param_assign_user_defined

    topoid = 0
    meshid = -1

    allocate(minsub(ndim,1),maxsub(ndim,1),sub2proc(1),cost(1),stat=info)

    cost(1) = 1.0_mk

    do i=1,ndim
        minsub(i,1)=min_phys(i)
        maxsub(i,1)=max_phys(i)
    enddo

    sub2proc(1) = 0

    call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=1,user_sub2proc=sub2proc)

    allocate(field_wp(nspec,(1-ghostsize(1)):(ndata(1,1)+ghostsize(1)),         &
    &        (1-ghostsize(2)):(ndata(2,1)+ghostsize(2)),1),stat=info) ! 2d
    !allocate(field_wp(nspec,ndata(1,1),ndata(2,1),ndata(3,1),1),stat=info) ! 3d

    field_wp = 0.0_mk

    do i=1,ndim
        h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
    enddo
    !----------------
    ! setup mesh data
    !----------------

    do i=1,ndim
        mu(i) = 0.0_mk !(max_phys(i) - min_phys(i))/2.0_mk
    enddo
    sigma = 2.5_mk

    do j=1,ndata(2,1)
        do i=1,ndata(1,1)
            field_x(1) = min_phys(1) + h(1)*real(i-1,mk)
            field_x(2) = min_phys(2) + h(2)*real(j-1,mk)
            field_wp(1,i,j,1) = 0.0001_mk/(2.0_mk*pi*sigma(1)*sigma(2))*&
   &                 exp(-0.5_mk*(((field_x(1)-mu(1))**2/sigma(1)**2)+  &
   &                              ((field_x(2)-mu(2))**2/sigma(2)**2)))
        enddo
    enddo

    print *, 'finished initialization'

   call ppm_map_field_ghost_get(topoid,meshid,ghostsize,info)
   call ppm_map_field_push(topoid,meshid,field_wp,1,info)
   call ppm_map_field_send(info)
   call ppm_map_field_pop(topoid,meshid,field_wp,1,ghostsize,info)

    !----------------
    ! m --> p
    !----------------
    call ppm_interp_m2p(topoid,meshid,xp,np,wp,1,kernel,ghostsize,field_wp,info)

    !----------------
    ! test m --> p
    !----------------
    f_moments = 0.0_mk
    p_moments = 0.0_mk
    do p_i = 1,np
        x = xp(:,p_i)
        if (xp(1,p_i).ge.len_phys(1)) x(1)  = xp(1,p_i) - len_phys(1)
        if (xp(2,p_i).ge.len_phys(2)) x(2)  = xp(2,p_i) - len_phys(2)
        do aj = 2,nmom
            p_moments(aj) = p_moments(aj) + &
 &                      p_h(1)*p_h(2)*wp(1,p_i)*x(1)**alpha(1,aj)* &
 &                      x(2)**alpha(2,aj)
        enddo
        p_moments(1) = p_moments(1) + wp(1,p_i)*p_h(1)*p_h(2)
    enddo
    do j=1,ndata(2,1)-1
        do i=1,ndata(1,1)-1
               field_x(1) = min_phys(1) + h(1)*real(i-1,mk)
               field_x(2) = min_phys(2) + h(2)*real(j-1,mk)
               x = field_x
               if (field_x(1).ge.len_phys(1)) x(1)  = field_x(1) - len_phys(1)
               if (field_x(2).ge.len_phys(2)) x(2)  = field_x(2) - len_phys(2)
               do aj = 2,nmom
                   f_moments(aj) = f_moments(aj) + &
 &                      h(1)*h(2)*field_wp(1,i,j,1)*x(1)**alpha(1,aj)* &
 &                      x(2)**alpha(2,aj)
            enddo
            f_moments(1) = f_moments(1) + field_wp(1,i,j,1)*h(1)*h(2)
        enddo
    enddo

    do aj = 1,6
        if (abs(p_moments(aj) - f_moments(aj)) .GT. tol) then
            print *, 'failed at moment: ', aj
            print *, 'particle moments: ',p_moments
            print *, 'field moments: ',   f_moments
            stop 'ERROR: m2p interpolation: moments not conserved.'
        endif
    enddo


    !----------------
    ! teardown...
    !----------------
    print *, 'cleaning up'
8000 continue

    call ppm_finalize(info)

#ifdef __MPI
    call MPI_Finalize(info)
#endif
    
    deallocate(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

end program ppm_test_interp_m2p
