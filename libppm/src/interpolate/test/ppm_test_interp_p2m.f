      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_interp_p2m
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
program ppm_test_interp_p2m

    use ppm_module_typedef
    use ppm_module_mktopo
    use ppm_module_topo_get
    use ppm_module_interp_m2p
    use ppm_module_interp_p2m
    use ppm_module_init
    use ppm_module_finalize
    use ppm_module_map

    implicit none

#include "../../ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif

    integer, parameter              :: debug = 0
    integer, parameter              :: MK = ppm_kind_double
    integer                         :: ndim,nspec
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    real(MK)                        :: tol
    integer                         :: info
    integer                         :: topoid,meshid
    real(MK),dimension(:,:),pointer :: xp,wp
    real(MK),dimension(:  ),pointer :: min_phys,max_phys,h
    integer, dimension(:  ),pointer :: ghostsize
    integer                         :: i,j,p_i,ai,aj
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
    type(ppm_t_topo), pointer       :: topo
    real(MK)                        :: maxm3
    integer, parameter              :: nmom = 10
    integer, dimension(2,nmom)      :: alpha
    real(MK),dimension(nmom)        :: f_moments, p_moments
    !----------------
    ! setup
    !----------------
    tol = 10.0_mk*epsilon(1.0_mk)
    tolexp = int(log10(epsilon(1.0_mk)))
    ndim = 2
    nspec = 1
    data ((alpha(ai,aj), ai=1,2), aj=1,nmom) /0,0, 1,0, 0,1, 2,0, 0,2, &
   &                                          1,1, 3,0, 0,3, 2,1, 1,2/

    allocate(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),field_x(ndim),stat=info)

    do i=1,ndim
        min_phys(i) = 0.0_mk
        max_phys(i) = 1.0_mk
        ghostsize(i) = 2
    enddo
    bcdef(1:6) = ppm_param_bcdef_freespace

    nullify(xp)
    nullify(wp)

#ifdef __MPI
    call MPI_Init(info)
#endif

    call ppm_init(ndim,MK,tolexp,0,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = 400
    mp = 0

    allocate(xp(ndim,np),wp(nspec,np),stat=info)
    !call random_seed(put=(/17,42/))
    call random_number(xp)
    wp = 0.0_mk


    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_user_defined
    assig = ppm_param_assign_user_defined

    topoid = 0
    meshid = -1

    allocate(minsub(ndim,1),maxsub(ndim,1),sub2proc(1),cost(1),stat=info)

    cost(1) = 1.0_mk

    do i=1,ndim
        minsub(i,1)=min_phys(i)
        maxsub(i,1)=max_phys(i)
    enddo

    sub2proc(1) = 0

    allocate(nm(ndim),stat=info)
    do i=1,ndim
        nm(i) = 32
    enddo

    call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=1,user_sub2proc=sub2proc)

    allocate(field_wp(nspec,(1-ghostsize(1)):(ndata(1,1)+ghostsize(1)),         &
    &        (1-ghostsize(2)):(ndata(2,1)+ghostsize(2)),1),stat=info) ! 2d
    !allocate(field_wp(nspec,ndata(1,1),ndata(2,1),ndata(3,1),1),stat=info) ! 3d

    do i=1,ndim
        h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
    enddo

    print *, 'finished initialization'

    call ppm_map_part_global(topoid,xp,np,info) ! positions
    call ppm_map_part_push(wp,nspec,np,info)    ! strengths
    call ppm_map_part_send(np,mp,info)          ! send
    call ppm_map_part_pop(wp,nspec,np,mp,info)  ! strengths
    call ppm_map_part_pop(xp,ndim,np,mp,info)   ! positions
    if (info .NE. 0) then
       print *, 'Failed to do global mapping'
       goto 8000
    endif
    np = mp

    maxm3 = 0.0_mk
    do p_i=1,np
        !----------------
        ! p --> m
        !----------------
        wp(1,p_i) = 1.0_mk
        call ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,field_wp,info)

        !----------------
        ! test p --> m
        !----------------
        f_moments = 0.0_mk
        p_moments = 0.0_mk
	    do j = 1-ghostsize(2), ndata(2,1)+ghostsize(2)
	        do i = 1-ghostsize(1),ndata(1,1)+ghostsize(1)
	            field_x(1) = min_phys(1) + h(1)*real(i-1,mk)
	            field_x(2) = min_phys(2) + h(2)*real(j-1,mk)
	            do aj = 1,nmom
	               f_moments(aj) = f_moments(aj) + field_wp(1,i,j,1)* &
	&                            field_x(1)**alpha(1,aj)*field_x(2)**alpha(2,aj)
	            enddo
	        enddo
	    enddo
	    do aj = 1,nmom
	       p_moments(aj) = xp(1,p_i)**alpha(1,aj)*xp(2,p_i)**alpha(2,aj)
	    enddo
	    do aj = 1,6
            if (abs(f_moments(aj) - p_moments(aj)) .GT. tol) then
                print *, 'particle pos:',     xp(:,p_i)
                print *, 'failed at moment: ', aj
                print *, 'field moments: ',   f_moments
                print *, 'particle moments: ',p_moments
                stop 'ERROR: p2m interpolation: moments not conserved.'
            endif
	    enddo
        do aj = 7,10
            if (abs(f_moments(aj) - p_moments(aj)) .GT. maxm3) then
                maxm3 = abs(f_moments(aj) - p_moments(aj))
            endif
        enddo

	    wp(1,p_i) = 0.0_mk
	enddo

	print *, 'Maximum 3rd moment diff / h^3', maxm3/h**3
    !----------------
    ! m --> p
    !----------------


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

end program ppm_test_interp_p2m
