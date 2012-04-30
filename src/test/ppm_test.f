      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_CASE
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
program ppm_test_CASE

    use ppm_module_core
    use ppm_module_typedef
    use ppm_module_mktopo
    use ppm_module_init
    use ppm_module_finalize

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
    integer,  parameter             :: ngrid = 513
    integer,  parameter             :: npgrid = 1025
    real(MK),dimension(:,:),pointer :: xp,wp
    real(MK),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
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
    type(ppm_t_topo), pointer       :: topo
    real(MK)                        :: maxm3
    integer, parameter              :: nmom = 10
    integer, dimension(2,nmom)      :: alpha
    real(MK),dimension(nmom)        :: f_moments, p_moments
    real(MK),dimension(2)           :: sigma,mu
    !================
    ! setup
    !================
    tol = 10.0_mk*epsilon(1.0_mk)
    tolexp = int(log10(epsilon(1.0_mk)))
    ndim = 2
    nspec = 1

    allocate(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),p_h(ndim),field_x(ndim),stat=info)

    do i=1,ndim
        min_phys(i) = 0.0_mk
        max_phys(i) = 1.0_mk
        ghostsize(i) = 2
    enddo
    bcdef(1:6) = -1 ! BOUNDARY CONDITION

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

    !----------------
    ! make topology
    !----------------

    topoid = 0
    meshid = -1

    call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info)

    !----------------
    ! setup mesh data
    !----------------


    print *, 'finished initialization'
    !================
    ! body
    !================

    !---------------------------
    ! call routines to be tested
    !---------------------------

    !----------------
    ! test output
    !----------------




    !================
    ! teardown...
    !================
    print *, 'cleaning up'
8000 continue

    call ppm_finalize(info)

#ifdef __MPI
    call MPI_Finalize(info)
#endif
    
    deallocate(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

end program ppm_test_CASE
