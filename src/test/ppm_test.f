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

    USE ppm_module_core
    USE ppm_module_typedef
    USE ppm_module_mktopo
    USE ppm_module_init
    USE ppm_module_finalize
    USE ppm_module_mpi
    implicit none

#include "../../ppm_define.h"

    integer, parameter              :: debug = 0
    integer, parameter              :: MK = ppm_kind_double
    REAL(MK),parameter              :: pi = 3.1415926535897931_mk
    integer                         :: ndim,nspec
    integer                         :: decomp
    integer                         :: assig
    integer                         :: tolexp
    REAL(MK)                        :: tol
    integer                         :: info
    integer                         :: topoid,meshid
    integer,  parameter             :: ngrid = 513
    integer,  parameter             :: npgrid = 1025
    REAL(MK),dimension(:,:),pointer :: xp,wp
    REAL(MK),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
    integer, dimension(:  ),pointer :: ghostsize
    integer                         :: i,j,ai,aj,p_i
    integer, dimension(6)           :: bcdef
    REAL(MK),dimension(:  ),pointer :: cost
    integer, dimension(:,:),pointer :: istart,ndata
    integer, dimension(:  ),pointer :: nm
    REAL(MK),dimension(:,:),pointer :: minsub,maxsub
    integer, dimension(:  ),pointer :: sub2proc
    integer                         :: np,mp
    integer, parameter              :: kernel = ppm_param_rmsh_kernel_mp4
    REAL(MK),dimension(:,:,:,:  ), pointer :: field_wp ! 2d  field_up(ldn,i,j,isub)
    !REAL(MK),dimension(:,:,:,:,:), pointer :: field_up ! 3d  field_up(ldn,i,j,k,isub)
    REAL(MK),dimension(:  ),pointer :: field_x
    type(ppm_t_topo), pointer       :: topo
    REAL(MK)                        :: maxm3
    integer, parameter              :: nmom = 10
    integer, dimension(2,nmom)      :: alpha
    REAL(MK),dimension(nmom)        :: f_moments, p_moments
    REAL(MK),dimension(2)           :: sigma,mu
    !================
    ! setup
    !================
    tol = 10.0_MK*EPSILON(1.0_MK)
    tolexp = int(log10(EPSILON(1.0_MK)))
    ndim = 2
    nspec = 1

    ALLOCATE(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),p_h(ndim),field_x(ndim),stat=info)

    do i=1,ndim
        min_phys(i) = 0.0_MK
        max_phys(i) = 1.0_MK
        ghostsize(i) = 2
    enddo
    bcdef(1:6) = -1 ! BOUNDARY CONDITION

    NULLIFY(xp)
    NULLIFY(wp)

#ifdef __MPI
    call MPI_Init(info)
#endif

    call ppm_init(ndim,MK,tolexp,0,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = npgrid**2
    mp = 0

    ALLOCATE(xp(ndim,np),wp(nspec,np),stat=info)
    xp = 0.0_MK
    wp = 0.0_MK


    ALLOCATE(nm(ndim),stat=info)
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

    DEALLOCATE(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

end program ppm_test_CASE
