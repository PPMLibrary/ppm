!-------------------------------------------------------------------------
!     Test Case   :                   ppm_test_map_part_ghost
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

program ppm_test_map_part_ghost
!-------------------------------------------------------------------------
!     Test Case   :  ghost mappings, pushs, pops and multi-push/pops 
!-------------------------------------------------------------------------

use ppm_module_typedef
use ppm_module_mktopo
use ppm_module_topo_get
use ppm_module_init
use ppm_module_finalize
use ppm_module_core

implicit none
#include "../../ppm_define.h"
#ifdef __MPI
INCLUDE 'mpif.h'
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = ppm_kind_double
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer,parameter               :: pdim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank
integer                         :: topoid
integer                         :: np = 100000
integer                         :: mp
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp
real(mk),dimension(:  ),pointer :: rcp
real(mk),dimension(:,:),pointer :: wp
real(mk),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: ghostlayer
integer, dimension(:  ),pointer :: ghostsize
integer                         :: i,j,k,sum1,sum2
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
integer, dimension(:  ),pointer :: nm
type(ppm_t_topo), pointer       :: topo
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3

!----------------
! setup
!----------------
tol = 10.0_mk*epsilon(1.0_mk)
tolexp = int(log10(epsilon(1.0_mk)))
min_rcp = 0.01_mk
max_rcp = 0.1_mk

allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &         ghostsize(ndim),ghostlayer(2*ndim),&
    &         nm(ndim),h(ndim),p_h(ndim),stat=info)

min_phys(1:ndim) = 0.0_mk
max_phys(1:ndim) = 1.0_mk
len_phys(1:ndim) = max_phys-min_phys
ghostsize(1:ndim) = 2
ghostlayer(1:2*ndim) = max_rcp
bcdef(1:6) = ppm_param_bcdef_periodic

nullify(xp,rcp,wp)

#ifdef __MPI
comm = mpi_comm_world
call mpi_init(info)
call mpi_comm_rank(comm,rank,info)
#else
rank = 0
#endif
call ppm_init(ndim,mk,tolexp,0,debug,info,99)

call random_seed(size=seedsize)
allocate(seed(seedsize))
allocate(randnb((1+ndim)*np),stat=info)
do i=1,seedsize
    seed(i)=10+i*i*(rank+1)
enddo
call random_seed(put=seed)
call random_number(randnb)

!----------------
! create particles
!----------------

allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)
xp = 0.0_mk
rcp = 0.0_mk

!p_h = len_phys / real(npgrid,mk)
!do j=1,npgrid
!    do i=1,npgrid
!        p_i = i + (j-1)*npgrid
!        xp(1,p_i) = min_phys(1)+real(i-1,mk)*p_h(1)
!        xp(2,p_i) = min_phys(2)+real(j-1,mk)*p_h(2)
!        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
!        do k=1,pdim
!            wp(k,i) = rcp(i)*REAL(k,MK)
!        enddo
!    enddo
!enddo
do i=1,np
    do j=1,ndim
        xp(j,i) = min_phys(j)+&
            len_phys(j)*randnb((ndim+1)*i-(ndim-j))
    enddo
    rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
    do j=1,pdim
        wp(j,i) = rcp(i)*REAL(j,MK)
    enddo
enddo

!----------------
! make topology
!----------------
decomp = ppm_param_decomp_cuboid
!decomp = ppm_param_decomp_xpencil
assig  = ppm_param_assign_internal

topoid = 0

call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               max_rcp,cost,info)

print *,rank,'global map:'
call ppm_map_part_global(topoid,xp,np,info)
call ppm_map_part_push(rcp,np,info)
call ppm_map_part_push(wp,pdim,np,info)
call ppm_map_part_send(np,newnp,info)
call ppm_map_part_pop(wp,pdim,np,newnp,info)
call ppm_map_part_pop(rcp,np,newnp,info)
call ppm_map_part_pop(xp,ndim,np,newnp,info)
np=newnp
print *,rank,'done'

call ppm_topo_check(topoid,xp,np,ok,info)
if (.not. ok) write(*,*) '[',rank,'] topo_check failed'

call ppm_map_part_ghost_get(topoid,xp,ndim,np,isymm,max_rcp,info)
call ppm_map_part_push(rcp,np,info)
call ppm_map_part_push(wp,pdim,np,info)
call ppm_map_part_send(np,mp,info)
call ppm_map_part_pop(wp,pdim,np,mp,info)
call ppm_map_part_pop(rcp,np,mp,info)
call ppm_map_part_pop(xp,ndim,np,mp,info)
print *,rank,'np: ',np,' ghosts: ',mp -np
print *,rank,'done'
print *,rank,'is ghost_get set?',ppm_map_type_isactive(ppm_param_map_ghost_get)
call ppm_dbg_print(topoid,max_rcp,1,1,info,xp,np,mp)




!compare neighbour lists obtained by the 2 different routines


!if (rank.eq.0) then
!    write(*,'(A)') '!! --------------- !!'
!    write(*,'(A)') '!! SUCCESS'
!    write(*,'(2(A,I6),A)') '!!     Completed test for INL with ',np,&
!        ' real particles and ',mp-np,' ghosts'
!    write(*,'(2(A,E10.4))') '!!     Smallest cutoff was ', minval(rcp(1:np)),&
!        ' and largest ',maxval(rcp(1:np))
!    write(*,'(2(A,I6))') '!!     Smallest nb of neigh was ',&
!        minval(nvlist),' and largest ',maxval(nvlist)
!    write(*,'(A)') '!! --------------- !!'
!endif

!----------------
! teardown...
!----------------
8000 continue

call ppm_finalize(info)
#ifdef __MPI
call MPI_finalize(info)
#endif

deallocate(xp,rcp,wp,min_phys,max_phys,len_phys,ghostsize,nm)

if (rank.eq.0) print *, 'done.'


end program ppm_test_map_part_ghost
