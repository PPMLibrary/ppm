!-------------------------------------------------------------------------
!     Test Case   :                   ppm_test_inl
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
module old_inl
implicit none
contains

subroutine old_inl_vlist_2d(topo_id,topo,xp,rcp,cutoff,npart,mpart,&
        nvlist,vlist,info)


    USE ppm_module_typedef
    USE ppm_module_neighlist
    implicit none


    ! arguments
    integer, parameter              :: mk = ppm_kind_double
    integer,                         INTENT(IN   )   :: topo_id
    type(ppm_t_topo), pointer,       INTENT(IN   )   :: topo
    REAL(MK),dimension(:,:),pointer, INTENT(IN   )   :: xp
    REAL(MK),dimension(:),pointer,   INTENT(INOUT)   :: rcp
    REAL(MK),                        INTENT(IN   )   :: cutoff
    integer,                         INTENT(IN   )   :: npart
    integer,                         INTENT(IN   )   :: mpart
    integer, dimension(:),  pointer, INTENT(INOUT)   :: nvlist
    integer, dimension(:,:),pointer, INTENT(INOUT)   :: vlist
    integer,                         INTENT(  OUT)   :: info

    ! local variable
    integer                           :: ip,iq,isub,iinter
    integer                           :: ipart,jpart
    integer                           :: jbox,cbox
    integer                           :: n1,n2,i,j
    integer                           :: maxvlen
    integer                           :: ibegin,iend,jbegin,jend
    REAL(MK)                          :: dist2,cutoff2
    REAL(MK),dimension(2)          :: cellsize
    integer, dimension(:,:), pointer              :: ncells => null()
    type(ppm_type_ptr_to_clist), dimension(:), pointer  :: clist => null()
    integer,                     dimension(:,:),pointer :: ind => null()
    integer,                     dimension(:,:),pointer :: jnd => null()
    integer                                             :: nnd
    integer                                             :: ndim = 2


    !!-------------------------------------------------------------------------!
    !! initialize
    !!-------------------------------------------------------------------------!
    info = 0

    !!-------------------------------------------------------------------------!
    !! create cell lists with maximum cutoff
    !!-------------------------------------------------------------------------!
    cellsize = cutoff
    call ppm_neighlist_clist(topo_id,xp,mpart,cellsize,&
        .FALSE.,clist,ncells,info)

    !!-------------------------------------------------------------------------!
    !! create the index list of cell-cell interactons
    !!-------------------------------------------------------------------------!
    call ppm_neighlist_mkneighidx(.FALSE.,ind,jnd,nnd,info)

    !!-------------------------------------------------------------------------!
    !! run over cells and create verlet lists for the particles inside
    !!-------------------------------------------------------------------------!

    !reallocate nvlist only if the number of particles has changed
    if (associated(nvlist)) then
        if (size(nvlist).lt.npart) then
            DEALLOCATE(nvlist)
            ALLOCATE(nvlist(npart),stat=info)
        endif
    else
        ALLOCATE(nvlist(npart),stat=info)
    endif

    ! reset nvlist
    do ip = 1,npart
        nvlist(ip) = 0
    enddo

    !!-------------------------------------------------------------------------!
    !! determine size of verlet lists
    !!-------------------------------------------------------------------------!
    do isub=1,topo%nsublist
        n1  = ncells(1,isub)
        n2  = ncells(1,isub)*ncells(2,isub)
        ! loop over all cells
        do j=0,ncells(2,isub)-1
            do i=0,ncells(1,isub)-1

                ! index of the center box
                cbox = i + 1 + n1*j
                do iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
                    jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
                    !special treatment for ghost cells
                    if(jbox.le.0) cycle
                    if(jbox.gt.(product(ncells(:,isub))) ) cycle

                    !  get pointers to first and last particle 
                    ibegin = clist(isub)%lhbx(cbox)
                    iend   = clist(isub)%lhbx(cbox+1)-1
                    if (iend .lt. ibegin) cycle

                    ! within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    if (cbox .eq. jbox) then
                        do ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            if (ip .le. npart) then
                                do jpart=ibegin,iend

                                    if (jpart .ne. ipart) then
                                        ! translate to real particle index
                                        iq = clist(isub)%lpdx(jpart) 

                                        ! if not a ghost cell
                                        dist2=&
                                            sum((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 

                                        cutoff2 = min(rcp(iq),rcp(ip))**2
                                        if (dist2 .le. cutoff2) then
                                            ! add particle iq to 
                                            !list of particle ip
                                            nvlist(ip) = nvlist(ip) + 1
                                        endif
                                    endif
                                enddo
                            endif
                        enddo

                        !  for the other boxes check all particles
                    else
                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1
                        ! skip this iinter if empty
                        if (jend .lt. jbegin) cycle
                        ! loop over all particles inside this cell
                        do ipart=ibegin,iend 
                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 
                            if (ip .le. npart) then
                                ! check against all particles in the other cell
                                do jpart=jbegin,jend
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 
                                    dist2=&
                                        sum((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                    cutoff2 = min(rcp(iq),rcp(ip))**2
                                    if (dist2 .le. cutoff2) then
                                        !add particle 
                                        !iq to list of particle ip
                                        nvlist(ip) = nvlist(ip) + 1
                                    endif
                                enddo
                            endif
                        enddo

                    endif       ! cbox .eq. jbox
                enddo          ! iinter
            enddo             ! i
        enddo                ! j
    enddo                      ! isub
    !!-------------------------------------------------------------------------!
    !! allocate verlet list length
    !!-------------------------------------------------------------------------!
    maxvlen = maxval(nvlist)
    if (associated(vlist)) then
        if(size(vlist,1).lt.maxvlen .or. size(vlist,2).lt.npart) then
            DEALLOCATE(vlist)
            ALLOCATE(vlist(maxvlen,npart),stat=info)
        endif
    else
        ALLOCATE(vlist(maxvlen,npart),stat=info)
    endif

    ! reset nvlist
    do ip = 1,size(nvlist)
        nvlist(ip) = 0
    enddo

    !!-------------------------------------------------------------------------!
    !! fill verlet lists
    !!-------------------------------------------------------------------------!
    do isub=1,topo%nsublist
        n1  = ncells(1,isub)
        n2  = ncells(1,isub)*ncells(2,isub)
        ! loop over all real cells (the -2 at the end does this)
        do j=1,ncells(2,isub)-2
            do i=1,ncells(1,isub)-2
                ! index of the center box
                cbox = i + 1 + n1*j
                do iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
                    jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
                    !  get pointers to first and last particle 
                    ibegin = clist(isub)%lhbx(cbox)
                    iend   = clist(isub)%lhbx(cbox+1)-1
                    if (iend .lt. ibegin) cycle

                    ! within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    if (cbox .eq. jbox) then
                        do ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            if (ip .le. npart) then
                                do jpart=ibegin,iend

                                    if (jpart .ne. ipart) then
                                        ! translate to real particle index
                                        iq = clist(isub)%lpdx(jpart) 

                                        dist2=&
                                            sum((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                        cutoff2 = min(rcp(iq),rcp(ip))**2
                                        if (dist2 .le. cutoff2) then
                                            ! add particle iq to 
                                            ! list of particle ip
                                            nvlist(ip) = nvlist(ip) + 1
                                            vlist(nvlist(ip),ip) = iq
                                        endif
                                    endif
                                enddo
                            endif
                        enddo

                        !  for the other boxes check all particles
                    else

                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        if (jend .lt. jbegin) cycle

                        do ipart=ibegin,iend ! loop over all particles inside this cell

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            if (ip .le. npart) then
                                ! check against all particles in the other cell
                                do jpart=jbegin,jend

                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 

                                    dist2=&
                                        sum((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                    cutoff2 = min(rcp(iq),rcp(ip))**2
                                    if (dist2 .le. cutoff2) then
                                        ! add particle iq to list 
                                        ! of particle ip
                                        nvlist(ip) = nvlist(ip) + 1
                                        vlist(nvlist(ip),ip) = iq
                                    endif

                                enddo
                            endif
                        enddo
                    endif       ! cbox .eq. jbox
                enddo          ! iinter

            enddo             ! i
        enddo                ! j
    enddo                      ! isub
    9999 continue ! jump here upon error

end subroutine old_inl_vlist_2d

end module old_inl

program ppm_test_inl
!-------------------------------------------------------------------------
!     Test Case   :   Inhomogeneous neighbour lists - 2D
!                           Particles placed on a grid
!                           cutoff radii are uniformly random
!-------------------------------------------------------------------------

USE ppm_module_typedef
USE ppm_module_mktopo
USE ppm_module_topo_get
USE ppm_module_init
USE ppm_module_finalize
USE ppm_module_core
USE ppm_module_inl_vlist
use old_inl

implicit none
#include "../../ppm_define.h"
#ifdef __MPI
INCLUDE 'mpif.h'
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = ppm_kind_double
REAL(MK),parameter              :: pi = 3.1415926535897931_mk
REAL(MK),parameter              :: skin = 0.0_MK
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
REAL(MK)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank
integer                         :: topoid
integer                         :: np,npgrid = 123
integer                         :: mp
integer                         :: newnp
integer, dimension(:),  pointer :: p_i
REAL(MK),dimension(:,:),pointer :: xp
REAL(MK),dimension(:  ),pointer :: rcp
REAL(MK),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
REAL(MK),dimension(:  ),pointer :: len_phys
REAL(MK),dimension(:  ),pointer :: ghostlayer
REAL(MK),dimension(:,:),pointer :: minsub,maxsub
integer, dimension(:  ),pointer :: sub2proc
integer                         :: ii,i,j,sum1,sum2
integer, dimension(6)           :: bcdef
REAL(MK),dimension(:  ),pointer :: cost
type(ppm_t_topo), pointer       :: topo
integer                         :: seedsize
integer,  dimension(:),ALLOCATABLE :: seed
REAL(MK), dimension(:),ALLOCATABLE :: randnb
integer,dimension(:,:),pointer   :: vlist,vlist2
integer,dimension(:),  pointer   :: nvlist,nvlist2
integer                          :: isymm = 0
logical                          :: lsymm = .FALSE.,ok

!----------------
! setup
!----------------
tol = 10.0_MK*EPSILON(1.0_MK)
tolexp = int(log10(EPSILON(1.0_MK)))
min_rcp = 0.01_mk
max_rcp = 0.1_mk
np=npgrid*npgrid

ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &         ghostlayer(2*ndim),&
    &         h(ndim),p_i(ndim),p_h(ndim),stat=info)

min_phys = 0.0_MK
max_phys = 1.0_MK
len_phys = max_phys-min_phys
ghostlayer(1:2*ndim) = max_rcp
bcdef(1:6) = ppm_param_bcdef_periodic

NULLIFY(xp,rcp)

#ifdef __MPI
comm = MPI_COMM_WORLD
call mpi_init(info)
CALL MPI_Comm_rank(comm,rank,info)
#else
rank = 0
#endif
call ppm_init(ndim,mk,tolexp,0,debug,info,99)

call random_seed(size=seedsize)
ALLOCATE(seed(seedsize))
ALLOCATE(randnb((ndim+1)*np),stat=info)
do i=1,seedsize
    seed(i)=10+i*i*(rank+1)
enddo
call random_seed(put=seed)
call random_number(randnb)

!----------------
! create particles on a grid
!----------------

ALLOCATE(xp(ndim,np),rcp(np),stat=info)
xp = 0.0_MK
rcp = 0.0_MK

p_h = len_phys / REAL(npgrid,mk)

do j=1,npgrid
    do i=1,npgrid
        p_i = i + (j-1)*npgrid
        xp(1,p_i) = min_phys(1)+REAL(i-1,mk)*p_h(1)
        xp(2,p_i) = min_phys(2)+REAL(j-1,mk)*p_h(2)
        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
    enddo
enddo

!----------------
! make topology
!----------------
assig  = ppm_param_assign_internal
topoid = 0

decomp = ppm_param_decomp_user_defined
ALLOCATE(minsub(ndim,4),maxsub(ndim,4),sub2proc(4),cost(4),stat=info)

cost = 1.0_MK; sub2proc = 0
minsub(1:ndim,1)=0.0_MK; maxsub(1:ndim,1)=0.5_mk

minsub(1,2)=0.5_mk; minsub(2,2)=0.0_MK
maxsub(1,2)=1.0_MK; maxsub(2,2)=0.5_mk

minsub(1,3)=0.0_MK; minsub(2,3)=0.5_mk
maxsub(1,3)=0.5_mk; maxsub(2,3)=1.0_MK

minsub(1:ndim,4)=0.5_mk; maxsub(1:ndim,4)=1.0_MK

topoid = 0
call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               max_rcp,cost,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=4,user_sub2proc=sub2proc)
call ppm_topo_get(topoid,topo,info)

!do multiple runs; each time, particle positions are perturbed
!by a different amount (starting from 0)
do ii=1,5
    call random_number(randnb)
    do i=1,np
        xp(1:ndim,i)=xp(1:ndim,i)+&
            0.05_MK*REAL(ii-1,mk)*randnb(ndim*i-(ndim-1):ndim*i)
    enddo
    !apply periodic boundary conditions
    do i=1,np
        do j=1,ndim
            xp(j,i)=MOD(xp(j,i),1.0_MK)
        enddo
    enddo

    !map particles onto the topology
    call ppm_map_part_global(topoid,xp,np,info)
    call ppm_map_part_push(rcp,np,info)
    call ppm_map_part_send(np,newnp,info)
    call ppm_map_part_pop(rcp,np,newnp,info)
    call ppm_map_part_pop(xp,ndim,np,newnp,info)
    np=newnp

    call ppm_topo_check(topoid,xp,np,ok,info)
    if (.not. ok) write(*,*) '[',rank,'] topo_check failed'

    call ppm_map_part_ghost_get(topoid,xp,ndim,np,isymm,max_rcp,info)
    call ppm_map_part_push(rcp,np,info)
    call ppm_map_part_send(np,mp,info)
    call ppm_map_part_pop(rcp,np,mp,info)
    call ppm_map_part_pop(xp,ndim,np,mp,info)

    call ppm_dbg_print(topoid,max_rcp,1,1,info,xp,np,mp)


    call ppm_inl_vlist(topoid,xp,np,mp,rcp,skin, &
        & lsymm,ghostlayer,info,vlist,nvlist)

    call old_inl_vlist_2d(topoid,topo,xp,rcp,max_rcp,np,mp,&
        nvlist2,vlist2,info)

    !compare neighbour lists obtained by the 2 different routines
    do i=1,np
        if (nvlist(i).NE.nvlist2(i)) then
            write(*,'(A)') '!! --------------- !!'
            write(*,'(A)') '!! FAILED'
            print *, '!!    nvlist are not equal for ip = ',i
            write(*,'(A)') '!! --------------- !!'
            goto 8000
        endif
        sum1=0;sum2=0;
        do j=1,nvlist(i)
            sum1 = sum1 + vlist(j,i)**2
            sum2 = sum2 + vlist2(j,i)**2
        enddo
        if (sum1 .NE. sum2) then
            write(*,'(A)') '!! --------------- !!'
            write(*,'(A)') '!! FAILED'
            print *, '!!    vlist are not equal for ip = ',i
            write(*,'(A)') '!! --------------- !!'
            goto 8000
        endif
    enddo

enddo

if (rank.eq.0) then
    write(*,'(A)') '!! --------------- !!'
    write(*,'(A)') '!! SUCCESS'
    write(*,'(2(A,I6),A)') '!!     Completed test for INL with ',np,&
    ' real particles and ',mp-np,' ghosts'
    write(*,'(2(A,E10.4))') '!!     Smallest cutoff was ', minval(rcp(1:np)),&
    ' and largest ',maxval(rcp(1:np))
    write(*,'(2(A,I6))') '!!     Smallest nb of neigh was ',&
    minval(nvlist),' and largest ',maxval(nvlist)
    write(*,'(A)') '!! --------------- !!'
endif

!----------------
! teardown...
!----------------
8000 continue

call ppm_finalize(info)
#ifdef __MPI
call MPI_finalize(info)
#endif

DEALLOCATE(xp,rcp,min_phys,max_phys,len_phys,minsub,maxsub,cost,sub2proc)

if (rank.eq.0) print *, 'done.'

end program ppm_test_inl
