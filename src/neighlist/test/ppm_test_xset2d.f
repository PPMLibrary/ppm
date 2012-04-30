!-------------------------------------------------------------------------
!     test case   :                   ppm_test_xset
!-------------------------------------------------------------------------
! copyright (c) 2012 cse lab (eth zurich), mosaic group (eth zurich), 
!                    center for fluid dynamics (dtu)
!
!
! this file is part of the parallel particle mesh library (ppm).
!
! ppm is free software: you can redistribute it and/or modify
! it under the terms of the gnu lesser general public license 
! as published by the free software foundation, either 
! version 3 of the license, or (at your option) any later 
! version.
!
! ppm is distributed in the hope that it will be useful,
! but without any warranty; without even the implied warranty of
! merchantability or fitness for a particular purpose. see the
! gnu general public license for more details.
!
! you should have received a copy of the gnu general public license
! and the gnu lesser general public license along with ppm. if not,
! see <http://www.gnu.org/licenses/>.
!
! parallel particle mesh library (ppm)
! eth zurich
! ch-8092 zurich, switzerland
!-------------------------------------------------------------------------
module old_xset
implicit none
contains

subroutine old_xset_vlist_2d(topo_id,topo,xp_red,nred, &
        xp_blue,rcp_blue,nblue,mblue,cutoff,nvlist_cross,vlist_cross,info)

#undef __3D

    use ppm_module_typedef
    use ppm_module_neighlist
    implicit none


    ! arguments
    integer, parameter              :: mk = ppm_kind_double
    integer,                         intent(in   ) :: topo_id
    type(ppm_t_topo), pointer,       intent(in   ) :: topo
    real(mk),dimension(:,:),pointer, intent(in   ) :: xp_red
    integer,                         intent(in   ) :: nred
    real(mk),dimension(:,:),         intent(in   ) :: xp_blue
    real(mk),dimension(:),           intent(in   ) :: rcp_blue
    integer,                         intent(in   ) :: nblue
    integer,                         intent(in   ) :: mblue
    real(mk),                        intent(in   ) :: cutoff
    integer, dimension(:),  pointer, intent(inout) :: nvlist_cross
    integer, dimension(:,:),pointer, intent(inout) :: vlist_cross
    integer,                         intent(  out)   :: info

    ! local variable
    integer                                             :: ip,iq,isub,iinter
    integer                                             :: ipart,jpart
    integer                                             :: ibox,jbox,cbox
    integer                                             :: n1,n2,n3,i,j,k
    integer                                             :: maxvlen
    integer                                             :: ibegin,iend,jbegin,jend
    real(mk)                                            :: dist2,cutoff2
    real(mk),dimension(2)                               :: cellsize
    integer, dimension(:,:), pointer                    :: ncells => null()
    type(ppm_type_ptr_to_clist), dimension(:), pointer  :: clist => null()
    integer,                     dimension(:,:),pointer :: ind => null()
    integer,                     dimension(:,:),pointer :: jnd => null()
    integer                                             :: nnd
    integer                                             :: ndim
    real(mk), dimension(:,:),pointer                    :: all_xp=>NULL()

        !---------------------------------------------------------------------!
        ! initialize
        !---------------------------------------------------------------------!
        ndim = 2
        info = 0
        !---------------------------------------------------------------------!
        ! first fill an array with the new and old particles
        !---------------------------------------------------------------------!
        allocate(all_xp(ndim,mblue+nred),stat=info)
        if(associated(nvlist_cross)) deallocate(nvlist_cross)
        allocate(nvlist_cross(nred),stat=info)
        if(associated(vlist_cross)) deallocate(vlist_cross)

        all_xp(1:ndim,1:nred) = xp_red(1:ndim,1:nred)
        all_xp(1:ndim,nred+1:nred+mblue) = xp_blue(1:ndim,1:mblue)

        !!--------------------------------------------------------------------
        !! create cell lists with maximum cutoff
        !!--------------------------------------------------------------------
        cellsize = cutoff

        call ppm_neighlist_clist(topo_id,all_xp,mblue+nred,cellsize,&
            .false.,clist,ncells,info)

        !!--------------------------------------------------------------------
        !! create the index list of cell-cell interactons
        !!--------------------------------------------------------------------
        call ppm_neighlist_mkneighidx(.false.,ind,jnd,nnd,info)

        !!--------------------------------------------------------------------
        !! run over cells and create verlet lists for the particles inside
        !!--------------------------------------------------------------------

        ! initialise nvlist_cross to zero
        do ip = 1,nred
            nvlist_cross(ip) = 0
        enddo

        !!---------------------------------------------------------------------!
        !!  determine size of verlet lists
        !!---------------------------------------------------------------------!
        do isub=1,topo%nsublist
            n1  = ncells(1,isub)
            n2  = ncells(1,isub) * ncells(2,isub)
#ifdef __3D
            n3  = ncells(3,isub)
#endif
            ! loop over all cells
#ifdef __3D
            do k=0,n3-1
#endif
                do j=0,ncells(2,isub)-1
                    do i=0,ncells(1,isub)-1

                        ! index of the center box
#ifdef __3D
                        cbox = i + 1 + n1*j + n2*k
#else
                        cbox = i + 1 + n1*j
#endif
                    do iinter=1,nnd ! loop over all box-box interactions

                        ! determine box indices for this interaction
#ifdef __3D
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter) + &
                            n2*jnd(3,iinter))
#else
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
#endif


                        !special treatment for ghost cells
                        if(jbox.le.0) cycle
                        if(jbox.gt.(product(ncells(1:ndim,isub))) ) cycle

                        !  get pointers to first and last particle 
                        ibegin = clist(isub)%lhbx(cbox)
                        iend   = clist(isub)%lhbx(cbox+1)-1
                        if (iend .lt. ibegin) cycle

                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        if (jend .lt. jbegin) cycle

                        ! loop over all particles inside this cell
                        ipart_loop: do ipart=ibegin,iend 
                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 
                            ! only consider the new (real) particles as ip
                            if (ip .gt. nred) cycle
                            ! check against all particles in the other cell
                            jpart_loop: do jpart=jbegin,jend
                                if (jpart .ne. ipart) then
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 

                                    if (iq .le. nred) then
                                        ! iq belongs to the new set
                                    else
                                        ! iq belongs to the old set 
                                        dist2= sum((all_xp(1:ndim,ip)- &
                                            all_xp(1:ndim,iq))**2)
                                            cutoff2 = rcp_blue(iq-nred)**2

                                        if (dist2 .le. cutoff2) then
                                            ! populate nvlist_cross
                                            nvlist_cross(ip) = &
                                                nvlist_cross(ip) + 1
                                        endif
                                    endif
                                endif

                            enddo jpart_loop
                        enddo ipart_loop

                    enddo          ! iinter

                enddo             ! i
            enddo                ! j
#ifdef __3D
        enddo   ! k
#endif
        enddo                      ! isub

        !!--------------------------------------------------------------------!
        !! allocate verlet list length
        !!--------------------------------------------------------------------!
        maxvlen = maxval(nvlist_cross)
        allocate(vlist_cross(maxvlen,nred),stat=info)

        ! initialise nvlist_cross to zero
        do ip = 1,nred
            nvlist_cross(ip) = 0
        enddo

        !!--------------------------------------------------------------------!
        !! fill verlet lists
        !!--------------------------------------------------------------------!
        do isub=1,topo%nsublist
            n1  = ncells(1,isub)
            n2  = ncells(1,isub)*ncells(2,isub)
#ifdef __3D
            n3  = ncells(3,isub)
#endif
        ! loop over all real cells (the -2 at the end does this)
#ifdef __3D
        do k=1,n3-2
#endif
            do j=1,ncells(2,isub)-2
                do i=1,ncells(1,isub)-2
                    ! index of the center box
#ifdef __3D
                    cbox = i + 1 + n1*j + n2*k
#else
                    cbox = i + 1 + n1*j
#endif

                    do iinter=1,nnd ! loop over all box-box interactions

                        ! determine box indices for this interaction
#ifdef __3D
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter) + &
                            n2*jnd(3,iinter))
#else
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
#endif

                        !  get pointers to first and last particle 
                        ibegin = clist(isub)%lhbx(cbox)
                        iend   = clist(isub)%lhbx(cbox+1)-1
                        if (iend .lt. ibegin) cycle

                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        if (jend .lt. jbegin) cycle

                        ! loop over all particles inside this cell
                        ipart_loop2: do ipart=ibegin,iend 

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 
                            ! only consider the new (real) particles as ip
                            if (ip .gt. nred) cycle

                            ! check against all particles in the other cell
                            do jpart=jbegin,jend
                                if (jpart .ne. ipart) then
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 

                                    if (iq .le. nred) then
                                        ! iq belongs to the new set
                                    else
                                        ! iq belongs to the old set 
                                        dist2= sum((all_xp(1:ndim,ip)- &
                                            all_xp(1:ndim,iq))**2)
                                            cutoff2 = rcp_blue(iq-nred)**2

                                        if (dist2 .le. cutoff2) then
                                            !populate nvlist_cross
                                            nvlist_cross(ip) = &
                                                nvlist_cross(ip) + 1
                                            !vlist contains index of old particles
                                            vlist_cross(nvlist_cross(ip),ip) = &
                                                iq - nred
                                        endif
                                    endif
                                endif
                            enddo
                        enddo ipart_loop2
                    enddo          ! iinter
                enddo             ! i
            enddo                ! j
#ifdef __3D
        enddo                ! k
#endif
        enddo                      ! isub


        9999 continue ! jump here upon error
        deallocate(all_xp)

end subroutine old_xset_vlist_2d

end module old_xset

program ppm_test_xset
!-------------------------------------------------------------------------
!     test case   :   inhomogeneous neighbour lists - 2d
!                           particles placed uniformly randomly in a box
!                           cutoff radii are uniformly random
!-------------------------------------------------------------------------

use ppm_module_typedef
use ppm_module_mktopo
use ppm_module_topo_get
use ppm_module_init
use ppm_module_finalize
use ppm_module_core
use ppm_module_inl_xset_vlist
use old_xset

implicit none
#include "../../ppm_define.h"
#ifdef __MPI
include 'mpif.h'
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = ppm_kind_double
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank
integer                         :: topoid
integer                         :: nred = 10000
integer                         :: mred
integer                         :: nblue = 3000
integer                         :: mblue
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp_red=>NULL()
real(mk),dimension(:,:),pointer :: xp_blue=>NULL()
real(mk),dimension(:  ),pointer :: rcp_blue=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys,h,p_h=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
real(mk),dimension(:  ),pointer :: ghostlayer=>NULL()
integer, dimension(:  ),pointer :: ghostsize=>NULL()
integer                         :: i,j,sum1,sum2
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
integer, dimension(:  ),pointer :: nm=>NULL()
type(ppm_t_topo), pointer       :: topo
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer,dimension(:,:),pointer   :: vlist=>NULL(),vlist2=>NULL()
integer,dimension(:),  pointer   :: nvlist=>NULL(),nvlist2=>NULL()
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok

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
allocate(randnb((ndim)*nred),stat=info)
do i=1,seedsize
    seed(i)=10+i*i*(rank+1)
enddo
call random_seed(put=seed)
call random_number(randnb)

!----------------
! create particles
!----------------

allocate(xp_red(ndim,nred),stat=info)
allocate(xp_blue(ndim,nblue),rcp_blue(nblue),stat=info)
xp_red = 0.0_mk
xp_blue = 0.0_mk
rcp_blue = 0.0_mk

do i=1,nred
    do j=1,2
        xp_red(j,i) = min_phys(j)+&
            len_phys(j)*randnb(2*i-j+1)
    enddo
enddo

deallocate(randnb)
allocate(randnb((1+ndim)*nblue),stat=info)
call random_number(randnb)
do i=1,nblue
    do j=1,ndim
        xp_blue(j,i) = min_phys(j)+&
            len_phys(j)*randnb((ndim+1)*i-(ndim-j))
    enddo
    rcp_blue(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
enddo

!----------------
! make topology
!----------------
!decomp = ppm_param_decomp_cuboid
decomp = ppm_param_decomp_xpencil
assig  = ppm_param_assign_internal

topoid = 0

call ppm_mktopo(topoid,xp_red,nred,decomp,assig,min_phys,max_phys,bcdef, &
    &               max_rcp,cost,info)

call ppm_map_part_global(topoid,xp_red,nred,info)
call ppm_map_part_send(nred,newnp,info)
call ppm_map_part_pop(xp_red,ndim,nred,newnp,info)
nred=newnp

call ppm_map_part_global(topoid,xp_blue,nblue,info)
call ppm_map_part_push(rcp_blue,nblue,info)
call ppm_map_part_send(nblue,newnp,info)
call ppm_map_part_pop(rcp_blue,nblue,newnp,info)
call ppm_map_part_pop(xp_blue,ndim,nblue,newnp,info)
nblue=newnp

call ppm_topo_check(topoid,xp_red,nred,ok,info)
if (.not. ok) write(*,*) '[',rank,'] topo_check failed for red'
call ppm_topo_check(topoid,xp_blue,nblue,ok,info)
if (.not. ok) write(*,*) '[',rank,'] topo_check failed for blue'

call ppm_map_part_ghost_get(topoid,xp_red,ndim,nred,isymm,max_rcp,info)
call ppm_map_part_send(nred,mred,info)
call ppm_map_part_pop(xp_red,ndim,nred,mred,info)

call ppm_map_part_ghost_get(topoid,xp_blue,ndim,nblue,isymm,max_rcp,info)
call ppm_map_part_push(rcp_blue,nblue,info)
call ppm_map_part_send(nblue,mblue,info)
call ppm_map_part_pop(rcp_blue,nblue,mblue,info)
call ppm_map_part_pop(xp_blue,ndim,nblue,mblue,info)

call ppm_dbg_print(topoid,max_rcp,2,1,info,xp_red,nred,mred)

call ppm_dbg_print(topoid,max_rcp,2,1,info,xp_blue,nblue,mblue)

call ppm_inl_xset_vlist(topoid,xp_red,nred,mred,&
 &                   xp_blue,nblue,mblue,rcp_blue,skin,    &
 &                   ghostlayer,info,vlist,nvlist)

call ppm_topo_get(topoid,topo,info)
    
call old_xset_vlist_2d(topoid,topo,xp_red,nred, &
        xp_blue,rcp_blue,nblue,mblue,max_rcp,nvlist2,vlist2,info)

!compare neighbour lists obtained by the 2 different routines
do i=1,nred
    if (nvlist(i).ne.nvlist2(i)) then
        write(*,'(a)') '!! --------------- !!'
        write(*,'(a)') '!! failed'
        print *, '!!    nvlist are not equal for ip = ',i
        print * ,'!! nvlist:  ',nvlist(i)
        print * ,'!! nvlist2: ',nvlist2(i)
        write(*,'(a)') '!! --------------- !!'
        goto 8000
    endif
    sum1=0;sum2=0;
    do j=1,nvlist(i)
        sum1 = sum1 + vlist(j,i)**2
        sum2 = sum2 + vlist2(j,i)**2
    enddo
    if (sum1 .ne. sum2) then
        write(*,'(a)') '!! --------------- !!'
        write(*,'(a)') '!! failed'
        print *, '!!    vlist are not equal for ip = ',i
        write(*,'(a)') '!! --------------- !!'
        goto 8000
    endif
enddo


if (rank.eq.0) then
    write(*,'(a)') '!! --------------- !!'
    write(*,'(a)') '!! success'
    write(*,'(2(a,i6),a)') '!!     completed test for xset with ',nred,&
    ' red particles and ',nblue,' blue particles'
    write(*,'(2(a,e10.4))') '!!     smallest cutoff was ', &
        minval(rcp_blue(1:nblue)),&
        ' and largest ',maxval(rcp_blue(1:nblue))
    write(*,'(2(a,i6))') '!!     smallest nb of neigh was ',&
        minval(nvlist),' and largest ',maxval(nvlist)
    write(*,'(a)') '!! --------------- !!'
endif

!----------------
! teardown...
!----------------
8000 continue

call ppm_finalize(info)
#ifdef __MPI
call mpi_finalize(info)
#endif

!deallocate(xp_red,min_phys,max_phys,len_phys,ghostsize,nm)
!deallocate(xp_blue,rcp_blue,nvlist,nvlist2,vlist,vlist2)

if (rank.eq.0) print *, 'done.'


end program ppm_test_xset
