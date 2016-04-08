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
MODULE old_xset
IMPLICIT NONE
contains

SUBROUTINE old_xset_vlist_2d(topo_id,topo,xp_red,nred, &
        xp_blue,rcp_blue,nblue,mblue,cutoff,nvlist_cross,vlist_cross,info)

#undef __3D

    USE ppm_module_typedef
    USE ppm_module_neighlist
    IMPLICIT NONE


    ! arguments
    INTEGER, PARAMETER              :: mk = ppm_kind_double
    INTEGER,                         INTENT(IN   ) :: topo_id
    TYPE(ppm_t_topo), POINTER,       INTENT(IN   ) :: topo
    REAL(MK),DIMENSION(:,:),POINTER, INTENT(IN   ) :: xp_red
    INTEGER,                         INTENT(IN   ) :: nred
    REAL(MK),DIMENSION(:,:),         INTENT(IN   ) :: xp_blue
    REAL(MK),DIMENSION(:),           INTENT(IN   ) :: rcp_blue
    INTEGER,                         INTENT(IN   ) :: nblue
    INTEGER,                         INTENT(IN   ) :: mblue
    REAL(MK),                        INTENT(IN   ) :: cutoff
    INTEGER, DIMENSION(:),  pointer, INTENT(INOUT) :: nvlist_cross
    INTEGER, DIMENSION(:,:),POINTER, INTENT(INOUT) :: vlist_cross
    INTEGER,                         INTENT(  OUT)   :: info

    ! local variable
    INTEGER                                             :: ip,iq,isub,iinter
    INTEGER                                             :: ipart,jpart
    INTEGER                                             :: ibox,jbox,cbox
    INTEGER                                             :: n1,n2,n3,i,j,k
    INTEGER                                             :: maxvlen
    INTEGER                                             :: ibegin,iend,jbegin,jend
    REAL(MK)                                            :: dist2,cutoff2
    REAL(MK),DIMENSION(2)                               :: cellsize
    INTEGER, DIMENSION(:,:), POINTER                    :: ncells => null()
    TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER  :: clist => null()
    INTEGER,                     DIMENSION(:,:),POINTER :: ind => null()
    INTEGER,                     DIMENSION(:,:),POINTER :: jnd => null()
    INTEGER                                             :: nnd
    INTEGER                                             :: ndim
    REAL(MK), DIMENSION(:,:),POINTER                    :: all_xp=>NULL()

        !---------------------------------------------------------------------!
        ! initialize
        !---------------------------------------------------------------------!
        ndim = 2
        info = 0
        !---------------------------------------------------------------------!
        ! first fill an array with the new and old particles
        !---------------------------------------------------------------------!
        ALLOCATE(all_xp(ndim,mblue+nred),STAT=info)
        IF (ASSOCIATED(nvlist_cross)) DEALLOCATE(nvlist_cross)
        ALLOCATE(nvlist_cross(nred),STAT=info)
        IF (ASSOCIATED(vlist_cross)) DEALLOCATE(vlist_cross)

        all_xp(1:ndim,1:nred) = xp_red(1:ndim,1:nred)
        all_xp(1:ndim,nred+1:nred+mblue) = xp_blue(1:ndim,1:mblue)

        !!--------------------------------------------------------------------
        !! create cell lists with maximum cutoff
        !!--------------------------------------------------------------------
        cellsize = cutoff

        CALL ppm_neighlist_clist(topo_id,all_xp,mblue+nred,cellsize,&
            .FALSE.,clist,ncells,info)

        !!--------------------------------------------------------------------
        !! create the index list of cell-cell interactons
        !!--------------------------------------------------------------------
        CALL ppm_neighlist_MKneighidx(.FALSE.,ind,jnd,nnd,info)

        !!--------------------------------------------------------------------
        !! run over cells and create verlet lists for the particles inside
        !!--------------------------------------------------------------------

        ! initialise nvlist_cross to zero
        DO ip = 1,nred
            nvlist_cross(ip) = 0
        ENDDO

        !!---------------------------------------------------------------------!
        !!  determine size of verlet lists
        !!---------------------------------------------------------------------!
        DO isub=1,topo%nsublist
            n1  = ncells(1,isub)
            n2  = ncells(1,isub) * ncells(2,isub)
#ifdef __3D
            n3  = ncells(3,isub)
#endif
            ! loop over all cells
#ifdef __3D
            DO k=0,n3-1
#endif
                DO j=0,ncells(2,isub)-1
                    DO i=0,ncells(1,isub)-1

                        ! index of the center box
#ifdef __3D
                        cbox = i + 1 + n1*j + n2*k
#else
                        cbox = i + 1 + n1*j
#endif
                    DO iinter=1,nnd ! loop over all box-box interactions

                        ! determine box indices for this interaction
#ifdef __3D
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter) + &
                            n2*jnd(3,iinter))
#else
                        jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
#endif


                        !special treatment for ghost cells
                        if(jbox.le.0) CYCLE
                        if(jbox.GT.(product(ncells(1:ndim,isub))) ) CYCLE

                        !  get pointers to first and last particle
                        ibegin = clist(isub)%lhbx(cbox)
                        iend   = clist(isub)%lhbx(cbox+1)-1
                        IF (iend .LT. ibegin) CYCLE

                        ! get pointers to first and last particle
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        IF (jend .LT. jbegin) CYCLE

                        ! loop over all particles inside this cell
                        ipart_loop: DO ipart=ibegin,iend
                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart)
                            ! only consider the new (real) particles as ip
                            IF (ip .GT. nred) CYCLE
                            ! check against all particles in the other cell
                            jpart_loop: DO jpart=jbegin,jend
                                IF (jpart .ne. ipart) then
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart)

                                    IF (iq .le. nred) then
                                        ! iq belongs to the new set
                                    else
                                        ! iq belongs to the old set
                                        dist2= sum((all_xp(1:ndim,ip)- &
                                            all_xp(1:ndim,iq))**2)
                                            cutoff2 = rcp_blue(iq-nred)**2

                                        IF (dist2 .le. cutoff2) then
                                            ! populate nvlist_cross
                                            nvlist_cross(ip) = &
                                                nvlist_cross(ip) + 1
                                        ENDIF
                                    ENDIF
                                ENDIF

                            ENDDO jpart_loop
                        ENDDO ipart_loop

                    ENDDO          ! iinter

                ENDDO             ! i
            ENDDO                ! j
#ifdef __3D
        ENDDO   ! k
#endif
        ENDDO                      ! isub

        !!--------------------------------------------------------------------!
        !! allocate verlet list length
        !!--------------------------------------------------------------------!
        maxvlen = maxval(nvlist_cross)
        ALLOCATE(vlist_cross(maxvlen,nred),STAT=info)

        ! initialise nvlist_cross to zero
        DO ip = 1,nred
            nvlist_cross(ip) = 0
        ENDDO

        !!--------------------------------------------------------------------!
        !! fill verlet lists
        !!--------------------------------------------------------------------!
        DO isub=1,topo%nsublist
            n1  = ncells(1,isub)
            n2  = ncells(1,isub)*ncells(2,isub)
#ifdef __3D
            n3  = ncells(3,isub)
#endif
        ! loop over all real cells (the -2 at the end does this)
#ifdef __3D
        DO k=1,n3-2
#endif
            DO j=1,ncells(2,isub)-2
                DO i=1,ncells(1,isub)-2
                    ! index of the center box
#ifdef __3D
                    cbox = i + 1 + n1*j + n2*k
#else
                    cbox = i + 1 + n1*j
#endif

                    DO iinter=1,nnd ! loop over all box-box interactions

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
                        IF (iend .LT. ibegin) CYCLE

                        ! get pointers to first and last particle
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        IF (jend .LT. jbegin) CYCLE

                        ! loop over all particles inside this cell
                        ipart_loop2: DO ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart)
                            ! only consider the new (real) particles as ip
                            IF (ip .GT. nred) CYCLE

                            ! check against all particles in the other cell
                            DO jpart=jbegin,jend
                                IF (jpart .ne. ipart) then
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart)

                                    IF (iq .le. nred) then
                                        ! iq belongs to the new set
                                    else
                                        ! iq belongs to the old set
                                        dist2= sum((all_xp(1:ndim,ip)- &
                                            all_xp(1:ndim,iq))**2)
                                            cutoff2 = rcp_blue(iq-nred)**2

                                        IF (dist2 .le. cutoff2) then
                                            !populate nvlist_cross
                                            nvlist_cross(ip) = &
                                                nvlist_cross(ip) + 1
                                            !vlist contains index of old particles
                                            vlist_cross(nvlist_cross(ip),ip) = &
                                                iq - nred
                                        ENDIF
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDDO ipart_loop2
                    ENDDO          ! iinter
                ENDDO             ! i
            ENDDO                ! j
#ifdef __3D
        ENDDO                ! k
#endif
        ENDDO                      ! isub


        9999 CONTINUE ! jump here upon error
        DEALLOCATE(all_xp)

end subroutine old_xset_vlist_2d

end module old_xset

program ppm_test_xset
!-------------------------------------------------------------------------
!     test case   :   inhomogeneous neighbour lists - 2d
!                           particles placed uniformly randomly in a box
!                           cutoff radii are uniformly random
!-------------------------------------------------------------------------

USE ppm_module_typedef
USE ppm_module_mktopo
USE ppm_module_topo_get
USE ppm_module_init
USE ppm_module_finalize
USE ppm_module_core
USE ppm_module_inl_xset_vlist
USE old_xset

IMPLICIT NONE
#include "../../ppm_define.h"
#ifdef __MPI
INCLUDE 'mpif.h'
#endif

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: mk = ppm_kind_double
REAL(MK),parameter              :: pi = 3.1415926535897931_MK
REAL(MK),parameter              :: skin = 0.0_MK
INTEGER,parameter               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
REAL(MK)                        :: tol,min_rcp,max_rcp
INTEGER                         :: info,comm,rank
INTEGER                         :: topoid
INTEGER                         :: nred = 10000
INTEGER                         :: mred
INTEGER                         :: nblue = 3000
INTEGER                         :: mblue
INTEGER                         :: newnp
REAL(MK),DIMENSION(:,:),POINTER :: xp_red=>NULL()
REAL(MK),DIMENSION(:,:),POINTER :: xp_blue=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: rcp_blue=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: min_phys=>NULL(),max_phys,h,p_h=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: len_phys=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: ghostlayer=>NULL()
INTEGER, DIMENSION(:  ),POINTER :: ghostsize=>NULL()
INTEGER                         :: i,j,sum1,sum2
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost=>NULL()
INTEGER, DIMENSION(:  ),POINTER :: nm=>NULL()
TYPE(ppm_t_topo), POINTER       :: topo
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
REAL(MK), DIMENSION(:),ALLOCATABLE :: randnb
INTEGER,DIMENSION(:,:),POINTER   :: vlist=>NULL(),vlist2=>NULL()
INTEGER,DIMENSION(:),  pointer   :: nvlist=>NULL(),nvlist2=>NULL()
INTEGER                          :: isymm = 0
logical                          :: lsymm = .FALSE.,ok

!----------------
! setup
!----------------
tol = 10.0_MK*EPSILON(1.0_MK)
tolexp = INT(LOG10(EPSILON(1.0_MK)))
min_rcp = 0.01_MK
max_rcp = 0.1_MK

ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &         ghostsize(ndim),ghostlayer(2*ndim),&
    &         nm(ndim),h(ndim),p_h(ndim),STAT=info)

min_phys(1:ndim) = 0.0_MK
max_phys(1:ndim) = 1.0_MK
len_phys(1:ndim) = max_phys-min_phys
ghostsize(1:ndim) = 2
ghostlayer(1:2*ndim) = max_rcp
bcdef(1:6) = ppm_param_bcdef_periodic


#ifdef __MPI
comm = MPI_COMM_WORLD
CALL MPI_Init(info)
CALL MPI_Comm_rank(comm,rank,info)
#else
rank = 0
#endif
call ppm_init(ndim,mk,tolexp,0,debug,info,99)

CALL RANDOM_SEED(size=seedsize)
ALLOCATE(seed(seedsize))
ALLOCATE(randnb((ndim)*nred),STAT=info)
DO i=1,seedsize
    seed(i)=10+i*i*(rank+1)
ENDDO
CALL RANDOM_SEED(put=seed)
CALL RANDOM_NUMBER(randnb)

!----------------
! create particles
!----------------

ALLOCATE(xp_red(ndim,nred),STAT=info)
ALLOCATE(xp_blue(ndim,nblue),rcp_blue(nblue),STAT=info)
xp_red = 0.0_MK
xp_blue = 0.0_MK
rcp_blue = 0.0_MK

DO i=1,nred
    DO j=1,2
        xp_red(j,i) = min_phys(j)+&
            len_phys(j)*randnb(2*i-j+1)
    ENDDO
ENDDO

DEALLOCATE(randnb)
ALLOCATE(randnb((1+ndim)*nblue),STAT=info)
CALL RANDOM_NUMBER(randnb)
DO i=1,nblue
    DO j=1,ndim
        xp_blue(j,i) = min_phys(j)+&
            len_phys(j)*randnb((ndim+1)*i-(ndim-j))
    ENDDO
    rcp_blue(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
ENDDO

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
DO i=1,nred
    IF (nvlist(i).ne.nvlist2(i)) then
        write(*,'(a)') '!! --------------- !!'
        write(*,'(a)') '!! failed'
        print *, '!!    nvlist are not equal for ip = ',i
        print * ,'!! nvlist:  ',nvlist(i)
        print * ,'!! nvlist2: ',nvlist2(i)
        write(*,'(a)') '!! --------------- !!'
        goto 8000
    ENDIF
    sum1=0;sum2=0;
    DO j=1,nvlist(i)
        sum1 = sum1 + vlist(j,i)**2
        sum2 = sum2 + vlist2(j,i)**2
    ENDDO
    IF (sum1 .ne. sum2) then
        write(*,'(a)') '!! --------------- !!'
        write(*,'(a)') '!! failed'
        print *, '!!    vlist are not equal for ip = ',i
        write(*,'(a)') '!! --------------- !!'
        goto 8000
    ENDIF
ENDDO


if (rank.EQ.0) then
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
8000 CONTINUE

call ppm_finalize(info)
#ifdef __MPI
call mpi_finalize(info)
#endif

!DEALLOCATE(xp_red,min_phys,max_phys,len_phys,ghostsize,nm)
!DEALLOCATE(xp_blue,rcp_blue,nvlist,nvlist2,vlist,vlist2)

if (rank.EQ.0) print *, 'done.'


END PROGRAM ppm_test_xset
