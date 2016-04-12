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
IMPLICIT NONE
contains

subroutine old_inl_vlist_2d(topo_id,topo,xp,rcp,cutoff,npart,mpart,&
        nvlist,vlist,info)


    USE ppm_module_typedef
    USE ppm_module_neighlist
    IMPLICIT NONE


    ! arguments
    INTEGER, PARAMETER              :: mk = ppm_kind_double
    INTEGER,                         INTENT(IN   )   :: topo_id
    TYPE(ppm_t_topo), POINTER,       INTENT(IN   )   :: topo
    REAL(MK),DIMENSION(:,:),POINTER, INTENT(IN   )   :: xp
    REAL(MK),DIMENSION(:),POINTER,   INTENT(INOUT)   :: rcp
    REAL(MK),                        INTENT(IN   )   :: cutoff
    INTEGER,                         INTENT(IN   )   :: npart
    INTEGER,                         INTENT(IN   )   :: mpart
    INTEGER, DIMENSION(:),  pointer, INTENT(INOUT)   :: nvlist
    INTEGER, DIMENSION(:,:),POINTER, INTENT(INOUT)   :: vlist
    INTEGER,                         INTENT(  OUT)   :: info

    ! local variable
    INTEGER                           :: ip,iq,isub,iinter
    INTEGER                           :: ipart,jpart
    INTEGER                           :: jbox,cbox
    INTEGER                           :: n1,n2,i,j
    INTEGER                           :: maxvlen
    INTEGER                           :: ibegin,iend,jbegin,jend
    REAL(MK)                          :: dist2,cutoff2
    REAL(MK),DIMENSION(2)          :: cellsize
    INTEGER, DIMENSION(:,:), POINTER              :: ncells => null()
    TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER  :: clist => null()
    INTEGER,                     DIMENSION(:,:),POINTER :: ind => null()
    INTEGER,                     DIMENSION(:,:),POINTER :: jnd => null()
    INTEGER                                             :: nnd
    INTEGER                                             :: ndim = 2


    !!-------------------------------------------------------------------------!
    !! initialize
    !!-------------------------------------------------------------------------!
    info = 0

    !!-------------------------------------------------------------------------!
    !! create cell lists with maximum cutoff
    !!-------------------------------------------------------------------------!
    cellsize = cutoff
    CALL ppm_neighlist_clist(topo_id,xp,mpart,cellsize,&
        .FALSE.,clist,ncells,info)

    !!-------------------------------------------------------------------------!
    !! create the index list of cell-cell interactons
    !!-------------------------------------------------------------------------!
    CALL ppm_neighlist_MKneighidx(.FALSE.,ind,jnd,nnd,info)

    !!-------------------------------------------------------------------------!
    !! run over cells and create verlet lists for the particles inside
    !!-------------------------------------------------------------------------!

    !reallocate nvlist only if the number of particles has changed
    IF (ASSOCIATED(nvlist)) THEN
        IF (size(nvlist).LT.npart) THEN
            DEALLOCATE(nvlist)
            ALLOCATE(nvlist(npart),STAT=info)
        ENDIF
    else
        ALLOCATE(nvlist(npart),STAT=info)
    ENDIF

    ! reset nvlist
    DO ip = 1,npart
        nvlist(ip) = 0
    ENDDO

    !!-------------------------------------------------------------------------!
    !! determine size of verlet lists
    !!-------------------------------------------------------------------------!
    DO isub=1,topo%nsublist
        n1  = ncells(1,isub)
        n2  = ncells(1,isub)*ncells(2,isub)
        ! loop over all cells
        DO j=0,ncells(2,isub)-1
            DO i=0,ncells(1,isub)-1

                ! index of the center box
                cbox = i + 1 + n1*j
                DO iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
                    jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
                    !special treatment for ghost cells
                    if(jbox.LE.0) CYCLE
                    if(jbox.GT.(PRODUCT(ncells(:,isub))) ) CYCLE

                    !  get pointers to first and last particle 
                    ibegin = clist(isub)%lhbx(cbox)
                    iend   = clist(isub)%lhbx(cbox+1)-1
                    IF (iend .LT. ibegin) CYCLE

                    ! within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    IF (cbox .EQ. jbox) THEN
                        DO ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            IF (ip .LE. npart) THEN
                                DO jpart=ibegin,iend

                                    IF (jpart .NE. ipart) THEN
                                        ! translate to real particle index
                                        iq = clist(isub)%lpdx(jpart) 

                                        ! if not a ghost cell
                                        dist2=&
                                            SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 

                                        cutoff2 = MIN(rcp(iq),rcp(ip))**2
                                        IF (dist2 .LE. cutoff2) THEN
                                            ! add particle iq to 
                                            !list of particle ip
                                            nvlist(ip) = nvlist(ip) + 1
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO

                        !  for the other boxes check all particles
                    else
                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1
                        ! skip this iinter if empty
                        IF (jend .LT. jbegin) CYCLE
                        ! loop over all particles inside this cell
                        DO ipart=ibegin,iend 
                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 
                            IF (ip .LE. npart) THEN
                                ! check against all particles in the other cell
                                DO jpart=jbegin,jend
                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 
                                    dist2=&
                                        SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                    cutoff2 = MIN(rcp(iq),rcp(ip))**2
                                    IF (dist2 .LE. cutoff2) THEN
                                        !add particle 
                                        !iq to list of particle ip
                                        nvlist(ip) = nvlist(ip) + 1
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO

                    ENDIF       ! cbox .EQ. jbox
                ENDDO          ! iinter
            ENDDO             ! i
        ENDDO                ! j
    ENDDO                      ! isub
    !!-------------------------------------------------------------------------!
    !! allocate verlet list length
    !!-------------------------------------------------------------------------!
    maxvlen = MAXVAL(nvlist)
    IF (ASSOCIATED(vlist)) THEN
        if(size(vlist,1).LT.maxvlen .OR. size(vlist,2).LT.npart) THEN
            DEALLOCATE(vlist)
            ALLOCATE(vlist(maxvlen,npart),STAT=info)
        ENDIF
    else
        ALLOCATE(vlist(maxvlen,npart),STAT=info)
    ENDIF

    ! reset nvlist
    DO ip = 1,size(nvlist)
        nvlist(ip) = 0
    ENDDO

    !!-------------------------------------------------------------------------!
    !! fill verlet lists
    !!-------------------------------------------------------------------------!
    DO isub=1,topo%nsublist
        n1  = ncells(1,isub)
        n2  = ncells(1,isub)*ncells(2,isub)
        ! loop over all real cells (the -2 at the end does this)
        DO j=1,ncells(2,isub)-2
            DO i=1,ncells(1,isub)-2
                ! index of the center box
                cbox = i + 1 + n1*j
                DO iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
                    jbox = cbox + (jnd(1,iinter) + n1*jnd(2,iinter))
                    !  get pointers to first and last particle 
                    ibegin = clist(isub)%lhbx(cbox)
                    iend   = clist(isub)%lhbx(cbox+1)-1
                    IF (iend .LT. ibegin) CYCLE

                    ! within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    IF (cbox .EQ. jbox) THEN
                        DO ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            IF (ip .LE. npart) THEN
                                DO jpart=ibegin,iend

                                    IF (jpart .NE. ipart) THEN
                                        ! translate to real particle index
                                        iq = clist(isub)%lpdx(jpart) 

                                        dist2=&
                                            SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                        cutoff2 = MIN(rcp(iq),rcp(ip))**2
                                        IF (dist2 .LE. cutoff2) THEN
                                            ! add particle iq to 
                                            ! list of particle ip
                                            nvlist(ip) = nvlist(ip) + 1
                                            vlist(nvlist(ip),ip) = iq
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDIF
                        ENDDO

                        !  for the other boxes check all particles
                    else

                        ! get pointers to first and last particle 
                        jbegin = clist(isub)%lhbx(jbox)
                        jend   = clist(isub)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        IF (jend .LT. jbegin) CYCLE

                        DO ipart=ibegin,iend ! loop over all particles inside this cell

                            ! translate to real particle index
                            ip = clist(isub)%lpdx(ipart) 

                            IF (ip .LE. npart) THEN
                                ! check against all particles in the other cell
                                DO jpart=jbegin,jend

                                    ! translate to real particle index
                                    iq = clist(isub)%lpdx(jpart) 

                                    dist2=&
                                        SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
                                    cutoff2 = MIN(rcp(iq),rcp(ip))**2
                                    IF (dist2 .LE. cutoff2) THEN
                                        ! add particle iq to list 
                                        ! of particle ip
                                        nvlist(ip) = nvlist(ip) + 1
                                        vlist(nvlist(ip),ip) = iq
                                    ENDIF

                                ENDDO
                            ENDIF
                        ENDDO
                    ENDIF       ! cbox .EQ. jbox
                ENDDO          ! iinter

            ENDDO             ! i
        ENDDO                ! j
    ENDDO                      ! isub
    9999 CONTINUE ! jump here upon error

END SUBROUTINE old_inl_vlist_2d

END MODULE old_inl

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
INTEGER                         :: np,npgrid = 123
INTEGER                         :: mp
INTEGER                         :: newnp
INTEGER, DIMENSION(:),  pointer :: p_i
REAL(MK),DIMENSION(:,:),POINTER :: xp
REAL(MK),DIMENSION(:  ),POINTER :: rcp
REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,h,p_h
REAL(MK),DIMENSION(:  ),POINTER :: len_phys
REAL(MK),DIMENSION(:  ),POINTER :: ghostlayer
REAL(MK),DIMENSION(:,:),POINTER :: minsub,maxsub
INTEGER, DIMENSION(:  ),POINTER :: sub2proc
INTEGER                         :: ii,i,j,sum1,sum2
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost
TYPE(ppm_t_topo), POINTER       :: topo
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
REAL(MK), DIMENSION(:),ALLOCATABLE :: randnb
INTEGER,DIMENSION(:,:),POINTER   :: vlist,vlist2
INTEGER,DIMENSION(:),  pointer   :: nvlist,nvlist2
INTEGER                          :: isymm = 0
logical                          :: lsymm = .FALSE.,ok

!----------------
! setup
!----------------
tol = 10.0_MK*EPSILON(1.0_MK)
tolexp = INT(LOG10(EPSILON(1.0_MK)))
min_rcp = 0.01_MK
max_rcp = 0.1_MK
np=npgrid*npgrid

ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
    &         ghostlayer(2*ndim),&
    &         h(ndim),p_i(ndim),p_h(ndim),STAT=info)

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

CALL RANDOM_SEED(SIZE=seedsize)
ALLOCATE(seed(seedsize))
ALLOCATE(randnb((ndim+1)*np),STAT=info)
DO i=1,seedsize
    seed(i)=10+i*i*(rank+1)
ENDDO
CALL RANDOM_SEED(PUT=seed)
CALL RANDOM_NUMBER(randnb)

!----------------
! create particles on a grid
!----------------

ALLOCATE(xp(ndim,np),rcp(np),STAT=info)
xp = 0.0_MK
rcp = 0.0_MK

p_h = len_phys / REAL(npgrid,mk)

DO j=1,npgrid
    DO i=1,npgrid
        p_i = i + (j-1)*npgrid
        xp(1,p_i) = min_phys(1)+REAL(i-1,mk)*p_h(1)
        xp(2,p_i) = min_phys(2)+REAL(j-1,mk)*p_h(2)
        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
    ENDDO
ENDDO

!----------------
! make topology
!----------------
assig  = ppm_param_assign_internal
topoid = 0

decomp = ppm_param_decomp_user_defined
ALLOCATE(minsub(ndim,4),maxsub(ndim,4),sub2proc(4),cost(4),STAT=info)

cost = 1.0_MK; sub2proc = 0
minsub(1:ndim,1)=0.0_MK; maxsub(1:ndim,1)=0.5_MK

minsub(1,2)=0.5_MK; minsub(2,2)=0.0_MK
maxsub(1,2)=1.0_MK; maxsub(2,2)=0.5_MK

minsub(1,3)=0.0_MK; minsub(2,3)=0.5_MK
maxsub(1,3)=0.5_MK; maxsub(2,3)=1.0_MK

minsub(1:ndim,4)=0.5_MK; maxsub(1:ndim,4)=1.0_MK

topoid = 0
call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               max_rcp,cost,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=4,user_sub2proc=sub2proc)
call ppm_topo_get(topoid,topo,info)

!do multiple runs; each time, particle positions are perturbed
!by a different amount (starting from 0)
DO ii=1,5
    CALL RANDOM_NUMBER(randnb)
    DO i=1,np
        xp(1:ndim,i)=xp(1:ndim,i)+&
            0.05_MK*REAL(ii-1,mk)*randnb(ndim*i-(ndim-1):ndim*i)
    ENDDO
    !apply periodic boundary conditions
    DO i=1,np
        DO j=1,ndim
            xp(j,i)=MOD(xp(j,i),1.0_MK)
        ENDDO
    ENDDO

    !map particles onto the topology
    CALL ppm_map_part_global(topoid,xp,np,info)
    CALL ppm_map_part_push(rcp,np,info)
    CALL ppm_map_part_send(np,newnp,info)
    CALL ppm_map_part_pop(rcp,np,newnp,info)
    CALL ppm_map_part_pop(xp,ndim,np,newnp,info)
    np=newnp

    CALL ppm_topo_check(topoid,xp,np,ok,info)
    IF (.NOT. ok) WRITE(*,*) '[',rank,'] topo_check failed'

    CALL ppm_map_part_ghost_get(topoid,xp,ndim,np,isymm,max_rcp,info)
    CALL ppm_map_part_push(rcp,np,info)
    CALL ppm_map_part_send(np,mp,info)
    CALL ppm_map_part_pop(rcp,np,mp,info)
    CALL ppm_map_part_pop(xp,ndim,np,mp,info)

    CALL ppm_dbg_print(topoid,max_rcp,1,1,info,xp,np,mp)


    CALL ppm_inl_vlist(topoid,xp,np,mp,rcp,skin, &
        & lsymm,ghostlayer,info,vlist,nvlist)

    CALL old_inl_vlist_2d(topoid,topo,xp,rcp,max_rcp,np,mp,&
        nvlist2,vlist2,info)

    !compare neighbour lists obtained by the 2 different routines
    DO i=1,np
        IF (nvlist(i).NE.nvlist2(i)) THEN
            WRITE(*,'(A)') '!! --------------- !!'
            WRITE(*,'(A)') '!! FAILED'
            print *, '!!    nvlist are not equal for ip = ',i
            WRITE(*,'(A)') '!! --------------- !!'
            GOTO 8000
        ENDIF
        sum1=0;sum2=0;
        DO j=1,nvlist(i)
            sum1 = sum1 + vlist(j,i)**2
            sum2 = sum2 + vlist2(j,i)**2
        ENDDO
        IF (sum1 .NE. sum2) THEN
            WRITE(*,'(A)') '!! --------------- !!'
            WRITE(*,'(A)') '!! FAILED'
            print *, '!!    vlist are not equal for ip = ',i
            WRITE(*,'(A)') '!! --------------- !!'
            GOTO 8000
        ENDIF
    ENDDO

ENDDO

if (rank.EQ.0) THEN
    WRITE(*,'(A)') '!! --------------- !!'
    WRITE(*,'(A)') '!! SUCCESS'
    WRITE(*,'(2(A,I6),A)') '!!     Completed test for INL with ',np,&
    ' real particles and ',mp-np,' ghosts'
    WRITE(*,'(2(A,E10.4))') '!!     Smallest cutoff was ', MINVAL(rcp(1:np)),&
    ' and largest ',MAXVAL(rcp(1:np))
    WRITE(*,'(2(A,I6))') '!!     Smallest nb of neigh was ',&
    MINVAL(nvlist),' and largest ',MAXVAL(nvlist)
    WRITE(*,'(A)') '!! --------------- !!'
endif

!----------------
! teardown...
!----------------
8000 CONTINUE

call ppm_finalize(info)
#ifdef __MPI
call MPI_finalize(info)
#endif

DEALLOCATE(xp,rcp,min_phys,max_phys,len_phys,minsub,maxsub,cost,sub2proc)

if (rank.EQ.0) print *, 'done.'

END PROGRAM ppm_test_inl
