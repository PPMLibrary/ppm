!!!----------------------------------------------------------------------------!
!!! Routine still in testing mode
!!!!
!!! creates Verlet lists for particles with different cutoffs
!!! WARNING: the lists are re-allocated only when they get bigger
!!! this can look like a small memory leak...
!!!
!!! The preprocessor flag sym_lists can take the values:
!!! 1 : neighbour lists are symmetric; cutoff = MAX(rcp(ip),rcp(iq))
!!! 2 : neighbour lists are NOT symmetric; cutoff = rcp(ip)
!!! 3 : neighbour lists are NOT symmetric; cutoff = rcp(iq)
!!! 4 : neighbour lists are symmetric; cutoff = MIN(rcp(ip),rcp(iq))
!!!
!!!
!!! remark: expect cutoff to be up-to-date
!!!
!!!
!!! If ensure_size is present and TRUE, then the neighbour lists are
!!! rebuilt until they are all large enough to allow accurate computation
!!! of the DC operators (not very efficient...)
!!!----------------------------------------------------------------------------!

    SUBROUTINE cnl_vlist(xp,rcp,Npart,Mpart,xmin,xmax,nvlist,vlist,ndim,&
 &                             info)


    USE ppm_module_typedef
    IMPLICIT NONE

#define sym_lists 4
!#define __3D

      INTEGER, PARAMETER :: MK = ppm_kind_double
    ! arguments
    REAL(MK), DIMENSION(:,:),POINTER, INTENT(IN   )   :: xp
    REAL(MK), DIMENSION(:),POINTER,   INTENT(INOUT)   :: rcp
    INTEGER,                          INTENT(IN   )   :: Npart
    INTEGER,                          INTENT(IN   )   :: Mpart
    REAL(MK), DIMENSION(ndim),        INTENT(IN   )   :: xmin,xmax
    INTEGER,  DIMENSION(:),  POINTER, INTENT(INOUT)   :: nvlist
    INTEGER,  DIMENSION(:,:),POINTER, INTENT(INOUT)   :: vlist
    INTEGER,                          INTENT(  OUT)   :: info
    INTEGER,                          INTENT(IN   )   :: ndim


    ! local variable
    INTEGER                           :: ip,iq,iinter,ineigh
    INTEGER                           :: ipart,jpart
    INTEGER                           :: ibox,jbox,cbox
    INTEGER                           :: n1,n2,n3,i,j,k
    INTEGER                           :: maxvlen
    REAL(MK)                          :: cellsize_1d
    INTEGER                           :: ibegin,iend,jbegin,jend
    CHARACTER(LEN=256)                :: cbuf
    REAL(MK)                          :: dist2,cutoff2,cutoff
    REAL(MK)                          :: verlet_skin
    REAL(MK),DIMENSION(ndim)          :: cellsize
    REAL(KIND(1.D0))                  :: t0
    REAL(MK)                          :: minrcp
    REAL(MK), DIMENSION(:),POINTER    :: cutoff_incr
    LOGICAL                           :: first_loop
    INTEGER, DIMENSION(:,:), POINTER       :: ind
    !!! First interaction partner (box which *interacts*).
    !!!
    !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
    !!! 2nd index: interaction number 1...nnd.
    INTEGER, DIMENSION(:,:), POINTER       :: jnd
    !!! Second interaction partner (box which *is interacted with*).
    !!! 
    !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
    !!! 2nd index: interaction number 1...nnd.
    INTEGER                                :: nnd
    !!! Number of box-box interactions to be performed.
    TYPE(t_clist), DIMENSION(1)  :: clist
    REAL(MK)                          :: ts,te
    REAL(MK)                          :: tcs,tce
    REAL(MK)                          :: tvs,tve

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
    NULLIFY(ind)
    NULLIFY(jnd)
    NULLIFY(cutoff_incr)
    
    !!-------------------------------------------------------------------------!
    !! Create cell lists with maximum cutoff
    !!-------------------------------------------------------------------------!
    verlet_skin = 1._MK
    first_loop=.TRUE.
    cutoff = maxval(rcp)
    cellsize = cutoff*verlet_skin
    CALL cnl_clist(xp,Mpart,xmin,xmax,cellsize,&
        .FALSE.,clist,ndim,info)
    IF (info .NE. 0) THEN
        print *, 'cnl_clist failed'
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Create the index list of cell-cell interactons
    !!-------------------------------------------------------------------------!
    CALL cnl_MkNeighIdx(.FALSE.,ind,jnd,nnd,ndim,info)
    IF (info .NE. 0) THEN
        print *, 'cnl_MkNeighIdx failed'
        info = -1
        GOTO 9999
    ENDIF
    !!-------------------------------------------------------------------------!
    !! Run over cells and create Verlet lists for the particles inside
    !!-------------------------------------------------------------------------!
    !reallocate nvlist only if the number of particles has changed
    IF (ASSOCIATED(nvlist)) THEN
        IF (SIZE(nvlist).LT.Npart) THEN
            DEALLOCATE(nvlist)
            ALLOCATE(nvlist(Npart),STAT=info)
        ENDIF
    ELSE
        ALLOCATE(nvlist(Npart),STAT=info)
    ENDIF
    IF (info .NE. 0) THEN
        print *,'Allocation of nvlist failed.'
        info = -1
        GOTO 9999
    ENDIF

    ! reset nvlist
    DO ip = 1,Npart
        nvlist(ip) = 0
    ENDDO


    !!-------------------------------------------------------------------------!
    !! Determine size of Verlet lists
    !!-------------------------------------------------------------------------!
    n1  = clist(1)%Nm(1)
    n2  = clist(1)%Nm(1)*clist(1)%Nm(2)
#ifdef __3D
    n3  = clist(1)%Nm(3)
#endif
    ! loop over ALL cells
#ifdef __3D
    DO k=0,n3-1
#endif
        DO j=0,clist(1)%Nm(2)-1
            DO i=0,clist(1)%Nm(1)-1

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
                    IF(jbox.LE.0) CYCLE
                    IF(jbox.GT.(PRODUCT(clist(1)%Nm(1:ndim))) ) CYCLE

                    !  get pointers to first and last particle 
                    ibegin = clist(1)%lhbx(cbox)
                    iend   = clist(1)%lhbx(cbox+1)-1
                    IF (iend .LT. ibegin) CYCLE

                    ! Within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    IF (cbox .EQ. jbox) THEN
                        DO ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(1)%lpdx(ipart) 

                            IF (ip .LE. Npart) THEN
                                DO jpart=ibegin,iend

                                    IF (jpart .NE. ipart) THEN
                                        ! translate to real particle index
                                        iq = clist(1)%lpdx(jpart) 

                                        ! if not a ghost cell
                                        dist2=&
                                            SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 

#if sym_lists==1
                                        cutoff2 = (verlet_skin*&
                                            MAX(rcp(iq),rcp(ip)))**2
#elif sym_lists==2
                                        cutoff2 = (verlet_skin*rcp(ip))**2
#elif sym_lists==3
                                        cutoff2 = (verlet_skin*rcp(iq))**2
#else
                                        IF(first_loop) THEN
                                            cutoff2 = (verlet_skin*&
                                                MIN(rcp(iq),rcp(ip)))**2
                                        ELSE
                                            cutoff2 = (cutoff_incr(ip)+&
                                                verlet_skin*&
                                                MIN(rcp(iq),rcp(ip)))**2
                                        ENDIF
#endif
                                        IF (dist2 .LE. cutoff2) THEN
                                            ! add particle iq to 
                                            !list of particle ip
                                            nvlist(ip) = nvlist(ip) + 1
                                        ENDIF
                                    ENDIF
                                ENDDO
                            ENDIF

                        ENDDO

                        !  For the other boxes check all particles
                    ELSE

                        ! get pointers to first and last particle 
                        jbegin = clist(1)%lhbx(jbox)
                        jend   = clist(1)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        If (jend .LT. jbegin) CYCLE

                        ! loop over all particles inside this cell
                        DO ipart=ibegin,iend 

                            ! translate to real particle index
                            ip = clist(1)%lpdx(ipart) 

                            IF (ip .LE. Npart) THEN
                                ! check against all particles in the other cell
                                DO jpart=jbegin,jend

                                    ! translate to real particle index
                                    iq = clist(1)%lpdx(jpart) 

                                    dist2=&
                                        SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
#if sym_lists==1
                                    cutoff2 = (verlet_skin*&
                                        MAX(rcp(iq),rcp(ip)))**2
#elif sym_lists==2
                                    cutoff2 = (verlet_skin*rcp(ip))**2
#elif sym_lists==3
                                    cutoff2 = (verlet_skin*rcp(iq))**2
#else
                                    IF(first_loop) THEN
                                        cutoff2 = (verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
                                    ELSE
                                        cutoff2 = (cutoff_incr(ip)+&
                                            verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
                                    ENDIF
#endif

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
#ifdef __3D
    ENDDO                ! k
#endif

    !!-------------------------------------------------------------------------!
    !! Allocate Verlet list length
    !!-------------------------------------------------------------------------!
    maxvlen = MAXVAL(nvlist)

    IF (ASSOCIATED(vlist)) THEN
        IF(SIZE(vlist,1).LT.maxvlen .OR. SIZE(vlist,2).LT.Npart) THEN
            print *,'reallocating vlist'
            DEALLOCATE(vlist)
            ALLOCATE(vlist(maxvlen,Npart),STAT=info)
        ENDIF
    ELSE
        ALLOCATE(vlist(maxvlen,Npart),STAT=info)
    ENDIF

    IF (info .NE. 0) THEN
        print *,'Allocation of vlist failed.'
        info = -1
        GOTO 9999
    ENDIF

    ! reset nvlist
    DO ip = 1,SIZE(nvlist)
        nvlist(ip) = 0
    ENDDO

    !!-------------------------------------------------------------------------!
    !! Fill Verlet lists
    !!-------------------------------------------------------------------------!
    n1  = clist(1)%Nm(1)
    n2  = clist(1)%Nm(1)*clist(1)%Nm(2)
#ifdef __3D
    n3  = clist(1)%Nm(3)
#endif
    ! loop over all REAL cells (the -2 at the end does this)
#ifdef __3D
    DO k=1,n3-2
#endif
        DO j=1,clist(1)%Nm(2)-2
            DO i=1,clist(1)%Nm(1)-2
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
                    ibegin = clist(1)%lhbx(cbox)
                    iend   = clist(1)%lhbx(cbox+1)-1
                    IF (iend .LT. ibegin) CYCLE

                    ! Within the box itself use symmetry and avoid adding 
                    ! the particle itself to its own list
                    IF (cbox .EQ. jbox) THEN
                        DO ipart=ibegin,iend

                            ! translate to real particle index
                            ip = clist(1)%lpdx(ipart) 

                            IF (ip .LE. Npart) THEN
                                DO jpart=ibegin,iend

                                    IF (jpart .NE. ipart) THEN
                                        ! translate to real particle index
                                        iq = clist(1)%lpdx(jpart) 

                                        dist2=SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 

#if sym_lists==1
                                        cutoff2 = (verlet_skin*MAX(rcp(iq),rcp(ip)))**2
#elif sym_lists==2
                                        cutoff2 = (verlet_skin*rcp(ip))**2
#elif sym_lists==3
                                        cutoff2 = (verlet_skin*rcp(iq))**2
#else
                                        IF(first_loop) THEN
                                            cutoff2 = (verlet_skin*&
                                                MIN(rcp(iq),rcp(ip)))**2
                                        ELSE
                                            cutoff2 = (cutoff_incr(ip)+&
                                                verlet_skin*MIN(rcp(iq),rcp(ip)))**2
                                        ENDIF
#endif
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

                        !  For the other boxes check all particles
                    ELSE

                        ! get pointers to first and last particle 
                        jbegin = clist(1)%lhbx(jbox)
                        jend   = clist(1)%lhbx(jbox+1)-1

                        ! skip this iinter if empty
                        If (jend .LT. jbegin) CYCLE

                        DO ipart=ibegin,iend ! loop over all particles inside this cell

                            ! translate to real particle index
                            ip = clist(1)%lpdx(ipart) 

                            IF (ip .LE. Npart) THEN
                                ! check against all particles in the other cell
                                DO jpart=jbegin,jend

                                    ! translate to real particle index
                                    iq = clist(1)%lpdx(jpart) 

                                    dist2=&
                                        SUM((xp(1:ndim,ip)-xp(1:ndim,iq))**2) 
#if sym_lists==1
                                    cutoff2 = (verlet_skin*&
                                        MAX(rcp(iq),rcp(ip)))**2
#elif sym_lists==2
                                    cutoff2 = (verlet_skin*rcp(ip))**2
#elif sym_lists==3
                                    cutoff2 = (verlet_skin*rcp(iq))**2
#else
                                    IF(first_loop) THEN
                                        cutoff2 = (verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
                                    ELSE
                                        cutoff2 = (cutoff_incr(ip)+&
                                            verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
                                    ENDIF
#endif
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
#ifdef __3D
    ENDDO                ! k
#endif

#if debug_verbosity > 1
    WRITE(cbuf,'(A,F8.2,2(A,I6),1X,L)')'Avg num of neighbors: ',&
        REAL(SUM(nvlist(1:Npart)),MK)/REAL(Npart,MK),&
        ' Max:',MAXVAL(nvlist(1:Npart)), ' Min:',MINVAL(nvlist(1:Npart)),&
        first_loop
    print *,cbuf

#endif

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
    DEALLOCATE(clist(1)%Nm)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE cnl_vlist
