#if __DIM == 2
SUBROUTINE DTYPE(cnl_vlist_build_2d)(xp,rcp,Npart,Mpart,xmin,xmax,nvlist,vlist,&
           ind,jnd,nnd,clist,ndim,verlet_skin,info)
#endif
#if __DIM == 3
SUBROUTINE DTYPE(cnl_vlist_build_3d)(xp,rcp,Npart,Mpart,xmin,xmax,nvlist,vlist,&
           ind,jnd,nnd,clist,ndim,verlet_skin,info)
#endif


    USE ppm_module_typedef
    IMPLICIT NONE

#define sym_lists 4

    DEFINE_MK()
    ! arguments
    REAL(MK), DIMENSION(:,:),POINTER, INTENT(IN   )   :: xp
    REAL(MK), DIMENSION(:),POINTER,   INTENT(IN   )   :: rcp
    INTEGER,                          INTENT(IN   )   :: Npart
    INTEGER,                          INTENT(IN   )   :: Mpart
    REAL(MK), DIMENSION(ndim),        INTENT(IN   )   :: xmin,xmax
    INTEGER,  DIMENSION(:),  POINTER, INTENT(INOUT)   :: nvlist
    INTEGER,  DIMENSION(:,:),POINTER, INTENT(INOUT)   :: vlist
    INTEGER, DIMENSION(:,:), POINTER, INTENT(IN   )   :: ind
    !!! First interaction partner (box which *interacts*).
    INTEGER, DIMENSION(:,:), POINTER, INTENT(IN   )   :: jnd
    !!! Second interaction partner (box which *is interacted with*).
    !!! 
    !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
    !!! 2nd index: interaction number 1...nnd.
    INTEGER,                          INTENT(IN   )   :: nnd
    !!! Number of box-box interactions to be performed.
    TYPE(t_clist), DIMENSION(1),      INTENT(IN   )   :: clist
    INTEGER,                          INTENT(IN   )   :: ndim
    REAL(MK),                         INTENT(IN   )   :: verlet_skin
    INTEGER,                          INTENT(  OUT)   :: info


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
    REAL(MK),DIMENSION(ndim)          :: cellsize
    REAL(KIND(1.D0))                  :: t0
    REAL(MK)                          :: minrcp
    !!!
    !!! 1st index: 1...3 (x,y,[z]) index shift.                              +
    !!! 2nd index: interaction number 1...nnd.
    REAL(MK)                          :: ts,te
    REAL(MK)                          :: tcs,tce
    REAL(MK)                          :: tvs,tve


    info = 0
    !!-------------------------------------------------------------------------!
    !! Determine size of Verlet lists
    !!-------------------------------------------------------------------------!
    n1  = clist(1)%Nm(1)
    n2  = clist(1)%Nm(1)*clist(1)%Nm(2)
#if __DIM == 3
    n3  = clist(1)%Nm(3)
#endif
    ! loop over ALL cells
#if __DIM == 3
    DO k=0,n3-1
#endif
        DO j=0,clist(1)%Nm(2)-1
            DO i=0,clist(1)%Nm(1)-1

                ! index of the center box
#if __DIM == 3
                cbox = i + 1 + n1*j + n2*k
#else
                cbox = i + 1 + n1*j
#endif
                DO iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
#if __DIM == 3
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
                                        cutoff2 = (verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
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
                                    cutoff2 = (verlet_skin*&
                                        MIN(rcp(iq),rcp(ip)))**2
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
#if __DIM == 3
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
#if __DIM == 3
    n3  = clist(1)%Nm(3)
#endif
    ! loop over all REAL cells (the -2 at the end does this)
#if __DIM == 3
    DO k=1,n3-2
#endif
        DO j=1,clist(1)%Nm(2)-2
            DO i=1,clist(1)%Nm(1)-2
                ! index of the center box
#if __DIM == 3
                cbox = i + 1 + n1*j + n2*k
#else
                cbox = i + 1 + n1*j
#endif

                DO iinter=1,nnd ! loop over all box-box interactions

                    ! determine box indices for this interaction
#if __DIM == 3
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
                                        cutoff2 = (verlet_skin*&
                                            MIN(rcp(iq),rcp(ip)))**2
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
                                    cutoff2 = (verlet_skin*&
                                        MIN(rcp(iq),rcp(ip)))**2
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
#if __DIM == 3
    ENDDO                ! k
#endif

    9999 CONTINUE ! jump here upon error

#if __DIM == 2
END SUBROUTINE DTYPE(cnl_vlist_build_2d)
#endif
#if __DIM == 3
END SUBROUTINE DTYPE(cnl_vlist_build_3d)
#endif

#undef __DIM

