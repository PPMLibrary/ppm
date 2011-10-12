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
    INTEGER                           :: ip
    INTEGER                           :: i,j,k
    REAL(MK)                          :: cellsize_1d
    CHARACTER(LEN=256)                :: cbuf
    REAL(MK)                          :: cutoff
    REAL(MK)                          :: verlet_skin
    REAL(MK),DIMENSION(ndim)          :: cellsize
    REAL(KIND(1.D0))                  :: t0
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
    
    !!-------------------------------------------------------------------------!
    !! Create cell lists with maximum cutoff
    !!-------------------------------------------------------------------------!
    verlet_skin = 1._MK
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


    IF (ndim .EQ. 2) THEN
        CALL cnl_vlist_build_2d(xp,rcp,Npart,Mpart,xmin,xmax,nvlist,vlist,&
            ind,jnd,nnd,clist,ndim,verlet_skin,info)
    ELSE
        CALL cnl_vlist_build_3d(xp,rcp,Npart,Mpart,xmin,xmax,nvlist,vlist,&
            ind,jnd,nnd,clist,ndim,verlet_skin,info)
    ENDIF
    IF (info .NE. 0) THEN
        print *,'cnl_vlist_build failed'
        info = -1
        GOTO 9999
    ENDIF

#if debug_verbosity > 1
    WRITE(cbuf,'(A,F8.2,2(A,I6),1X,L)')'Avg num of neighbors: ',&
        REAL(SUM(nvlist(1:Npart)),MK)/REAL(Npart,MK),&
        ' Max:',MAXVAL(nvlist(1:Npart)), ' Min:',MINVAL(nvlist(1:Npart))
    print *,cbuf

#endif

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
    DEALLOCATE(clist(1)%Nm)

    9999 CONTINUE ! jump here upon error

END SUBROUTINE cnl_vlist
