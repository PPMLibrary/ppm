!!!----------------------------------------------------------------------------!
!!!
!!! Delete particles that have too many close neighbours
!!!
!!!
!!! Uses ppm_alloc to grow/shrink arrays 
!!!
!!!
!!!
!!!----------------------------------------------------------------------------!
SUBROUTINE sop_fuse2_particles(Particles,opts,info,&
        level_fun,wp_fun,nb_fun,printp,nb_part_del)

    USE ppm_module_alloc, ONLY: ppm_alloc

    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    TYPE(sop_t_opts), POINTER,            INTENT(IN   )   :: opts
    INTEGER,                              INTENT(  OUT)   :: info

    OPTIONAL                                              :: wp_fun
    OPTIONAL                                              :: nb_fun
    OPTIONAL                                              :: level_fun
    !!! if level function is known analytically
    INTEGER, OPTIONAL                                     :: printp
    INTEGER, OPTIONAL                                     :: nb_part_del
    !!! printout particles that are deleted into file fort.(5000+printp)
    ! argument-functions need an interface
    INTERFACE
        FUNCTION wp_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)        :: pos
            REAL(MK)                                      :: wp_fun
        END FUNCTION wp_fun

        !Level function (usually known only during initialisation)
        FUNCTION level_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)        :: pos
            REAL(MK)                                      :: level_fun
        END FUNCTION level_fun

        !Function that returns the width of the narrow band
        FUNCTION nb_fun(kappa,scale_D)
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK)                             :: nb_fun
            REAL(MK),                INTENT(IN)  :: kappa
            REAL(MK),                INTENT(IN)  :: scale_D
        END FUNCTION nb_fun

    END INTERFACE

    ! local variables
    INTEGER                                :: ip,iq,ineigh,i,di
    CHARACTER(LEN=256)                     :: filename,cbuf
    CHARACTER(LEN=256)                     :: caller='sop_fuse2_particles'
    REAL(KIND(1.D0))                       :: t0
    REAL(MK)                               :: dist,lev
    INTEGER                                :: del_part
    REAL(MK)                               :: threshold
    REAL(MK),     DIMENSION(:,:),POINTER   :: xp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: rcp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: D => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: Dtilde => NULL()
    INTEGER                                :: Npart, close_neigh
    INTEGER,      DIMENSION(:),  POINTER   :: nvlist => NULL()
    INTEGER,      DIMENSION(:,:),POINTER   :: vlist => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: level => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: wp => NULL()

    INTEGER,DIMENSION(:),POINTER           :: fuse_part
    INTEGER,DIMENSION(:),POINTER           :: nb_neigh
    LOGICAL                                :: readonly
    INTEGER                                :: nb_neigh_max_thresh
    REAL(MK)                               :: rr

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
    IF (opts%del_parts) THEN
        IF (ppm_dim .eq. 2) THEN
            nb_neigh_max_thresh = 10
        ELSE
            nb_neigh_max_thresh = 15
        ENDIF
    ELSE
        nb_neigh_max_thresh = HUGE(1)
    ENDIF

#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    IF (opts%level_set) THEN
        IF (PRESENT(level_fun) .NEQV. PRESENT(wp_fun)) THEN
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'incompatible optional arguments',__LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF
    IF (.NOT. Particles%neighlists) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,&
            'Need neighbour lists: please compute them first',__LINE__,info)
        GOTO 9999
    ENDIF

    threshold = opts%fuse_radius
    
    nvlist => Particles%nvlist
    vlist => Particles%vlist
    Npart = Particles%Npart

    xp => Get_xp(Particles,with_ghosts=.TRUE.)
    D  => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
    Dtilde  => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
    rcp=> Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
    IF (opts%level_set) THEN
        IF(.NOT.PRESENT(level_fun)) &
            level => Get_wps(Particles,Particles%level_id,with_ghosts=.TRUE.)
        IF(.NOT.PRESENT(wp_fun)) &
            wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.TRUE.)
    ENDIF
    fuse_part  => Get_wpi(Particles,fuse_id,with_ghosts=.true.)
    nb_neigh  => Get_wpi(Particles,nb_neigh_id,with_ghosts=.true.)

    !!-------------------------------------------------------------------------!
    !! Mark particles for deletion (by changing nvlist to 999)
    !!-------------------------------------------------------------------------!
    particle_loop: DO ip=1,Npart

        IF (opts%level_set) THEN
            !kill particles that are too far away from the interface
            IF(PRESENT(level_fun)) THEN
                lev = level_fun(xp(1:ppm_dim,ip))
                IF (ABS(lev) .GT.  opts%nb_width_kill*&
                    nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) THEN
                    nvlist(ip)=999
                    CYCLE particle_loop
                ENDIF
            ELSE
                IF (ABS(level(ip)).GT.&
                    & opts%nb_width_kill*nb_fun(wp(ip),opts%scale_D)) THEN
                    nvlist(ip)=999
                    CYCLE particle_loop
                ENDIF
            ENDIF
        ENDIF
        IF (opts%remove_large_parts) THEN
            IF (D(ip).GE.MIN(opts%maximum_D,opts%scale_D)) THEN
                nvlist(ip)=999
                CYCLE particle_loop
            ENDIF
        ENDIF



        IF (nb_neigh(ip) .GT. nb_neigh_max_thresh) THEN
            close_neigh = 0
            check_close_neigh_loop: DO ineigh = 1,nb_neigh(ip)
                iq = vlist(ineigh,ip)

                !compute nearest-neighbour distance for each particle
                ! rescaled by the MIN of D(ip) and D(iq)
                ! (the definition of close_neigh in sop_close_neigh
                ! is slightly different to what is needed here and
                ! we have to update it again now)
                rr = SQRT(SUM((xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq))**2)) / &
                    MIN(D(ip),D(iq))

                IF (rr .LT. 1.0_mk) THEN
                    close_neigh = close_neigh + 1
                    !sort vlist (only for the close neighbours)
                    vlist(ineigh,ip) = vlist(close_neigh,ip)
                    vlist(close_neigh,ip) = iq 
                ENDIF
                nb_neigh(ip) = close_neigh

            ENDDO check_close_neigh_loop
        ENDIF
    ENDDO particle_loop


    particle_loop2: DO ip=1,Npart
        IF (nvlist(ip) .EQ. 999) CYCLE particle_loop2
        IF (nb_neigh(ip) .GT. nb_neigh_max_thresh) THEN
            DO ineigh = 1,nb_neigh(ip)
                iq=vlist(ineigh,ip)
                IF (iq.LE.Npart) THEN
                    IF (nb_neigh(iq).GT.0) nb_neigh(iq) = nb_neigh(iq)-1
                ENDIF
            ENDDO
            nvlist(ip) = 999
            cycle particle_loop2
        ENDIF
        neighbour_loop2: DO ineigh=1,nvlist(ip)
            iq=vlist(ineigh,ip)
            dist=SQRT(SUM((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2)) / &
                MIN(D(ip),D(iq))
           IF (dist .LT. threshold) THEN
                ! Delete the particle that has the larger D
                IF (D(ip).GT.D(iq)) THEN
                    nvlist(ip)=999
                    CYCLE particle_loop2
                ELSE IF(D(ip).EQ.(D(iq))) THEN

                    ! If equal D, then delete the one that has the "smallest position"
                    ! (any total ordering between particles would do)
                    DO di = 1,ppm_dim
                        IF (xp(di,ip) .LT. xp(di,iq)) THEN
                            !mark particle for deletion
                            nvlist(ip) = 999
                            CYCLE particle_loop2
                        ELSE IF (xp(di,ip) .EQ. xp(di,iq)) THEN
                            IF (ip.LT.iq) THEN
                                nvlist(ip) = 999
                                CYCLE particle_loop2
                            ENDIF
                        ELSE IF (xp(di,ip) .GT. xp(di,iq)) THEN
                            ! do nothing, this particle is going to stay
                            !  and its neighbour is going to be deleted
                             !CYCLE neighbour_loop
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
        ENDDO neighbour_loop2
    ENDDO particle_loop2

#if debug_verbosity > 2
        IF (PRESENT(printp)) THEN
            OPEN(5000+printp)
        ENDIF
#endif
    !!-------------------------------------------------------------------------!
    !! Delete particles that have been marked
    !!-------------------------------------------------------------------------!

    del_part = 0

    del_part_loop: DO ip=Npart,1,-1
        IF (nvlist(ip) .EQ. 999 ) THEN

#if debug_verbosity > 2
            IF(PRESENT(printp))THEN
                WRITE(5000+printp,*) xp(1:ppm_dim,ip)
            ENDIF
#endif

            IF (ip .EQ. Npart - del_part) THEN
                ! no need to copy anything
                ! this particle will be deleted when Npart is decreased
            ELSE    
                ! copying particles from the end of xp to the index that has
                ! to be removed
                xp(1:ppm_dim,ip) = xp(1:ppm_dim,Npart-del_part)
                D(ip) = D(Npart-del_part)
                Dtilde(ip) = Dtilde(Npart-del_part)
                rcp(ip) = rcp(Npart-del_part)
                nvlist(ip) = nvlist(Npart-del_part)
                fuse_part(ip) = 0
                nb_neigh(ip) = nb_neigh(Npart-del_part)
                IF (opts%level_set .AND. .NOT.PRESENT(level_fun)) THEN
                    level(ip) = level(Npart-del_part)
                    wp(ip) = wp(Npart-del_part)
                ENDIF
            ENDIF
            del_part = del_part+1
        ENDIF
    ENDDO del_part_loop

#if debug_verbosity > 2
        IF (PRESENT(printp)) THEN
            CLOSE(5000+printp)
        ENDIF
#endif
    !New number of particles, after deleting some

    
    Particles%Npart = Npart - del_part

#ifdef __MPI
    CALL MPI_Allreduce(del_part,del_part,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#endif

    readonly = (del_part.EQ.0) 

    xp => Set_xp(Particles,read_only=readonly)
    D  => Set_wps(Particles,Particles%D_id,read_only=readonly)
    Dtilde  => Set_wps(Particles,Particles%Dtilde_id,read_only=readonly)
    fuse_part => Set_wpi(Particles,fuse_id,read_only=readonly)
    nb_neigh  => Set_wpi(Particles,nb_neigh_id,read_only=readonly)
    rcp=> Set_wps(Particles,Particles%rcp_id,read_only=readonly)
    IF (opts%level_set) THEN
        IF(.NOT.PRESENT(level_fun)) &
            level=> Set_wps(Particles,Particles%level_id,read_only=readonly)
        IF(.NOT.PRESENT(wp_fun)) &
            wp => Set_wps(Particles,Particles%adapt_wpid,read_only=readonly)
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 1
    IF (ppm_rank .EQ.0) THEN
        WRITE(cbuf,'(A,I0,A,I0)') 'Deleting ', del_part,&
            ' particles - Np: ',Particles%Npart
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF
#endif
    IF(PRESENT(nb_part_del)) nb_part_del = del_part

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif
    9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_fuse2_particles

#undef __KIND
