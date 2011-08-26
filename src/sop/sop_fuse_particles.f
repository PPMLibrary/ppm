!!!----------------------------------------------------------------------------!
!!!
!!! Delete particles that are two close from each other
!!!
!!! i.e. when the distance between two particles is less than
!!!       threshold
!!!
!!! Uses ppm_alloc to grow/shrink arrays 
!!!
!!!
!!! Remark: One can use MEAN(D(ip),D(iq)) or MIN(D(ip),D(iq)) as the 
!!!         reference scale for particle fusion.
!!!
!!!----------------------------------------------------------------------------!
SUBROUTINE sop_fuse_particles(Particles,opts,info,wp_fun,printp)

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
    INTEGER, OPTIONAL                                     :: printp
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

    END INTERFACE

    ! local variables
    INTEGER                                :: ip,iq,ineigh,i,di,j
    CHARACTER(LEN=256)                     :: filename,cbuf
    CHARACTER(LEN=256)                     :: caller='sop_fuse_particles'
    REAL(KIND(1.D0))                       :: t0
    REAL(MK)                               :: dist,dist1,dist2,lev
    INTEGER                                :: del_part
    REAL(MK)                               :: threshold
    REAL(MK),     DIMENSION(:,:),POINTER   :: xp => NULL()
    REAL(MK),     DIMENSION(:,:),  POINTER :: inv => NULL()
    REAL(MK),     DIMENSION(:,:),  POINTER :: D => NULL()
    INTEGER                                :: Npart
    INTEGER,      DIMENSION(:),  POINTER   :: nvlist => NULL()
    INTEGER,      DIMENSION(:,:),POINTER   :: vlist => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: level => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: wp => NULL()

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

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
    D  => Get_wpv(Particles,Particles%D_id,with_ghosts=.TRUE.)
    inv=> Get_wpv(Particles,Particles%G_id,with_ghosts=.TRUE.)


    !removme
    !DO ip=1,Particles%Npart
        !write(107,'(4(E14.5,2X),I5)') xp(1:ppm_dim,ip),D(ip),&
            !level_fun(xp(1:ppm_dim,ip)),nvlist(ip)
    !ENDDO
    !write(*,*) opts%scale_D
    !stop
    !removme

    !!-------------------------------------------------------------------------!
    !! Mark particles for deletion (by changing nvlist to 999)
    !!-------------------------------------------------------------------------!
    ! haeckic: mark the right particles for deletion
    particle_loop: DO ip=1,Npart

        IF (nvlist(ip) .GT. opts%nneigh_critical) THEN
            DO ineigh=1,nvlist(ip)
                iq=vlist(ineigh,ip)
                CALL particles_anisotropic_distance(Particles,ip,iq,dist1,info)
                CALL particles_anisotropic_distance(Particles,iq,ip,dist2,info)
                dist = MIN(dist1,dist2)
                !dist=SQRT(SUM((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2)) / MIN(D(ip),D(iq))
                IF (dist .LT. 1D-8) THEN
                    WRITE(*,*) 'dist = ',dist, ' 2 particles very close'
                    WRITE(*,*) xp(1:ppm_dim,ip)
                    WRITE(*,*) xp(1:ppm_dim,iq)
                    !WRITE(*,*) D(ip),D(iq)
                    WRITE(*,*) 'something is probably wrong...'
                    info = -1
                    GOTO 9999
                ENDIF


!                      IF (ip.eq.14 .and. iq .eq. 1761 ) THEN
!                         write(*,*) 'disti ',ip,dist1,dist2
!                         write(*,*) 'inv ',ip, Particles%wpv(Particles%G_id)%vec(:,ip)
!                         write(*,*) 'inv ',iq, Particles%wpv(Particles%G_id)%vec(:,iq)
!                      ENDIF
!                      IF (ip.eq.1761 .and. iq .eq. 14 ) THEN
!                         write(*,*) 'disti ',ip,dist1,dist2
!                      ENDIF
                     IF (ip.eq.881 .and. iq .eq. 880 ) THEN
                        write(*,*) 'disti ',ip,dist1,dist2
                     ENDIF

                     IF (ip.eq.881 .or. ip .eq. 880 ) THEN
                        write(*,*) ip,iq
                     ENDIF



                     IF (ip.eq.880 .and. iq .eq. 881 ) THEN
                        write(*,*) 'disti ',ip,dist1,dist2
                     ENDIF


                IF (dist .LT. threshold) THEN

                     ! haeckic: a hack to make the comparison deterministic
                     dist1 = particles_longer_axis(Particles,ip)
                     dist2 = particles_longer_axis(Particles,iq)

!                      IF (ip.eq.14 .and. iq .eq. 1761 ) THEN
!                         write(*,*) 'disti ',ip,dist1-dist2
!                      ENDIF
!                      IF (ip.eq.1761 .and. iq .eq. 14 ) THEN
!                         write(*,*) 'disti ',ip,dist1-dist2
!                      ENDIF

                     IF (ip.eq.15 .and. iq .eq. 99 ) THEN
                        write(*,*) 'disti ',dist1-dist2
                     ENDIF
                     IF (ip.eq.99 .and. iq .eq. 15 ) THEN
                        write(*,*) 'disti ',dist1-dist2
                     ENDIF

                    ! Delete the particle that have the larger D
                    ! haeckic: decide which one to delete!?!? now just the one with bigger longer axis
                    IF (dist1.GT.dist2) THEN
                        nvlist(ip)=999
                        CYCLE particle_loop
                    ELSE IF(dist1.EQ.dist2) THEN
                                         IF (ip.eq.14 .and. iq .eq. 1761 ) THEN
                        write(*,*) 'disti ',ip,dist1-dist2
                     ENDIF
                     IF (ip.eq.1761 .and. iq .eq. 14 ) THEN
                        write(*,*) 'disti ',ip,dist1-dist2
                     ENDIF
                        ! If equal longer axis, then delete the one that has the "smallest position"
                        ! (any total ordering between particles would do)
                        DO di = 1,ppm_dim
                            IF (xp(di,ip) .LT. xp(di,iq)) THEN
                                !mark particle for deletion
                                nvlist(ip) = 999
                                CYCLE particle_loop
                            ELSE IF (xp(di,ip) .GT. xp(di,iq)) THEN
                                ! do nothing, this particle is going to stay
                                ! and its neighbour is going to be deleted
                                CYCLE particle_loop
                            ENDIF
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ENDIF
    ENDDO particle_loop

#if debug_verbosity > 1
        IF (PRESENT(printp)) THEN
            OPEN(5000+printp)
        ENDIF
#endif
    !!-------------------------------------------------------------------------!
    !! Delete particles that have been marked
    !!-------------------------------------------------------------------------!

    del_part = 0

    del_part_loop: DO ip=Npart,1,-1

            IF (ip.eq.881 .or. ip .eq. 880) THEN
               write(*,*) 'pos',ip, xp(1:ppm_dim,ip)
            ENDIF

        IF (nvlist(ip) .EQ. 999 ) THEN

#if debug_verbosity > 1
            IF(PRESENT(printp))THEN
                WRITE(5000+printp,*) xp(1:ppm_dim,ip)
            ENDIF
#endif

            IF (ip.eq.881 .or. ip .eq. 880) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip.eq.843 .or. ip .eq. 842) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip.eq.791 .or. ip .eq. 786) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip.eq.710 .or. ip .eq. 715) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip.eq.595 .or. ip .eq. 605) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

          IF (ip.eq.451 .or. ip .eq. 436) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip.eq.216 .or. ip .eq. 243) THEN
               write(*,*) 'here',ip, Npart - del_part
               write(*,*) 'pos',ip, xp(1:ppm_dim,Npart-del_part)
            ENDIF

            IF (ip .EQ. Npart - del_part) THEN
                ! no need to copy anything
                ! this particle will be deleted when Npart is decreased
            ELSE    
                ! copying particles from the end of xp to the index that has
                ! to be removed
                ! haeckic: here a correct copying
                xp(1:ppm_dim,ip) = xp(1:ppm_dim,Npart-del_part)
                DO j=1,Particles%tensor_length
                     D(j,ip)   = D(j,Npart-del_part)
                     inv(j,ip) = inv(j,Npart-del_part)
                ENDDO
!                 D(ip) = D(Npart-del_part)
!                 inv(ip) = inv(Npart-del_part)
                nvlist(ip) = nvlist(Npart-del_part)
            ENDIF
            del_part = del_part+1
        ENDIF
    ENDDO del_part_loop

#if debug_verbosity > 1
        IF (PRESENT(printp)) THEN
            CLOSE(5000+printp)
        ENDIF
#endif
    !New number of particles, after deleting some
    write(*,*) ' Particles fused: ', del_part
    Particles%Npart = Npart - del_part
    xp => Set_xp(Particles)
    D  => Set_wpv(Particles,Particles%D_id)
    inv=> Set_wpv(Particles,Particles%G_id)

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 1
#ifdef __MPI
    CALL MPI_Allreduce(del_part,del_part,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#endif
    IF (ppm_rank .EQ.0) THEN
        WRITE(cbuf,'(A,I8,A)') 'Deleting ', del_part,' particles'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF
#endif

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif
    9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_fuse_particles

#undef __KIND
