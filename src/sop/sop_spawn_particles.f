SUBROUTINE sop_spawn_particles(Particles,opts,info,&
        nneigh_threshold,level_fun,wp_fun,nb_fun,printp)
    !!!----------------------------------------------------------------------------!
    !!!
    !!! Insert particles around those with too few neighbours
    !!!
    !!! Uses ppm_alloc to grow/shrink arrays 
    !!!
    !!! Warning: on output, some particles may be outside the computational domain
    !!! Call impose_boundary_conditions to fix this.
    !!!----------------------------------------------------------------------------!

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

    INTEGER, OPTIONAL                                     :: nneigh_threshold
    OPTIONAL                                              :: wp_fun
    OPTIONAL                                              :: nb_fun
    OPTIONAL                                              :: level_fun
    !!! if level function is known analytically
    INTEGER, OPTIONAL                                     :: printp
    !!! printout particles that are created into file fort.(6000+printp)
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
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
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
    REAL(MK),     DIMENSION(:,:),POINTER   :: xp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: rcp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: D => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: wp => NULL()
    REAL(MK),     DIMENSION(:),  POINTER   :: level => NULL()
    INTEGER                                :: Npart
    INTEGER,      DIMENSION(:),  POINTER   :: nvlist => NULL()
    INTEGER                                :: ip,iq,ineigh,i,di
    CHARACTER(LEN=256)                     :: cbuf
    CHARACTER(LEN=256)                     :: caller='sop_spawn_particles'
    REAL(KIND(1.D0))                       :: t0
    REAL(MK)                               :: lev
    REAL(MK)                               :: theta1,theta2
    INTEGER                                :: nvlist_theoretical
    INTEGER                                :: add_part
    INTEGER,        DIMENSION(2)           :: lda
    INTEGER,        DIMENSION(1)           :: lda1
    INTEGER, PARAMETER                     :: nb_new_part = 3
    !!! number of new particles that are generated locally

    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    IF (opts%level_set) THEN
        IF (PRESENT(level_fun) .NEQV. PRESENT(wp_fun)) THEN
            CALL ppm_write(ppm_rank,caller,'incompatible optional arguments',info)
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    IF (PRESENT(nneigh_threshold)) THEN
        nvlist_theoretical = nneigh_threshold
    ELSE
        nvlist_theoretical = opts%nneigh_theo
    ENDIF

    add_part = 0

    nvlist => Particles%nvlist
    Npart = Particles%Npart

    !!-------------------------------------------------------------------------!
    !! Count number of new particles to insert
    !!-------------------------------------------------------------------------!
    IF (PRESENT(wp_fun)) THEN
        ! counting how many particles have to be added
        ! if not enough neighbours, add nb_new_part particle
        IF (opts%level_set) THEN
            xp => Get_xp(Particles)
            DO ip=1,Npart
                IF (nvlist(ip) .LT. nvlist_theoretical) THEN
                    lev = level_fun(xp(1:ppm_dim,ip))
                    IF (ABS(lev) .LT. opts%nb_width*&
                        nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) THEN
                        add_part = add_part + nb_new_part
                    ENDIF
                ENDIF
            ENDDO
            xp => Set_xp(Particles,read_only=.TRUE.)
        ELSE
            DO ip=1,Npart
                IF (nvlist(ip) .LT. nvlist_theoretical) &
                    add_part = add_part + nb_new_part
            ENDDO
        ENDIF
    ELSE
        IF (opts%level_set) THEN
            level => Get_wps(Particles,Particles%level_id)
            wp    => Get_wps(Particles,Particles%adapt_wpid)
            DO ip=1,Npart
                IF (nvlist(ip) .LT. nvlist_theoretical) THEN
                    IF (ABS(level(ip)).LT.opts%nb_width*nb_fun(wp(ip),opts%scale_D))&
                        add_part = add_part + nb_new_part
                ENDIF
            ENDDO
            level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
            wp    => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
        ELSE
            DO ip=1,Npart
                IF (nvlist(ip) .LT. nvlist_theoretical) &
                    add_part = add_part + nb_new_part
            ENDDO
        ENDIF
    ENDIF

    nvlist => NULL()

#if debug_verbosity > 1
        IF (PRESENT(printp)) THEN
            OPEN(6000+printp)
        ENDIF
#endif
    !!-------------------------------------------------------------------------!
    !! Re-allocate (grow) arrays if necessary
    !!-------------------------------------------------------------------------!
    IF (add_part .GT. 0) THEN
        lda = (/ppm_dim,Npart+add_part/)
        CALL ppm_alloc(Particles%xp,lda,ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'allocation failed',info)
            info = -1
            GOTO 9999
        ENDIF
        lda1 = Npart+add_part
        CALL ppm_alloc(Particles%wps(Particles%D_id)%vec,lda1,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'allocation failed',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_alloc(Particles%wps(Particles%rcp_id)%vec,lda1,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'allocation failed',info)
            info = -1
            GOTO 9999
        ENDIF

#ifdef __USE_RANDOMNUMBERS
    !!-------------------------------------------------------------------------!
    !! Draw random numbers to add noise on positions of new particles
    !!-------------------------------------------------------------------------!
        !generate new random numbers if we are running out of them
        IF((ppm_dim-1)*randnb_i+nb_new_part*ppm_dim*(add_part+1) .GE. SIZE(randnb)) THEN
            DEALLOCATE(randnb)
            ALLOCATE(randnb(nb_new_part*ppm_dim*(Npart+add_part)),STAT=info)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'allocation failed.',info)
                info = -1
                GOTO 9999
            ENDIF
            DO i=1,seedsize
                seed(i)=10+i*i*i + Npart
            ENDDO
            CALL RANDOM_SEED(PUT=seed)
            CALL RANDOM_NUMBER(randnb)
            randnb_i = 0
        ENDIF
#endif

    !!-------------------------------------------------------------------------!
    !! Add particles
    !! (either randomly, or using fake random numbers based on particles'
    !! positions. This is a quick hack to ensure exact consistency between
    !! simulations run an different number of processors)
    !!-------------------------------------------------------------------------!
    add_part = 0

    xp => Get_xp(Particles)
    D => Get_wps(Particles,Particles%D_id)
    rcp => Get_wps(Particles,Particles%rcp_id)
    nvlist => Particles%nvlist
    IF (opts%level_set) THEN
        IF (.NOT. PRESENT(level_fun)) THEN
            level => Get_wps(Particles,Particles%level_id)
            wp    => Get_wps(Particles,Particles%adapt_wpid)
        ENDIF
    ENDIF

    IF (ppm_dim .EQ. 2) THEN
        add_particles2d: DO ip=1,Npart

            IF (nvlist(ip) .LT. nvlist_theoretical) THEN
                !FOR LEVEL SETS ONLY
                IF (opts%level_set) THEN
                    IF (PRESENT(wp_fun)) THEN
                        lev = level_fun(xp(1:ppm_dim,ip))
                        IF (ABS(lev) .GE. opts%nb_width*&
                            nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) &
                            CYCLE add_particles2d
                    ELSE
                        IF (ABS(level(ip)).GE.opts%nb_width*nb_fun(wp(ip),opts%scale_D)) &
                            CYCLE add_particles2d
                    ENDIF
                ENDIF

                DO i=1,nb_new_part
                    add_part = add_part + 1

#ifdef __USE_RANDOMNUMBERS
                    randnb_i = randnb_i + 1
#endif

                    xp(1:ppm_dim,Npart + add_part) = xp(1:ppm_dim,ip) + &
                        0.7_MK*D(ip)&       !radius
#ifdef __USE_RANDOMNUMBERS
                        * (/COS(2.094_MK*i + randnb(randnb_i)),&
                            &   SIN(2.094_MK*i + randnb(randnb_i)) &
                            &   /) !direction (n 2pi/3) + random phase
#else
                        * (/COS(2.094_MK*i + 1000._MK*xp(1,ip)/D(ip)) , &
                            &   SIN(2.094_MK*i + 1000._MK*xp(1,ip)/D(ip)) &
                            &   /) !direction (n 2pi/3) + random phase
#endif

                    D(Npart + add_part)   = D(ip)
                    rcp(Npart + add_part) = rcp(ip)
#if debug_verbosity > 1
                    IF (PRESENT(printp)) THEN
                        write(6000+printp,*) xp(1:ppm_dim,Npart+add_part)
                    ENDIF
#endif
                ENDDO
            ENDIF
        ENDDO add_particles2d
    ELSE ! if ppm_dim .eq. 3
        add_particles3d: DO ip=1,Npart

            IF (nvlist(ip) .LT. nvlist_theoretical) THEN
                !FOR LEVEL SETS ONLY
                IF (Particles%level_set) THEN
                    IF (PRESENT(wp_fun)) THEN
                        lev = level_fun(xp(1:ppm_dim,ip))
                        IF (ABS(lev) .GE. opts%nb_width*&
                            nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) &
                            CYCLE add_particles3d
                    ELSE
                        IF (ABS(level(ip)).GE.opts%nb_width*nb_fun(wp(ip),opts%scale_D)) &
                            CYCLE add_particles3d
                    ENDIF
                ENDIF

                DO i=1,nb_new_part
                    add_part = add_part + 1

#ifdef __USE_RANDOMNUMBERS
                    randnb_i = randnb_i + 1
#endif

#ifdef __USE_RANDOMNUMBERS
                    theta1 = ACOS(1._MK - 2._MK*randnb(2*randnb_i))
                    theta2 = 2._MK * PI * randnb(2*randnb_i-1)
#else
                    theta1 = ACOS(SIN(1000._MK*xp(1,ip)/D(ip)))
                    theta2 = PI * (1._MK+COS(1000._MK*xp(2,ip)/D(ip)))
#endif
                    xp(1:ppm_dim,Npart + add_part) = xp(1:ppm_dim,ip) + &
                        !random 3D points on a sphere
                    0.7_MK*D(ip)&       !radius
                        * (/COS(theta1), &
                        &   SIN(theta1) * COS(theta2), &
                        &   SIN(theta1) * SIN(theta2)  &
                        &   /) 
                    D(Npart + add_part)   = D(ip)
                    rcp(Npart + add_part) = rcp(ip)
#if debug_verbosity > 1
                    IF (PRESENT(printp)) THEN
                        write(6000+printp,*) xp(1:ppm_dim,Npart+add_part)
                    ENDIF
#endif
                ENDDO
            ENDIF
        ENDDO add_particles3d
    ENDIF

    IF (opts%level_set) THEN
        IF (.NOT. PRESENT(level_fun)) THEN
            level => Set_wps(Particles,Particles%level_id)
            wp    => Set_wps(Particles,Particles%adapt_wpid)
        ENDIF
    ENDIF
    !!-------------------------------------------------------------------------!
    !!new size Npart
    !!-------------------------------------------------------------------------!

    xp => Set_xp(Particles)
    D => Set_wps(Particles,Particles%D_id)
    rcp => Set_wps(Particles,Particles%rcp_id)
    nvlist => NULL()
    Particles%Npart = Npart + add_part

    CALL particles_updated_nb_part(Particles,info,&
        preserve_wps=(/Particles%D_id,Particles%rcp_id/),&
        preserve_wpv=(/ INTEGER :: /))
    IF (info .NE.0) THEN
        CALL ppm_write(ppm_rank,caller,'particles_updated_nb_part failed',info)
        info = -1
        GOTO 9999
    ENDIF

ENDIF !add_part .NE.0

#if debug_verbosity > 1
IF (PRESENT(printp)) THEN
    CLOSE(6000+printp)
ENDIF
#endif

!!-------------------------------------------------------------------------!
!! Finalize
!!-------------------------------------------------------------------------!
#if debug_verbosity > 1
#ifdef __MPI
CALL MPI_Allreduce(add_part,add_part,1,MPI_INTEGER,MPI_SUM,ppm_comm,info)
#endif
IF (ppm_rank .EQ.0) THEN
    WRITE(cbuf,'(A,I8,A)') 'Adding ', add_part,' particles'
    CALL ppm_write(ppm_rank,caller,cbuf,info)
ENDIF
#endif

#if debug_verbosity > 0
CALL substop(caller,t0,info)
#endif
9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_spawn_particles

#undef __KIND
