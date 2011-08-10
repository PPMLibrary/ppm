SUBROUTINE sop_spawn_particles(Particles,opts,info,nb_part_added,&
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

    INTEGER, OPTIONAL,                    INTENT(  OUT)   :: nb_part_added
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
    INTEGER, PARAMETER                     :: nb_new_part = 1
    REAL(MK)                               :: angle
#ifdef __USE_RANDOMNUMBERS
    LOGICAL                                :: alloc_rand
#endif
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
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &    'incompatible optional arguments',__LINE__,info)
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
    ! counting how many particles have to be added
    ! if not enough neighbours, add nb_new_part particle
    IF (opts%level_set) THEN
        IF (PRESENT(wp_fun)) THEN
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
            level => Get_wps(Particles,Particles%level_id)
            wp    => Get_wps(Particles,Particles%adapt_wpid)
            DO ip=1,Npart
                IF (nvlist(ip) .LT. nvlist_theoretical) THEN
                    IF (ABS(level(ip)).LT.opts%nb_width*nb_fun(wp(ip),&
                        opts%scale_D))&
                        add_part = add_part + nb_new_part
                ENDIF
            ENDDO
            level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
            wp    => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
        ENDIF
    ELSE
        !DO ip=1,Npart
            !IF (nvlist(ip) .LT. nvlist_theoretical) &
                !add_part = add_part + nb_new_part
        !ENDDO
        CALL check_quadrants(Particles,info)
        add_part = COUNT(Particles%nvlist .GT. 0) * nb_new_part
    ENDIF

    nvlist => NULL()

#if debug_verbosity > 1
        IF (PRESENT(printp)) THEN
            OPEN(6000+printp)
            OPEN(7000+printp)
        ENDIF
#endif
    !!-------------------------------------------------------------------------!
    !! Re-allocate (grow) arrays if necessary
    !!-------------------------------------------------------------------------!
    IF (add_part .GT. 0) THEN
        lda = (/ppm_dim,Npart+add_part/)
        CALL ppm_alloc(Particles%xp,lda,ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of xp failed',__LINE__,info)
            GOTO 9999
        ENDIF
        lda1 = Npart+add_part
        CALL ppm_alloc(Particles%wps(Particles%D_id)%vec,lda1,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of D failed',__LINE__,info)
            GOTO 9999
        ENDIF
        CALL ppm_alloc(Particles%wps(Particles%rcp_id)%vec,lda1,&
            ppm_param_alloc_grow_preserve,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &    'allocation of rcp failed',__LINE__,info)
            GOTO 9999
        ENDIF

#ifdef __USE_RANDOMNUMBERS
    !!-------------------------------------------------------------------------!
    !! Draw random numbers to add noise on positions of new particles
    !!-------------------------------------------------------------------------!
        !generate new random numbers if we are running out of them
        alloc_rand = .FALSE.
        IF (.NOT. ALLOCATED(randnb)) THEN
            alloc_rand = .TRUE.
        ELSE IF((ppm_dim-1)*randnb_i+nb_new_part*ppm_dim*(add_part+1) .GE. SIZE(randnb)) THEN
            alloc_rand = .TRUE.
            DEALLOCATE(randnb)
        ENDIF
        IF(alloc_rand) THEN
            ALLOCATE(randnb(nb_new_part*ppm_dim*(Npart+add_part)),STAT=info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,   &
                    &    'allocation of randnb failed',__LINE__,info)
                GOTO 9999
            ENDIF
            IF (.NOT.ASSOCIATED(ppm_particles_seed)) THEN
                CALL RANDOM_SEED(SIZE=ppm_particles_seedsize)
                ldc(1) = ppm_particles_seedsize
                CALL ppm_alloc(ppm_particles_seed,ldc(1:1),ppm_param_alloc_fit,info)
                IF (info .NE. 0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,   &
                        &            'allocation failed',__LINE__,info)
                    GOTO 9999
                ENDIF
                DO i=1,ppm_particles_seedsize
                    ppm_particles_seed(i)=ppm_rank*i*i*i*i
                ENDDO
                CALL RANDOM_SEED(PUT=ppm_particles_seed)
            ENDIF
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

    xp => Particles%xp !Cannot use get_xp, because we need access to elements
    ! of the array xp that are beyond Particles%Npart
    D => Particles%wps(Particles%D_id)%vec      !same reason
    rcp => Particles%wps(Particles%rcp_id)%vec  !same reason
    nvlist => Particles%nvlist
    IF (opts%level_set) THEN
        IF (.NOT. PRESENT(level_fun)) THEN
            level => Get_wps(Particles,Particles%level_id)
            wp    => Get_wps(Particles,Particles%adapt_wpid)
        ENDIF
    ENDIF

    IF (ppm_dim .EQ. 2) THEN
        add_particles2d: DO ip=1,Npart

            !IF (nvlist(ip) .LT. nvlist_theoretical) THEN
            IF (Particles%nvlist(ip).GT.0) THEN
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
                    angle = 0._MK
#ifdef __USE_RANDOMNUMBERS
                    randnb_i = randnb_i + 1
                    angle = -PI/4._mk*(randnb(randnb_i)+0.5_mk) 
#endif
                    angle = angle + PI/2._mk * &
                        REAL(Particles%nvlist(ip),MK)

                    xp(1:ppm_dim,Npart + add_part) = xp(1:ppm_dim,ip) + &
                        0.723_MK*D(ip)&       !radius
                        * (/COS(angle),SIN(angle)/) !direction 

                    D(Npart + add_part)   = D(ip)
                    rcp(Npart + add_part) = rcp(ip)
#if debug_verbosity > 1
                    IF (PRESENT(printp)) THEN
                        write(6000+printp,*) xp(1:ppm_dim,Npart+add_part)
                        write(7000+printp,'(4(E12.4,2X))') xp(1:ppm_dim,ip),&
                            xp(1:ppm_dim,ip)- xp(1:ppm_dim,Npart+add_part)
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
        preserve_wpv= (/ (i, i=1,0) /)) !F90-friendly way to init an empty array
        !preserve_wpv=(/ INTEGER :: /)) !valid only in F2003
    IF (info .NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'particles_updated_nb_part failed',__LINE__,info)
        GOTO 9999
    ENDIF

ENDIF !add_part .NE.0

#if debug_verbosity > 1
IF (PRESENT(printp)) THEN
    CLOSE(6000+printp)
    CLOSE(7000+printp)
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
IF (PRESENT(nb_part_added)) THEN
    nb_part_added = add_part
ENDIF

#if debug_verbosity > 0
CALL substop(caller,t0,info)
#endif
9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_spawn_particles

SUBROUTINE check_quadrants(Particles,info)
    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    INTEGER,                              INTENT(  OUT)   :: info
    !local variables

    REAL(MK),DIMENSION(:,:), POINTER                      :: xp=>NULL()
    INTEGER,DIMENSION(:),    POINTER                      :: nvlist=>NULL()
    INTEGER,DIMENSION(:,:),  POINTER                      :: vlist=>NULL()
    INTEGER                                               :: ip,iq,ineigh
    LOGICAL,DIMENSION(4)                                  :: needs_neigh_l

    info = 0
    nvlist => Particles%nvlist
    vlist => Particles%vlist

    xp => Get_xp(Particles,with_ghosts=.TRUE.)


    DO ip=1,Particles%Npart
        ineigh = 1
        needs_neigh_l = .TRUE.
        IF (nvlist(ip).GT.0) THEN
            DO WHILE(ANY(needs_neigh_l) .AND. ineigh.LE.nvlist(ip))
                iq = vlist(ineigh,ip)
                IF (xp(1,iq) .GT. xp(1,ip)) THEN
                    IF (xp(2,iq) .GT. xp(2,ip)) THEN
                        needs_neigh_l(1) = .FALSE.
                    ELSE
                        needs_neigh_l(4) = .FALSE.
                    ENDIF
                ELSE
                    IF (xp(2,iq) .GT. xp(2,ip)) THEN
                        needs_neigh_l(2) = .FALSE.
                    ELSE
                        needs_neigh_l(3) = .FALSE.
                    ENDIF
                ENDIF
                ineigh = ineigh + 1
            ENDDO
        ENDIF

        IF (ANY(needs_neigh_l)) THEN
            loop_quadrants: DO iq=1,4
                IF(needs_neigh_l(iq)) THEN
                    nvlist(ip) = iq
                    EXIT loop_quadrants
                ENDIF
            ENDDO loop_quadrants
        ELSE
            nvlist(ip) = 0
        ENDIF

    ENDDO
    xp => Set_xp(Particles,read_only=.TRUE.)
    vlist => NULL()
    nvlist => NULL()

END SUBROUTINE check_quadrants


SUBROUTINE check_duplicates(Particles)
    IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    !local variables

    REAL(MK),DIMENSION(:,:), POINTER                      :: xp=>NULL()
    INTEGER,DIMENSION(:,:),  POINTER                      :: vlist=>NULL()
    INTEGER                                               :: ip,iq,ineigh
    INTEGER                                               :: info
    
    info = 0

    CALL particles_apply_bc(Particles,Particles%active_topoid,info)
    CALL particles_mapping_partial(Particles,Particles%active_topoid,info)
    CALL particles_mapping_ghosts(Particles,Particles%active_topoid,info)
    CALL particles_neighlists(Particles,Particles%active_topoid,info)

    vlist => Particles%vlist
    xp => Get_xp(Particles,with_ghosts=.TRUE.)

    DO ip=1,Particles%Npart
        DO ineigh = 1,Particles%nvlist(ip)
            iq = vlist(ineigh,ip)
            IF (SUM((xp(1:ppm_dim,iq) - xp(1:ppm_dim,ip))**2).LT.1E-10) THEN
                write(*,*) 'duplicate particles'
                write(*,*) 'ip = ',ip,' iq = ',iq
            stop
            ENDIF
        ENDDO
    ENDDO
    xp => Set_xp(Particles,read_only=.TRUE.)
    vlist => NULL()

END SUBROUTINE check_duplicates



#undef __KIND
