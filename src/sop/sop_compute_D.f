!!!----------------------------------------------------------------------------!
!!! Computes resolution field D
!!!----------------------------------------------------------------------------!

SUBROUTINE sop_compute_D(Particles,D_fun,opts,info,     &
        wp_fun,wp_grad_fun,level_fun,level_grad_fun,&
        D_needs_gradients,nb_fun)

    USE ppm_module_error
    USE ppm_module_dcops

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
    !!! particles
    !REAL(MK),  DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)   :: eta
    TYPE(sop_t_opts),  POINTER,           INTENT(IN   )   :: opts
    !!! options
    INTEGER,                              INTENT(  OUT)   :: info

    !optional arguments

    OPTIONAL                                              :: wp_fun
    !!! if field is known analytically
    OPTIONAL                                              :: wp_grad_fun
    !!! if field gradients are known analytically
    OPTIONAL                                              :: level_fun
    !!! if level function is known analytically
    OPTIONAL                                              :: level_grad_fun
    !!! if level gradients are known analytically
    LOGICAL, OPTIONAL                                     :: D_needs_gradients
    !!! if resolution depends on function values only (not on its derivatives)

    ! argument-functions need an interface
    INTERFACE
        !Monitor function
        FUNCTION D_fun(f1,dfdx,opts,f2)
            USE ppm_module_sop_typedef
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK)                               :: D_fun
            REAL(MK),                   INTENT(IN) :: f1
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN) :: dfdx
            TYPE(sop_t_opts),POINTER,   INTENT(IN) :: opts
            REAL(MK),OPTIONAL,          INTENT(IN) :: f2
        END FUNCTION D_fun

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

        !Field function (usually known only during initialisation)
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

        !Gradient of the field func. (usually known only during initialisation)
        FUNCTION wp_grad_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim)                      :: wp_grad_fun
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
        END FUNCTION wp_grad_fun

        !Gradient of the level func. (usually known only during initialisation)
        FUNCTION level_grad_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim)                      :: level_grad_fun
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
        END FUNCTION level_grad_fun
    END INTERFACE

    ! local variables
    INTEGER                                    :: i,ip,ineigh,iq
    CHARACTER(LEN = 64)                        :: myformat
    CHARACTER(LEN = 256)                       :: cbuf
    CHARACTER(LEN = 256)    :: caller='sop_compute_D'
    REAL(KIND(1.D0))                           :: t0

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: D_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: rcp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level_old => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad_old => NULL()

    TYPE(ppm_t_particles), POINTER             :: Particles_old

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: rcp => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: D => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: Dtilde => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: wp_grad => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad => NULL()


    REAL(MK),     DIMENSION(:,:), POINTER      :: eta => NULL()
    REAL(MK)                                   :: min_D
    LOGICAL                                    :: need_derivatives
    LOGICAL                                    :: D_needs_grad
    REAL(MK),    DIMENSION(ppm_dim)            :: dummy_grad
    INTEGER                                    :: topo_id

    !-------------------------------------------------------------------------!
    ! Initialize
    !-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif
    dummy_grad=0._MK
    topo_id = Particles%active_topoid

    !-------------------------------------------------------------------------!
    ! Checks consistency of parameters
    !-------------------------------------------------------------------------!
    IF (PRESENT(wp_grad_fun) .AND. .NOT.PRESENT(wp_fun)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'provided analytical gradients but not analytical &
            &   function values. This case is not yet implemented',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
        IF (Particles%level_set) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,   &
                & 'provided analytical function values but no &
                & function gradients. With level sets, &
                & this case is not yet implemented',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    !-------------------------------------------------------------------------!
    ! Perform consistency checks
    !-------------------------------------------------------------------------!
    !check data structure exists
    IF (.NOT.ASSOCIATED(Particles).OR..NOT.ASSOCIATED(Particles%xp)) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check that we are dealing with adaptive particles 
    IF (.NOT.Particles%adaptive) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'these particles have not been specified as adaptive',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check all particles are inside the computational domain
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (opts%level_set) THEN
        IF (.NOT. Particles%level_set) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Need to enable level-set for Particles first',&
            &  __LINE__,info)
        GOTO 9999
        ENDIF
        IF (PRESENT(wp_fun) .AND..NOT.PRESENT(level_fun)) THEN
        info = ppm_error_fatal
        CALL ppm_error(999,caller,   &
            &  'Need to provide analytical level function',&
            &  __LINE__,info)
        GOTO 9999
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------!
    ! Checks if we need gradients to compute D
    !-------------------------------------------------------------------------!
    D_needs_grad = .FALSE.
    IF (PRESENT(D_needs_gradients)) THEN
        D_needs_grad = D_needs_gradients
    ENDIF

    !-------------------------------------------------------------------------!
    ! Determines whether we will need to approximate derivatives
    !-------------------------------------------------------------------------!
    need_derivatives=.FALSE.
    IF (.NOT.PRESENT(wp_grad_fun)) THEN
        IF (Particles%level_set .OR. D_needs_grad) & 
            need_derivatives=.TRUE.
    ENDIF

    ! Check that the scalar field on which particles are supposed to adapt
    ! has been defined or is provided by an analytical function
    IF (.NOT.PRESENT(wp_fun)) THEN
        !check if a scalar property has already been specified as the argument
        !for the resolution function (i.e. the particles will adapt to resolve
        !this property well)
        IF (Particles%adapt_wpid.EQ.0) THEN
            info = ppm_error_fatal
            CALL ppm_error(999,caller,   &
                &  'need to define adapt_wpid first',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    IF (opts%level_set) THEN
        ! Check that the level set has been defined or
        ! is provided by an analytical function
        IF (PRESENT(level_fun)) THEN
            IF (.NOT.PRESENT(level_grad_fun)) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'need to analytical gradients for level function',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            !check if a scalar property has already been specified as the argument
            !for the resolution function (i.e. the particles will adapt to resolve
            !this property well)
            IF (Particles%level_id.EQ.0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'need to define level_id first',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%level_grad_id .EQ. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(999,caller,   &
                    &  'need to define level_grad_id first',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF


    !if the resolution depends on the gradient of wp, determines
    ! where this gradient is allocated
    IF (D_needs_grad) THEN
        !if so, checks whether we need to compute this gradient
        !and have an array allocated for it
        IF (PRESENT(wp_grad_fun)) THEN
            !no need to allocate an array for wp_grad
        ELSE
            IF (Particles%adapt_wpgradid.EQ.0) THEN
                !no array has already been specified for wp_grad
                !need to allocate one
                CALL particles_allocate_wpv(Particles,Particles%adapt_wpgradid,&
                    ppm_dim,info,with_ghosts=.FALSE.)
                IF (info.NE.0) THEN
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ELSE
                IF (Particles%wpv_m(Particles%adapt_wpgradid).NE.1) THEN
                    CALL particles_allocate_wpv(Particles,Particles%adapt_wpgradid,&
                        ppm_dim,info,with_ghosts=.FALSE.,&
                        iopt=ppm_param_alloc_grow_preserve)
                    info = ppm_error_fatal
                    CALL ppm_error(999,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    ! crash if not enough neighbours
    IF (need_derivatives .AND. Particles%nneighmin.LT.opts%nneigh_critical) THEN
        xp => Get_xp(Particles)
        rcp => Get_wps(Particles,Particles%rcp_id)
        WRITE(cbuf,*) 'Verlet lists should be up-to-date on entry. ',&
            'Neighbours are needed to compute nn2'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,'(2(A,I5,2X))') 'nneigh_critical ',opts%nneigh_critical,&
            'nneigh_toobig ',opts%nneigh_toobig
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,'(A,I5)') 'nneighmin: ',Particles%nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,*) 'Writing debug data to file fort.123'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(myformat,'(I1,A)') ppm_dim+1,'(E30.16,2X),I4'
        DO ip=1,Particles%Npart
            WRITE(123,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
        ENDDO
        DO ip=Particles%Npart+1,Particles%Mpart
            WRITE(125,myformat) xp(1:ppm_dim,ip),rcp(ip),0
        ENDDO
        CALL ppm_write(ppm_rank,caller,&
            'Calling neighlists anyway, but crashing just after that',info)
        CALL particles_neighlists(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_neighlists failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        WRITE(cbuf,'(A,I5,1X,I5)') 'nneighmin is now: ',Particles%nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,*) 'Writing debug data to file fort.124'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        DO ip=1,Particles%Npart
            WRITE(124,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
        ENDDO
        DO ip=Particles%Npart+1,Particles%Mpart
            WRITE(126,myformat) xp(1:ppm_dim,ip),rcp(ip),0
        ENDDO
        xp=>NULL()
        rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
        WRITE(cbuf,'(2(A,I6),2(A,E30.20))') 'Npart=',&
            Particles%Npart,' Mpart=',Particles%Mpart,&
            ' cutoff=',Particles%cutoff,' skin=',Particles%skin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Compute D (desired resolution)
    !!-------------------------------------------------------------------------!
    ! re-activate Dtilde_id if it had already been used
    CALL particles_allocate_wps(Particles,Particles%Dtilde_id,info,&
        with_ghosts=.TRUE.,zero=.TRUE.,iopt=ppm_param_alloc_grow_preserve)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'allocation failed',info)
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Case where we need to approximate derivatives
    !!-------------------------------------------------------------------------!
    if_needs_derivatives: IF (need_derivatives) THEN
        IF (.NOT. Particles%neighlists) THEN
            CALL ppm_write(ppm_rank,caller,'trying to use neighb &
                & lists that dont exist',info)
            info = -1
            GOTO 9999
        ENDIF
        ! Compute gradients using PSE kernels
        ! NOTE: not the most efficient way of computing first order
        ! derivatives (if one does not use the same kernels to also 
        ! compute the Laplacian 
        ! FIXME
        CALL ppm_part_dcops(Particles,Particles%eta_id,opts%c,info,&
            islaplacian=.TRUE.)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_part_dcops failed',info)
            info = -1
            GOTO 9999
        ENDIF
        !rcp => Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
        !CALL sop_PSE_secondorder(Particles%xp,rcp,Particles%Npart,&
            !Particles%Mpart,Particles%nvlist,Particles%vlist, &
            !eta,c,Particles%nneighmin,Particles%nneighmax,info)
        !rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)

    ENDIF if_needs_derivatives

    !-------------------------------------------------------------------------!
    ! Compute D_tilde
    !-------------------------------------------------------------------------!

    if_D_needs_grad: IF (D_needs_grad) THEN
        IF (opts%level_set) THEN
            CALL ppm_write(ppm_rank,caller,'D_needs_grad &
                & with level_sets: not supported',info)
            info = -1
            GOTO 9999
        ENDIF

        IF (.NOT.PRESENT(wp_grad_fun)) THEN

            xp=>Get_xp(Particles,with_ghosts=.TRUE.)
            wp_grad => Get_wpv(Particles,Particles%adapt_wpgradid)
            DO ip = 1,Particles%Npart 
                wp_grad(1:ppm_dim,ip) = 0._MK
            ENDDO


            IF (.NOT. Particles%neighlists) THEN
                CALL ppm_write(ppm_rank,caller,'trying to use neighb &
                    & lists that dont exist',info)
                info = -1
                GOTO 9999
            ENDIF

            eta => Get_wpv(Particles,Particles%eta_id)
            IF (PRESENT(wp_fun)) THEN
                DO ip=1,Particles%Npart
                    DO ineigh=1,Particles%nvlist(ip)
                        iq=Particles%vlist(ineigh,ip)
                        wp_grad(1:ppm_dim,ip) = wp_grad(1:ppm_dim,ip)+   &
                            0.5_MK*(wp_fun(xp(1:ppm_dim,ip))+            &
                            &       wp_fun(xp(1:ppm_dim,iq)))  *         &
                            (xp(1:ppm_dim,iq)-xp(1:ppm_dim,ip))*eta(ineigh,ip)
                    ENDDO
                ENDDO
            ELSE 

                wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.TRUE.)
                DO ip=1,Particles%Npart
                    DO ineigh=1,Particles%nvlist(ip)
                        iq=Particles%vlist(ineigh,ip)
                        wp_grad(1:ppm_dim,ip) = wp_grad(1:ppm_dim,ip)+        &
                            0.5_MK*(wp(ip)+wp(iq))*                   &
                            (xp(1:ppm_dim,iq)-xp(1:ppm_dim,ip))*eta(ineigh,ip)
                    ENDDO

                ENDDO
#if debug_verbosity > 0
                DO ip=1,Particles%Npart
                    IF(ANY(ISNAN(wp_grad(1:ppm_dim,ip)))) THEN
                        rcp=>Particles%wps(Particles%rcp_id)%vec
                        WRITE(cbuf,*) 'wp_grad NAN ',ip
                        CALL ppm_write(ppm_rank,caller,cbuf,info)
                        WRITE(cbuf,*) wp_grad(1:ppm_dim,ip)
                        CALL ppm_write(ppm_rank,caller,cbuf,info)
                        DO iq=1,Particles%Npart
                            IF (PRESENT(wp_fun)) THEN
                                write(127,*) xp(1:ppm_dim,iq), rcp(iq), &
                                    wp_fun(xp(1:ppm_dim,iq))
                            ELSE
                                write(127,*) xp(1:ppm_dim,iq), rcp(iq), wp(iq)
                            ENDIF
                        ENDDO
                        DO iq=Particles%Npart+1,Particles%Mpart
                            IF (PRESENT(wp_fun)) THEN
                                write(128,*) xp(1:ppm_dim,iq), rcp(iq), &
                                    wp_fun(xp(1:ppm_dim,iq))
                            ELSE
                                write(128,*) xp(1:ppm_dim,iq), rcp(iq), wp(iq)
                            ENDIF
                        ENDDO
                        rcp=>NULL()
                        info = -1
                        GOTO 9999
                    ENDIF
                    IF(MAXVAL(ABS(wp_grad(1:ppm_dim,ip))) .GT. 1E20 ) THEN
                        WRITE(cbuf,*) 'wp_grad big ', wp_grad
                        CALL ppm_write(ppm_rank,caller,cbuf,info)
                        info = -1
                        GOTO 9999
                    ENDIF
                ENDDO
#endif
                wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
            ENDIF

            eta => Set_wpv(Particles,Particles%eta_id,read_only=.TRUE.)
            wp_grad => Set_wpv(Particles,Particles%adapt_wpgradid,read_only=.FALSE.)
            xp => Set_xp(Particles,read_only=.TRUE.)

        ELSE !IF (PRESENT(wp_grad_fun)) 
            Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
            xp => Get_xp(Particles,with_ghosts=.TRUE.)
            DO ip=1,Particles%Mpart
                Dtilde(ip) = D_fun(wp_fun(xp(:,ip)),wp_grad_fun(xp(:,ip)),opts)
            ENDDO
            Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.FALSE.)
            xp => Set_xp(Particles,read_only=.TRUE.)
        ENDIF

        !Compute D on real particles

        xp => Get_xp(Particles,with_ghosts=.FALSE.)
        Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.FALSE.)

        IF (PRESENT(wp_fun)) THEN
            DO ip=1,Particles%Npart
                Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),&
                    wp_grad_fun(xp(1:ppm_dim,ip)),opts)
            ENDDO
        ELSE
            wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.FALSE.)
            DO ip=1,Particles%Npart
                Dtilde(ip) = D_fun(wp(ip),wp_grad(:,ip),opts)
            ENDDO
            wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
        ENDIF

        Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.FALSE.)
        xp => Set_xp(Particles,read_only=.TRUE.)

        ! Get ghosts for D_tilde
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_mapping_ghosts failed',info)
            info = -1
            GOTO 9999
        ENDIF

    ELSE ! .NOT. D_needs_grad

        Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.FALSE.)
        IF (PRESENT(wp_fun)) THEN
            xp => Get_xp(Particles,with_ghosts=.TRUE.)
            IF (opts%level_set) THEN
                DO ip=1,Particles%Mpart
                    Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),dummy_grad,&
                        opts,level_fun(xp(1:ppm_dim,ip)))
                ENDDO
            ELSE
                DO ip=1,Particles%Mpart
                    Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),dummy_grad,opts)
                ENDDO
            ENDIF
            xp => Set_xp(Particles,read_only=.TRUE.)
        ELSE
            wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.TRUE.)
            IF (opts%level_set) THEN
                level => Get_wps(Particles,Particles%level_id,with_ghosts=.TRUE.)
                DO ip=1,Particles%Mpart
                    Dtilde(ip) = D_fun(wp(ip),dummy_grad,opts,level(ip))
                ENDDO
                level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
            ELSE
                DO ip=1,Particles%Mpart
                    Dtilde(ip) = D_fun(wp(ip),dummy_grad,opts)
                ENDDO
            ENDIF
            wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
        ENDIF
        Dtilde => Set_wps(Particles,Particles%Dtilde_id,&
            read_only=.FALSE.,ghosts_ok=.TRUE.)

    ENDIF if_D_needs_grad


    !-------------------------------------------------------------------------!
    ! Rescale D_tilde
    !-------------------------------------------------------------------------!
    Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
    DO ip = 1,Particles%Mpart
        Dtilde(ip) = MIN(MAX(Dtilde(ip),opts%minimum_D),opts%maximum_D)
    ENDDO
    Dtilde => Set_wps(Particles,Particles%Dtilde_id,&
        read_only=.FALSE.,ghosts_ok=.TRUE.)

    !---------------------------------------------------------------------!
    ! Increase the desired resolution where it is needed
    ! D^(n+1) = Min(D^(n),D_tilde^(n+1))
    !---------------------------------------------------------------------!
    IF (Particles%D_id .EQ. 0 ) THEN
        CALL particles_allocate_wps(Particles,Particles%D_id,&
            info,with_ghosts=.TRUE.)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_allocate_wps failed',info)
            info = -1
            GOTO 9999
        ENDIF
        D => Get_wps(Particles,Particles%D_id,with_ghosts=.FALSE.)
        Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
        DO ip=1,Particles%Mpart
            D(ip) = Dtilde(ip)
        ENDDO
    ELSE
        D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
        Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
        DO ip=1,Particles%Mpart
            D(ip) = MIN(D(ip),Dtilde(ip))
        ENDDO
    ENDIF

    !removme
    IF (ANY(ISNAN(D(1:Particles%Mpart)))) &
        WRITE(*,*) 'NAN in D'
    !removme

    D => Set_wps(Particles,Particles%D_id,read_only=.FALSE.,ghosts_ok=.TRUE.)
    Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)


    !---------------------------------------------------------------------!
    ! Update cutoff radii
    !---------------------------------------------------------------------!
    rcp => Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
    D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
    DO ip=1,Particles%Mpart
        rcp(ip) = opts%rcp_over_D * D(ip)
    ENDDO
    CALL particles_updated_cutoff(Particles,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'particles_updated_cutoff failed',info)
        info = -1
        GOTO 9999
    ENDIF
    rcp => Set_wps(Particles,Particles%rcp_id,read_only=.FALSE.)
    D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)

    !---------------------------------------------------------------------!
    ! Update ghosts
    !---------------------------------------------------------------------!
    CALL particles_mapping_ghosts(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'particles_mapping_ghosts failed.',info)
        info = -1
        GOTO 9999
    ENDIF

write(*,*) 'ok 2.6'
    !---------------------------------------------------------------------!
    ! Update neighbour lists
    !---------------------------------------------------------------------!
    CALL particles_neighlists(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'particles_neighlists failed.',info)
        info = -1
        GOTO 9999
    ENDIF

write(*,*) 'ok 2.7'
    !---------------------------------------------------------------------!
    ! D^(n+1) = min(D_tilde^(n+1)(iq)) over all neighbours iq
    !---------------------------------------------------------------------!
    D      => Get_wps(Particles,Particles%D_id,     with_ghosts=.FALSE.)
    Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
    DO ip=1,Particles%Npart
        D(ip)=Dtilde(ip)
        DO ineigh=1,Particles%nvlist(ip)
            iq=Particles%vlist(ineigh,ip)
            IF (Dtilde(iq).LT.D(ip)) &
                D(ip)=Dtilde(iq)
        ENDDO
    ENDDO
    D      => Set_wps(Particles,Particles%D_id,     read_only=.FALSE.)
    Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)

#if debug_verbosity > 0
    D => Get_wps(Particles,Particles%D_id,with_ghosts=.FALSE.)
    CALL MPI_Allreduce(MINVAL(D(1:Particles%Npart)),min_D,1,&
        ppm_mpi_kind,MPI_MIN,ppm_comm,info)
    IF (ppm_rank .EQ.0) THEN
        WRITE(cbuf,'(A,E12.4)') 'Min D = ',min_D
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF
    D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
#endif

write(*,*) 'ok 2.8'

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error


END SUBROUTINE sop_compute_D

#undef __KIND
