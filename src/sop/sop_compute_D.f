!!!----------------------------------------------------------------------------!
!!! Computes resolution field D
!!!----------------------------------------------------------------------------!

SUBROUTINE sop_compute_D(Particles,D_fun,opts,info,     &
        wp_fun,wp_grad_fun,level_fun,level_grad_fun,nb_fun,&
        only_D_tilde,stats)

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
    OPTIONAL                                              :: nb_fun
    !!! if narrow-band function is known analytically
    LOGICAL, OPTIONAL                                     :: only_D_tilde
    !!! only compute D_tilde, then exits (no ghosts, no neighlists, no D)
    TYPE(sop_t_stats),  POINTER,OPTIONAL,  INTENT(  OUT)  :: stats
    !!! statistics on output

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
    CHARACTER(LEN = 256)                       :: caller='sop_compute_D'
    REAL(KIND(1.D0))                           :: t0

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: D_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: rcp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level_old => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad_old => NULL()

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
    REAL(MK),     DIMENSION(ppm_dim)           :: dummy_grad
    INTEGER                                    :: topo_id,eta_id
    REAL(MK),     DIMENSION(ppm_dim)           :: coeffs
    INTEGER,      DIMENSION(ppm_dim)           :: order
    INTEGER,      DIMENSION(ppm_dim*ppm_dim)   :: degree
    REAL(MK),     DIMENSION(ppm_dim)           :: wp_grad_fun0
    REAL(MK)                                   :: alpha

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
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'provided analytical gradients but not analytical &
            &   function values. This case is not yet implemented',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
        IF (Particles%level_set) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
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
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'Particles structure had not been defined. Call allocate first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check that we are dealing with adaptive particles 
    IF (.NOT.Particles%adaptive) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'These particles have not been declared as adaptive',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check all particles are inside the computational domain
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (opts%level_set) THEN
        IF (.NOT. Particles%level_set) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Need to enable level-set for Particles first',&
            &  __LINE__,info)
        GOTO 9999
        ENDIF
        IF (PRESENT(wp_fun) .AND..NOT.PRESENT(level_fun)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Need to provide analytical level function',&
            &  __LINE__,info)
        GOTO 9999
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------!
    ! Determines whether we will need to approximate derivatives
    !-------------------------------------------------------------------------!
    need_derivatives=.FALSE.
    IF (.NOT.PRESENT(wp_grad_fun)) THEN
        IF (Particles%level_set .OR. opts%D_needs_gradients) & 
            need_derivatives=.TRUE.
    ENDIF

    ! Check that the scalar field on which particles are supposed to adapt
    ! has been defined or is provided by an analytical function
    IF (.NOT.PRESENT(wp_fun)) THEN
        !check if a scalar property has already been specified as the argument
        !for the resolution function (i.e. the particles will adapt to resolve
        !this property well)
        IF (Particles%adapt_wpid.EQ.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
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
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,   &
                    &  'need to analytical gradients for level function',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ELSE
            !check if a scalar property has already been specified as the argument
            !for the resolution function (i.e. the particles will adapt to resolve
            !this property well)
            IF (Particles%level_id.EQ.0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,   &
                    &  'need to define level_id first',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%level_grad_id .EQ. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_argument,caller,   &
                    &  'need to define level_grad_id first',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF
    ENDIF


    !if the resolution depends on the gradient of wp, determines
    ! where this gradient is allocated
    IF (opts%D_needs_gradients) THEN
        !if so, checks whether we need to compute this gradient
        !and have an array allocated for it
        IF (PRESENT(wp_grad_fun)) THEN
            !no need to allocate an array for wp_grad
        ELSE
            IF (adapt_wpgradid.EQ.0) THEN
                !no array has already been specified for wp_grad
                !need to allocate one
                CALL particles_allocate_wpv(Particles,adapt_wpgradid,&
                    ppm_dim,info,with_ghosts=.FALSE.,name='adapt_wpgrad')
                IF (info.NE.0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'particles_allocate_wpv failed', __LINE__,info)
                    GOTO 9999
                ENDIF
            ELSE
                IF (.NOT.Particles%wpv(adapt_wpgradid)%is_mapped) THEN
                    CALL particles_allocate_wpv(Particles,adapt_wpgradid,&
                        ppm_dim,info,with_ghosts=.FALSE.,&
                        iopt=ppm_param_alloc_grow,name='adapt_wpgrad')
                    IF (info.NE.0) THEN
                        info = ppm_error_error
                        CALL ppm_error(ppm_err_alloc,caller,&
                            'particles_allocate_wpv failed',__LINE__,info)
                        GOTO 9999
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    ! crash if not enough neighbours
    IF (need_derivatives .AND. Particles%nneighmin.LT.opts%nneigh_critical) THEN
        xp => Get_xp(Particles)
        rcp => Get_wps(Particles,Particles%rcp_id)
        WRITE(cbuf,*) 'Not enough neighbours'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,'(2(A,I5,2X))') 'nneigh_critical ',opts%nneigh_critical,&
            'nneigh_toobig ',opts%nneigh_toobig
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,'(A,I5)') 'We have nneighmin = ',Particles%nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,*) 'Writing debug data to file fort.1230+rank'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(myformat,'(A,I1,A)') '(',ppm_dim+1,'(E30.16,2X),I4)'
        DO ip=1,Particles%Npart
            WRITE(1230+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
        ENDDO
        DO ip=Particles%Npart+1,Particles%Mpart
            WRITE(1250+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),0
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
        WRITE(cbuf,*) 'Writing debug data to file fort.1240+rank'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        DO ip=1,Particles%Npart
            WRITE(1240+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
        ENDDO
        DO ip=Particles%Npart+1,Particles%Mpart
            WRITE(1260+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),0
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
    ! (re)allocate Dtilde
    CALL particles_allocate_wps(Particles,Particles%Dtilde_id,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_grow,name='D_tilde')
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'particles_allocate_wps failed', __LINE__,info)
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Case where we need to approximate derivatives
    !!-------------------------------------------------------------------------!
    if_needs_derivatives: IF (need_derivatives) THEN
        IF (.NOT. Particles%neighlists) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'need neighbour lists to be uptodate', __LINE__,info)
            GOTO 9999
        ENDIF
        ! Compute gradients using DC operators
        coeffs=1._MK; order=3; degree = 0
        FORALL(i=1:ppm_dim) degree((i-1)*ppm_dim+i)=1 !Gradient
        eta_id = 0
        CALL particles_dcop_define(Particles,eta_id,coeffs,degree,&
            order,ppm_dim,info,name="gradient",vector=.TRUE.)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_define failed', __LINE__,info)
            GOTO 9999
        ENDIF
        IF (opts%check_dcops .AND. PRESENT(stats)) THEN
            CALL particles_dcop_compute(Particles,eta_id,info,&
                c=opts%c,min_sv=stats%min_sv)
            WRITE(cbuf,*) 'Smallest singular value for gradient was ',stats%min_sv
            CALL ppm_write(ppm_rank,caller,cbuf,info)
        ELSE
            CALL particles_dcop_compute(Particles,eta_id,info,c=opts%c)
        ENDIF
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_compute failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF if_needs_derivatives

    !-------------------------------------------------------------------------!
    ! Compute D_tilde
    !-------------------------------------------------------------------------!

    if_D_needs_grad: IF (opts%D_needs_gradients) THEN
        IF (opts%level_set) THEN
            CALL ppm_write(ppm_rank,caller,&
                'D_needs_gradients with level_sets: not supported',info)
            info = -1
            GOTO 9999
        ENDIF

        IF (.NOT.PRESENT(wp_grad_fun)) THEN
            IF (.NOT. Particles%neighlists) THEN
                CALL ppm_write(ppm_rank,caller,&
                    'trying to use neighb lists that dont exist',info)
                info = -1
                GOTO 9999
            ENDIF
            xp=>Get_xp(Particles,with_ghosts=.TRUE.)
            wp_grad => Get_wpv(Particles,adapt_wpgradid)
            eta => Get_dcop(Particles,eta_id)
            wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.TRUE.)

            DO ip = 1,Particles%Npart 
                wp_grad(1:ppm_dim,ip) = 0._MK
            ENDDO
            DO ip=1,Particles%Npart
                DO ineigh=1,Particles%nvlist(ip)
                    iq=Particles%vlist(ineigh,ip)
                    wp_grad(1:ppm_dim,ip) = wp_grad(1:ppm_dim,ip)+   &
                        (wp(iq)-wp(ip))*&
                        eta(1+(ineigh-1)*ppm_dim:ineigh*ppm_dim,ip)
                ENDDO
            ENDDO
            wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
            eta => Set_dcop(Particles,eta_id)
            wp_grad => Set_wpv(Particles,adapt_wpgradid)
            xp => Set_xp(Particles,read_only=.TRUE.)

            CALL particles_dcop_free(Particles,eta_id,info)
            IF (info.ne.0) THEN
                info = ppm_error_error
                CALL ppm_error(999,caller,&
                    'D_needs_gradients with level_sets: not supported',&
                    __LINE__,info)
                GOTO 9999
            ENDIF
        ENDIF

        !Compute Dtilde on real particles

        xp => Get_xp(Particles)
        Dtilde => Get_wps(Particles,Particles%Dtilde_id)

        IF (PRESENT(wp_fun)) THEN
            DO ip=1,Particles%Npart
                wp_grad_fun0= wp_grad_fun(xp(1:ppm_dim,ip))
                Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun0,opts) 
            ENDDO
        ELSE
            wp_grad => Get_wpv(Particles,adapt_wpgradid)
            wp => Get_wps(Particles,Particles%adapt_wpid)
            DO ip=1,Particles%Npart
                Dtilde(ip) = D_fun(wp(ip),wp_grad(1:ppm_dim,ip),opts)
            ENDDO
            wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
            wp_grad => Set_wpv(Particles,adapt_wpgradid,read_only=.TRUE.)
        ENDIF

        Dtilde => Set_wps(Particles,Particles%Dtilde_id)
        xp => Set_xp(Particles,read_only=.TRUE.)

        ! Get ghosts for D_tilde
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_mapping_ghosts failed',info)
            info = -1
            GOTO 9999
        ENDIF

    ELSE ! .NOT. D_needs_grad

        Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
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
                level => Get_wps(Particles,Particles%level_id,&
                    with_ghosts=.TRUE.)
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

    IF (PRESENT(only_D_tilde)) THEN
        IF (only_D_tilde) GOTO 8000
    ENDIF
    !---------------------------------------------------------------------!
    ! Increase the desired resolution where it is needed
    ! D^(n+1) = Min(D^(n),D_tilde^(n+1))
    !---------------------------------------------------------------------!
    !IF (Particles%D_id .EQ. 0 ) THEN
        !CALL particles_allocate_wps(Particles,Particles%D_id,&
            !info,name='D')
        !IF (info .NE. 0) THEN
            !info = ppm_error_error
            !CALL ppm_error(ppm_err_alloc,caller,&
                !'particles_allocate_wps failed',__LINE__,info)
            !GOTO 9999
        !ENDIF
        !D => Get_wps(Particles,Particles%D_id)
        !Dtilde => Get_wps(Particles,Particles%Dtilde_id)
        !DO ip=1,Particles%Npart
            !D(ip) = Dtilde(ip)
        !ENDDO
    !ELSE
        !D => Get_wps(Particles,Particles%D_id)
        !Dtilde => Get_wps(Particles,Particles%Dtilde_id)
        !DO ip=1,Particles%Npart
            !D(ip) = MIN(D(ip),Dtilde(ip))
        !ENDDO
    !ENDIF

    !D => Set_wps(Particles,Particles%D_id)
    !Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)

    !---------------------------------------------------------------------!
    ! Update cutoff radii
    !---------------------------------------------------------------------!
    rcp => Get_wps(Particles,Particles%rcp_id)
    Dtilde => Get_wps(Particles,Particles%Dtilde_id)
    DO ip=1,Particles%Npart
        !rcp(ip) = opts%rcp_over_D * Dtilde(ip)
        !TESTING THIS:
        rcp(ip) = Dtilde(ip)
    ENDDO
    rcp => Set_wps(Particles,Particles%rcp_id)
    Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)

    CALL particles_updated_cutoff(Particles,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'particles_updated_cutoff failed',__LINE__,info)
        GOTO 9999
    ENDIF

    !---------------------------------------------------------------------!
    ! Update ghosts
    !---------------------------------------------------------------------!
    CALL particles_mapping_ghosts(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'particles_mapping_ghosts failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !---------------------------------------------------------------------!
    ! Update neighbour lists
    !---------------------------------------------------------------------!
    CALL particles_neighlists(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'particles_neighlists failed',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%D_id .EQ. 0 ) THEN
        CALL particles_allocate_wps(Particles,Particles%D_id,&
            info,name='D')
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'particles_allocate_wps failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    !---------------------------------------------------------------------!
    ! D^(n+1) = min(D_tilde^(n+1)(iq)) over all neighbours iq
    !---------------------------------------------------------------------!
    D      => Get_wps(Particles,Particles%D_id)
    Dtilde => Get_wps(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
    xp => get_xp(Particles,with_ghosts=.true.)
    DO ip=1,Particles%Npart
        D(ip)=Dtilde(ip)
        DO ineigh=1,Particles%nvlist(ip)
            iq=Particles%vlist(ineigh,ip)

        !either this ....
            IF (Dtilde(iq).GE.Dtilde(ip)) CYCLE
            alpha = (sqrt(sum((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2))/Dtilde(iq)-1._mk) / &
                (opts%rcp_over_D-1._mk)
            if(alpha.le.0) then
                D(ip)=MIN(D(ip),Dtilde(iq))
            else
                D(ip)=MIN(D(ip),sqrt(alpha)*Dtilde(ip)+(1._mk-sqrt(alpha))*Dtilde(iq))
            endif

        !.... or this:
            !IF (Dtilde(iq).LT.D(ip)) &
                !D(ip)=Dtilde(iq)
        ENDDO
    ENDDO
    xp => set_xp(Particles,read_only=.true.)
    D      => Set_wps(Particles,Particles%D_id)
    Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)

#if debug_verbosity > 1
    D => Get_wps(Particles,Particles%D_id)
#ifdef __MPI
    CALL MPI_Allreduce(MINVAL(D(1:Particles%Npart)),min_D,1,&
        ppm_mpi_kind,MPI_MIN,ppm_comm,info)
#else
    min_D =MINVAL(D(1:Particles%Npart))
#endif
    IF (ppm_rank .EQ.0) THEN
        WRITE(cbuf,'(A,E12.4)') 'Min D = ',min_D
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF
    D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
#endif

    8000 CONTINUE


#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error


END SUBROUTINE sop_compute_D

#undef __KIND
