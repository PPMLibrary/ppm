!!!----------------------------------------------------------------------------!
!!! This is the core routine of the adaptive particle method
!!!
!!! It adapts the particles to the current field wp (or to an function provided
!!! analytically) using particle-particle interpolation (except if wp_fun is)
!!! provided and makes the corresponding changes to the 
!!! properties carried by particles
!!!
!!! Requires:
!!!      up-to-date cutoff radii
!!!      up-to-date Verlet lists
!!!      up-to-date field values (incl. ghosts)
!!!
!!! Returns:
!!!      new particle positions
!!!      updated field values (interpolated from old ones, only REAL particles)
!!!      updated cutoff radii
!!!      updated Verlet lists
!!!      updated field derivatives (OPTIONAL)
!!!
!!! Optional arguments:
!!!      wp_fun: a function that returns the value of wp at any given position
!!!                   (if absent, the values are interpolated and 
!!!                    stored in the array wp)
!!!
!!!      wp_grad_fun: a function that returns the gradient of wp.
!!!                   (if absent, the gradients are approximated and 
!!!                    stored in wp_grad)
!!!                   (requires: wp_fun)
!!!
!!!      wp_lap: pointer to an array (not necessarily allocated on input)
!!!              On output, contains laplacian of wp computed using data
!!!              of wp_old
!!!
!!!
!!! Do the gradient descent on xp, using interpolation from {xp_old} to
!!! {xp} at each step (requires cross--neighbour-list)
!!!
!!! FIXME: optimize the number of (re)-allocation needed (D, D_old,
!!! xp_old,etc..)
!!! FIXME: optimize how many times the ghosts have to be reconstructed
!!!----------------------------------------------------------------------------!

SUBROUTINE sop_adapt_particles(topo_id,Particles,D_fun,opts,info,     &
        wp_fun,wp_grad_fun,wplap_id,return_grad,level_fun,level_grad_fun,&
        smallest_sv,D_needs_gradients,nb_fun,stats)

    USE ppm_module_error
    USE ppm_module_map_part

    IMPLICIT NONE
#ifdef __MPI
    INCLUDE 'mpif.h'
#endif

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER  :: MK = ppm_kind_double
#endif
    ! arguments
    INTEGER,                              INTENT(IN   )   :: topo_id
    !!! topology
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    !!! particles
    !REAL(MK),  DIMENSION(:,:),ALLOCATABLE,INTENT(  OUT)   :: eta
    TYPE(sop_t_opts),  POINTER,           INTENT(INOUT)   :: opts
    !!! options
    INTEGER,                              INTENT(  OUT)   :: info

    !optional arguments

    LOGICAL, OPTIONAL                                     :: return_grad
    !!! return field gradients on output
    INTEGER,                      OPTIONAL,INTENT(INOUT)  :: wplap_id
    !!! if field laplacian is needed, it can be computed cheaply
    OPTIONAL                                              :: wp_fun
    !!! if field is known analytically
    OPTIONAL                                              :: wp_grad_fun
    !!! if field gradients are known analytically
    OPTIONAL                                              :: level_fun
    !!! if level function is known analytically
    OPTIONAL                                              :: level_grad_fun
    !!! if level gradients are known analytically
    OPTIONAL                                              :: nb_fun
    !!! function that describes the narrow band
    REAL(MK),OPTIONAL,                    INTENT(INOUT)   :: smallest_sv
    !!! smallest singular value of the Vandermonde matrices (DC-PSE)
    LOGICAL, OPTIONAL                                     :: D_needs_gradients
    !!! if resolution depends on function values only (not on its derivatives)
    TYPE(sop_t_stats),  POINTER,OPTIONAL,  INTENT(  OUT)  :: stats
    !!! statistics on output

    ! argument-functions need an interface
    INTERFACE
        !Monitor function
        FUNCTION D_fun(f1,dfdx,opts,f2)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
            USE ppm_module_sop_typedef
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
    INTEGER                                    :: iunit
    CHARACTER(LEN = 256)                       :: filename,cbuf
    CHARACTER(LEN = 256)    :: caller='sop_adapt_particles'
    REAL(KIND(1.D0))                           :: t0
    REAL(MK)                                   :: dist2
    REAL(MK)                                   :: Psi_threshold
    REAL(MK)                                   :: dummy_real
    INTEGER                                    :: num_it

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


    REAL(MK),     DIMENSION(:),   POINTER      :: nn2 => NULL()
    INTEGER                                    :: nn2_id
    INTEGER,      DIMENSION(:),   POINTER      :: nvlist_cross => NULL()
    INTEGER,      DIMENSION(:,:), POINTER      :: vlist_cross => NULL()
    REAL(MK),     DIMENSION(:,:), ALLOCATABLE  :: eta_interp
    INTEGER                                    :: nneighmax_cross
    INTEGER                                    :: nneighmin_cross
    REAL(MK)                                   :: tempvar
    LOGICAL                                    :: need_derivatives
    LOGICAL                                    :: D_needs_grad
    REAL(MK),    DIMENSION(ppm_dim)            :: dummy_grad

    !-------------------------------------------------------------------------!
    ! Initialize
    !-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif
    dummy_grad=0._MK

    !-------------------------------------------------------------------------!
    ! Checks if we need gradients to compute D
    !-------------------------------------------------------------------------!
    D_needs_grad = .FALSE.
    IF (PRESENT(D_needs_gradients)) THEN
        D_needs_grad = D_needs_gradients
    ENDIF

    !-------------------------------------------------------------------------!
    ! Checks consistency of parameters
    !-------------------------------------------------------------------------!
    IF (topo_id .NE. Particles%active_topoid) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Particles need to be on the same topology as &
            &   the one refered by the topo_id argument',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (PRESENT(wp_grad_fun) .AND. .NOT.PRESENT(wp_fun)) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'provided analytical gradients but not analytical &
            &   function values. This case is not yet implemented',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
        IF (D_needs_grad) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,   &
                & 'provided analytical function values but no &
                & function gradients. This case is not implemented',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        IF (Particles%level_set) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,   &
                & 'provided analytical function values but no &
                & function gradients. With level sets, &
                & this case is not yet implemented',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    IF (PRESENT(wp_fun) .AND. PRESENT(wplap_id)) THEN
        !if the field is known analytically, then we should not need to compute
        ! its laplacian using interpolating kernels.
        WRITE(cbuf,*) 'Calling adapt_particles with conflicting options'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        WRITE(cbuf,*) 'Warning: on exit, wp_lap will have nonsense values'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
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
        CALL ppm_error(999,caller,   &
            &  'these particles have not been specified as adaptive',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check all particles are inside the computational domain
    IF (.NOT.Particles%areinside) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'some particles may be outside the domain. Apply BC first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !check that neighbour lists have been computed
    IF (.NOT.Particles%neighlists) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'neighbour lists are not up-to-date. Call neighlists first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT. ASSOCIATED(opts)) THEN
        !provide defaults for opts
        CALL sop_init_opts(opts,info)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &  'allocation failed for opts',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    IF (opts%level_set) THEN
        IF (.NOT. Particles%level_set) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,   &
            &  'Need to enable level-set for Particles first',&
            &  __LINE__,info)
        GOTO 9999
        ENDIF
        IF (PRESENT(level_fun) .NEQV. PRESENT(wp_fun)) THEN
            CALL ppm_write(ppm_rank,caller,'incompatible optional arguments',info)
            info = -1
            GOTO 9999
        ENDIF
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
            info = ppm_error_error
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
                info = ppm_error_error
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
                info = ppm_error_error
                CALL ppm_error(999,caller,   &
                    &  'need to define level_id first',&
                    &  __LINE__,info)
                GOTO 9999
            ENDIF
            IF (Particles%level_grad_id .EQ. 0) THEN
                info = ppm_error_error
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
                    info = ppm_error_error
                    CALL ppm_error(999,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ELSE
                IF (Particles%wpv_m(Particles%adapt_wpgradid).NE.1) THEN
                    CALL particles_allocate_wpv(Particles,Particles%adapt_wpgradid,&
                        ppm_dim,info,with_ghosts=.FALSE.,&
                        iopt=ppm_param_alloc_grow_preserve)
                    info = ppm_error_error
                    CALL ppm_error(999,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    ! Stopping criterion for gradient descent
    ! (depends on the potential that is used)
    Psi_threshold = opts%adaptivity_criterion

    !-------------------------------------------------------------------------!
    ! Start the actual computations
    !-------------------------------------------------------------------------!

    !!-------------------------------------------------------------------------!
    !! Compute D (desired resolution)
    !!-------------------------------------------------------------------------!

    CALL sop_compute_D(Particles,D_fun,opts,info,     &
        wp_fun=wp_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
        level_grad_fun=level_grad_fun,D_needs_gradients=D_needs_gradients,&
        nb_fun=nb_fun)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'sop_compute_D failed',info)
        info = -1
        GOTO 9999
    ENDIF

    !REMOVME???
    rcp => Get_wps(Particles,Particles%rcp_id)
    D => Get_wps(Particles,Particles%D_id)

    DO ip=1,Particles%Npart
        rcp(ip) = opts%rcp_over_D * D(ip)
    ENDDO
            IF (opts%level_set) THEN
                IF (PRESENT(wp_fun)) THEN
                    xp => Get_xp(Particles)
                    DO ip=1,Particles%Npart
                        IF (ABS(level_fun(xp(1:ppm_dim,ip))) .GT. &
                &   opts%nb_width2*nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) THEN
                            rcp(ip) = 1.3_MK * rcp(ip)
                        ENDIF
                    ENDDO
                    xp => Set_xp(Particles,read_only=.TRUE.)
                ELSE
                    level => Get_wps(Particles,Particles%level_id)
                    wp => Get_wps(Particles,Particles%adapt_wpid)
                    DO ip=1,Particles%Npart
                        IF (ABS(level(ip)) .GT. &
                            &   opts%nb_width2*nb_fun(wp(ip),opts%scale_D)) THEN
                            rcp(ip) = 1.3_MK * rcp(ip)
                        ENDIF
                    ENDDO
                    level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
                    wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
                ENDIF
            ENDIF

    CALL particles_updated_cutoff(Particles,info)

    rcp => Set_wps(Particles,Particles%rcp_id,read_only=.FALSE.)
    D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)

    CALL particles_mapping_ghosts(Particles,topo_id,info)

    CALL particles_neighlists(Particles,topo_id,info)
    !REMOVME???


    !need the ghosts for D_old to be correct
    CALL particles_mapping_ghosts(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'particles_mapping_ghosts failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Find the nearest neighbour for each real particle
    !! (used in the definition of the primitive functions, for
    !! interpolation. Not needed if the functions are known analytically)
    !!-------------------------------------------------------------------------!
    IF (.NOT.PRESENT(wp_fun)) THEN
        nn2_id = 0
        CALL particles_allocate_wps(Particles,nn2_id,info)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'particles_allocate_wps failed',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF

        nn2=> Get_wps(Particles,nn2_id)
        xp => Get_xp(Particles,with_ghosts=.TRUE.)
        DO ip=1,Particles%Npart
            nn2(ip) = HUGE(1._MK)
            DO ineigh=1,Particles%nvlist(ip)
                iq=Particles%vlist(ineigh,ip)
                dist2 = SUM( (xp(1:ppm_dim,iq)-xp(1:ppm_dim,ip))**2 )
                nn2(ip) = MIN(nn2(ip),dist2)
            ENDDO
        ENDDO
        xp => Set_xp(Particles,read_only=.TRUE.)
        nn2=> Set_wps(Particles,nn2_id)

        !!---------------------------------------------------------------------!
        !! Get ghosts of nn2
        !! TODO: do that more cleanly
        !!---------------------------------------------------------------------!
        CALL ppm_map_part_push(Particles%wps(nn2_id)%vec,Particles%Npart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (push) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_send(Particles%Npart,Particles%Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (send) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        CALL ppm_map_part_pop(Particles%wps(nn2_id)%vec,&
            Particles%Npart,Particles%Mpart,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_map_part_ghost (pop) failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        Particles%wps_g(nn2_id) = 1
    ENDIF

    !-------------------------------------------------------------------------!
    ! Copy particle positions, field values and D
    !-------------------------------------------------------------------------!

    Particles_old => Particles
    Particles => NULL()

    !Particles is allocated with size Mpart - ie WITH the ghosts particles
    CALL ppm_alloc_particles(Particles,Particles_old%Mpart,&
        ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'ppm_alloc_particles failed',info)
        info = -1
        GOTO 9999
    ENDIF

    ! Soft copy
    xp => Particles%xp
    Particles = Particles_old
    Particles%xp => xp

    ! Set all arrays to unmapped
    Particles%wps => NULL()
    !DO i=1,Particles%max_wpsid
        !NULLIFY(Particles%wps(i)%vec)
    !ENDDO
    Particles%wps_m => NULL()
    Particles%wps_g => NULL()
    IF (Particles%nwpv.GT.0) THEN
        !DO i=1,Particles%max_wpvid
            !NULLIFY(Particles%wpv(i)%vec)
        !ENDDO
        Particles%wpv => NULL()
        Particles%wpv_m => NULL()
        Particles%wpv_g => NULL()
        Particles%wpv_s => NULL()
    ENDIF
    Particles%nwps = 0
    Particles%nwpv = 0
    Particles%max_wpsid = 0
    Particles%max_wpvid = 0

    !transfer nvlist and vlist to Particles
    !Particles%neighlists = Particles_old%neighlists
    Particles%nvlist => Particles_old%nvlist
    Particles%vlist => Particles_old%vlist
    Particles_old%nvlist => NULL()
    Particles_old%vlist => NULL()
    Particles_old%neighlists = .FALSE.

    CALL particles_allocate_wps(Particles,Particles%D_id,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_fit)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,'particles_allocate_wps failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    CALL particles_allocate_wps(Particles,Particles%rcp_id,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_fit)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(999,caller,'particles_allocate_wps failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    Particles%nwps=2

    xp => Get_xp(Particles,with_ghosts=.TRUE.)
    D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
    rcp => Get_wps(Particles,Particles%rcp_id,with_ghosts=.TRUE.)
    xp_old => Get_xp(Particles_old,with_ghosts=.TRUE.)
    D_old => Get_wps(Particles_old,Particles_old%D_id,with_ghosts=.TRUE.)
    rcp_old => Get_wps(Particles_old,Particles_old%rcp_id,with_ghosts=.TRUE.)
    DO ip=1,Particles%Mpart
        xp(1:ppm_dim,ip) = xp_old(1:ppm_dim,ip)
        D(ip) = D_old(ip)
        rcp(ip) = rcp_old(ip)
    ENDDO

    xp  => Set_xp(Particles,read_only=.TRUE.)
    D   => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
    rcp => Set_wps(Particles,Particles%rcp_id)
    xp_old  => Set_xp(Particles_old,read_only=.TRUE.)
    D_old   => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)
    rcp_old => Set_wps(Particles_old,Particles_old%rcp_id,read_only=.TRUE.)

    !If we use a level-set method with narrow band
    ! we have to keep track of the distance function
    ! DURING adaptation. We thus have to copy the 
    ! level function (and its gradient)
    IF (Particles_old%level_set .AND. .NOT.PRESENT(level_fun)) THEN
        CALL particles_allocate_wps(Particles,Particles%level_id,&
            info,zero=.TRUE.,iopt=ppm_param_alloc_fit)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'particles_allocate_wps failed',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        Particles%nwps = Particles%nwps + 1
        CALL particles_allocate_wpv(Particles,Particles%level_grad_id,&
            ppm_dim,info,zero=.TRUE.,iopt=ppm_param_alloc_fit)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'particles_allocate_wpv failed',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        Particles%nwpv = Particles%nwpv + 1

        level => Get_wps(Particles,Particles%level_id,with_ghosts=.FALSE.)
        level_grad => Get_wpv(Particles,&
            Particles%level_grad_id,with_ghosts=.FALSE.)
        level_old => Get_wps(Particles_old,Particles_old%level_id)
        level_grad_old => Get_wpv(Particles_old,&
            Particles_old%level_grad_id)
        DO ip=1,Particles%Npart
            level(ip) = level_old(ip)
            level_grad(1:ppm_dim,ip) = level_grad_old(1:ppm_dim,ip)
        ENDDO

        level => Set_wps(Particles,Particles%level_id,read_only=.FALSE.)
        level_grad => Set_wpv(Particles,&
            Particles%level_grad_id,read_only=.FALSE.)
        level_old => Set_wps(Particles_old,&
            Particles_old%level_id,read_only=.TRUE.)
        level_grad_old => Set_wpv(Particles_old,&
            Particles_old%level_grad_id,read_only=.TRUE.)
    ENDIF

    IF (.NOT.PRESENT(wp_fun)) THEN
        CALL particles_allocate_wps(Particles,Particles%adapt_wpid,&
            info,iopt=ppm_param_alloc_fit)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(999,caller,'particles_allocate_wps failed',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF
        Particles%nwps = Particles%nwps + 1

        wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.FALSE.)
        wp_old => Get_wps(Particles_old,Particles_old%adapt_wpid,with_ghosts=.FALSE.)
        DO ip=1,Particles%Npart
            wp(ip)= wp_old(ip)
        ENDDO
        wp => Set_wps(Particles,Particles%adapt_wpid)
        wp_old => Set_wps(Particles_old,Particles_old%adapt_wpid,read_only=.TRUE.)
    ENDIF

    ! Now, work on xp and keep xp_old for interpolation

    !!-------------------------------------------------------------------------!
    !! Move particles
    !! (Do a gradient descent until stopping criterion is met
    !! On output, vlist_cross may NOT have enough elements to 
    !! allow for computing corrected kernels)
    !! On output, particles MUST be within the computational domain
    !! On output, vlist SHOULD be correct
    !!-------------------------------------------------------------------------!

    ! The Intel compiler seems to accept the following call to 
    ! a subroutine with optional arguments, but this might be an issue 
    ! with some compilers (who knows...)
    ! (e.g. wp_grad_fun below may or may not be present)

    IF (Particles%level_set) THEN
        CALL sop_gradient_descent_ls(Particles_old,Particles, &
            nvlist_cross,vlist_cross,                &
            nneighmin_cross,nneighmax_cross,num_it,opts,info,wp_fun=wp_fun,&
            D_fun=D_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
            level_grad_fun=level_grad_fun,threshold=Psi_threshold,&
            need_deriv=need_derivatives,nb_fun=nb_fun)
    ELSE
        CALL sop_gradient_descent(Particles_old,Particles, &
            nvlist_cross,vlist_cross,                &
            nneighmin_cross,nneighmax_cross,num_it,opts,info,wp_fun=wp_fun,&
            D_fun=D_fun,wp_grad_fun=wp_grad_fun,threshold=Psi_threshold,&
            need_deriv=need_derivatives)
    ENDIF
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'sop_gradient_descent failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    IF (PRESENT(stats)) &
        stats%nb_grad_desc_steps = stats%nb_grad_desc_steps + num_it

#if debug_verbosity > 1
    CALL sop_dump_debug(Particles%xp,ppm_dim,Particles%Npart,1004,info)
    CALL sop_dump_debug(Particles%wps(Particles%D_id)%vec,Particles%Npart,1005,info)
#endif

    !!-------------------------------------------------------------------------!
    !! Compute field values at new particle locations
    !!-------------------------------------------------------------------------!

    IF (PRESENT(wp_fun)) THEN
        !!---------------------------------------------------------------------!
        !! If the field is known as a function, simply evaluate at new positions
        !!---------------------------------------------------------------------!

        ! DO NOTHING

        !CALL particles_allocate_wps(Particles,Particles%adapt_wpid,&
            !info,with_ghosts=.FALSE.)
        !IF (info.NE.0) THEN
            !info = ppm_error_error
            !CALL ppm_error(999,caller,'particles_allocate_wps failed',&
                !&  __LINE__,info)
            !GOTO 9999
        !ENDIF

        !xp => Get_xp(Particles)
        !wp => Get_wps(Particles,Particles%adapt_wpid)
        !DO ip=1,Particles%Npart
            !wp(ip)=wp_fun(xp(1:ppm_dim,ip))
        !ENDDO
        !wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.FALSE.)
        !xp => Set_xp(Particles,read_only=.TRUE.)
    ELSE
        !!---------------------------------------------------------------------!
        !! Otherwise, use particle-to-particle interpolation
        !!---------------------------------------------------------------------!
        CALL sop_interpolate(Particles_old,Particles,opts,nn2_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'sop_interpolate failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        
    ENDIF

#if debug_verbosity > 2
    CALL sop_dump_debug(Particles%xp,ppm_dim,Particles%Npart,1100,info)
    wp => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.FALSE.)
    CALL sop_dump_debug(wp,Particles%Npart,1101,info)
    wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
#endif

    !IF (ppm_rank .EQ. 0) THEN
        !WRITE(filename,'(A,A)') TRIM(debugdir),'Stats.dat'
        !iunit=20
        !OPEN(UNIT=iunit,FILE=filename,FORM='FORMATTED',&
            !ACCESS='sequential',POSITION='APPEND',IOSTAT=info)
        !WRITE(iunit,'(I8,2X,I8,2X)') Particles%Npart,nb_grad_desc_steps
        !CLOSE(iunit)
        !IF (info .NE. 0) THEN
            !CALL ppm_write(ppm_rank,caller,&
                !'I/O error, writeout of Stats.dat failed ',info)
            !info = -1
            !GOTO 9999
        !ENDIF
    !ENDIF

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
    !de-activate D_id (instead of de-allocating it)
    !Particles%wps_m(Particles%D_id) = -1

    !CALL particles_allocate_wps(Particles,D_id,info,iopt=ppm_param_dealloc)
    !IF (info .NE. 0) THEN
        !CALL ppm_write(ppm_rank,caller,'particles_allocate_wps (dealloc) failed',info)
        !info = -1
        !GOTO 9999
    !ENDIF

    CALL ppm_alloc_particles(Particles_old,Particles%Npart,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'ppm_alloc_particles(dealloc) failed',info)
        info = -1
        GOTO 9999
    ENDIF

    IF(ASSOCIATED(nvlist_cross)) DEALLOCATE(nvlist_cross)
    IF(ASSOCIATED(vlist_cross)) DEALLOCATE(vlist_cross)
    IF(ALLOCATED(eta_interp)) DEALLOCATE(eta_interp)

    IF(ASSOCIATED(xp))    CALL ppm_write(ppm_rank,caller,'forgot to Set xp',info)
    IF(ASSOCIATED(xp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set xp_old',info)
    IF(ASSOCIATED(D))     CALL ppm_write(ppm_rank,caller,'forgot to Set D',info)
    IF(ASSOCIATED(Dtilde))CALL ppm_write(ppm_rank,caller,'forgot to Set D_tilde',info)
    IF(ASSOCIATED(D_old)) CALL ppm_write(ppm_rank,caller,'forgot to Set D_old',info)
    IF(ASSOCIATED(xp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set xp_old',info)
    IF(ASSOCIATED(wp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set wp_old',info)
    IF(ASSOCIATED(level)) CALL ppm_write(ppm_rank,caller,'forgot to Set level',info)
    IF(ASSOCIATED(level_old))CALL ppm_write(ppm_rank,caller,'forgot to Set level_old',info)
    IF(ASSOCIATED(level_grad))CALL ppm_write(ppm_rank,caller,'forgot to Set level_grad',info)
    IF(ASSOCIATED(level_grad_old)) &
        CALL ppm_write(ppm_rank,caller,'forgot to Set level_grad_old',info)

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error


END SUBROUTINE sop_adapt_particles

#undef __KIND
