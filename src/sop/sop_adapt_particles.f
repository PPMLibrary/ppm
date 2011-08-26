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
        wp_fun,wp_grad_fun,wplap_id,return_grad,smallest_sv,stats)

    USE ppm_module_error
    USE ppm_module_map_part
    USE ppm_module_io_vtk

    IMPLICIT NONE

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
    REAL(MK),OPTIONAL,                    INTENT(INOUT)   :: smallest_sv
    !!! smallest singular value of the Vandermonde matrices (DC-PSE)
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

    END INTERFACE

    ! local variables
    INTEGER                                    :: i,ip,ineigh,iq
    INTEGER                                    :: iunit
    CHARACTER(LEN=256)                         :: filename,cbuf
    CHARACTER(LEN=256)                         :: caller='sop_adapt_particles'
    REAL(KIND(1.D0))                           :: t0
    REAL(MK)                                   :: dist2
    REAL(MK)                                   :: Psi_threshold
    INTEGER                                    :: num_it

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
    REAL(MK),     DIMENSION(:,:),   POINTER      :: D_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level_old => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad_old => NULL()

    TYPE(ppm_t_particles), POINTER             :: Particles_old

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp => NULL()
    REAL(MK),     DIMENSION(:,:),   POINTER    :: D => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: wp_grad => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad => NULL()


    REAL(MK),     DIMENSION(:),   POINTER      :: nn2 => NULL()
    INTEGER,      DIMENSION(:),   POINTER      :: nvlist_cross => NULL()
    INTEGER,      DIMENSION(:,:), POINTER      :: vlist_cross => NULL()
    INTEGER                                    :: nneighmax_cross
    INTEGER                                    :: nneighmin_cross
    LOGICAL                                    :: need_derivatives
    
    INTEGER                                    :: tensor_length
    REAL(MK),     DIMENSION(:,:), POINTER      :: inv => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: inv_old => NULL()

    !-------------------------------------------------------------------------!
    ! Initialize
    !-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    adapt_wpgradid = 0
    !FIXME: make it possible to pass this id as an argument (to avoid
    ! reallocating it every time, and to enable re-using the values
    ! outside of the routine)

    !-------------------------------------------------------------------------!
    ! Checks consistency of parameters
    !-------------------------------------------------------------------------!
    IF (topo_id .NE. Particles%active_topoid) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Particles need to be on the same topology as &
            &   the one refered by the topo_id argument',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (PRESENT(wp_grad_fun) .AND. .NOT.PRESENT(wp_fun)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'provided analytical gradients but not analytical &
            &   function values. This case is not yet implemented',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT. ASSOCIATED(opts)) THEN
        !provide defaults for opts
        CALL sop_init_opts(opts,info)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,   &
                &  'Allocation failed for opts',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF
    IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
        IF (opts%D_needs_gradients) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                & 'provided analytical function values but no &
                & function gradients. This case is not implemented',&
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
    IF (.NOT.Particles%anisotropic) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'These particles have not been declared as anisotropic',&
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

    !check that neighbour lists have been computed
    IF (.NOT.Particles%neighlists) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,   &
            &  'Neighbour lists are not up-to-date. Call neighlists first',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------!
    ! Determines whether we will need to approximate derivatives
    !-------------------------------------------------------------------------!
    need_derivatives=.FALSE.
    IF (.NOT.PRESENT(wp_grad_fun)) THEN
        IF (Particles%level_set .OR. opts%D_needs_gradients) & 
            need_derivatives=.TRUE.
    ENDIF

    ! if we give stats then init the minimum singular value
    IF (PRESENT(stats)) stats%min_sv=HUGE(1._mk)

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
        IF (.NOT.Particles%wps(Particles%adapt_wpid)%has_ghosts) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,   &
                &  'need to get the ghosts for adapt_wpid',&
                &  __LINE__,info)
            GOTO 9999
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
                    CALL ppm_error(ppm_err_alloc,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ELSE
                IF (.NOT.Particles%wpv(adapt_wpgradid)%is_mapped) THEN
                    CALL particles_allocate_wpv(Particles,adapt_wpgradid,&
                        ppm_dim,info,with_ghosts=.FALSE.,&
                        iopt=ppm_param_alloc_grow_preserve,name='adapt_wpgrad')
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,'particles_allocate_wpv failed',&
                        &  __LINE__,info)
                    GOTO 9999
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    ! Stopping criterion for gradient descent
    ! (depends on the potential that is used)
    ! this is the threshold for the inverse distance, i.e. 2 -> all need to be farer away than 0.5
    Psi_threshold = opts%adaptivity_criterion

    !-------------------------------------------------------------------------!
    ! Start the actual computations
    !-------------------------------------------------------------------------!
    !!-------------------------------------------------------------------------!
    !! Compute the anisotropic requirements
    !!-------------------------------------------------------------------------!
    CALL sop_compute_D(Particles,D_fun,opts,info,     &
        wp_fun=wp_fun,wp_grad_fun=wp_grad_fun,stats=stats)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'sop_computed_D failed',__LINE__,info)
        GOTO 9999
    ENDIF

    ! real tensors have now Dtilde
    ! Dtilde and D are stored in wpv 

    IF (.FALSE.) THEN

      ! haeckic: add the neighbor consideration option

      ! If we use neighbor consideration, then make this
!     rcp => Get_wps(Particles,Particles%rcp_id)
!     D => Get_wps(Particles,Particles%D_id)
! 
!     DO ip=1,Particles%Npart
!         rcp(ip) = opts%rcp_over_D * D(ip)
!     ENDDO
! 
!     rcp => Set_wps(Particles,Particles%rcp_id,read_only=.FALSE.)
!     D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)

    ENDIF


#if debug_verbosity > 1
    WRITE(filename,'(A,I0)') 'P_beforegraddesc_',Particles%itime
    CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif

    ! Problem: what if cutoff is now bigger than initially???? should not happen...
    CALL particles_updated_cutoff(Particles,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'particles_updated_cutoff failed',__LINE__,info)
        GOTO 9999
    ENDIF

    ! Problem: The ghost mapping needs to consider D_old cutoff, because they are used!!
    ! need the ghosts for D_old to be correct
    CALL particles_mapping_ghosts(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'particles_mapping_ghosts failed',__LINE__,info)
        GOTO 9999
    ENDIF

    CALL particles_neighlists(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'particles_neighlists failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !REMOVME???

    !removme
    !CALL sop_compute_D(Particles,D_fun,opts,info,     &
        !wp_fun=wp_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
        !level_grad_fun=level_grad_fun,nb_fun=nb_fun)
    !WRITE(filename,'(A,I0)') 'P_beforegraddescbis_',Particles%itime
    !CALL ppm_vtk_particle_cloud(filename,Particles,info)
    !stop
    !removme

    !!-------------------------------------------------------------------------!
    !! Find the nearest neighbour for each real particle
    !! (used in the definition of the primitive functions, for
    !! interpolation. Not needed if the functions are known analytically)
    !!-------------------------------------------------------------------------!
    IF (.NOT.PRESENT(wp_fun)) THEN
        CALL particles_allocate_wps(Particles,Particles%nn_sq_id,info,name='nn2')
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,'particles_allocate_wps failed',&
                &  __LINE__,info)
            GOTO 9999
        ENDIF

        nn2=> Get_wps(Particles,Particles%nn_sq_id)
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
        nn2=> Set_wps(Particles,Particles%nn_sq_id)

        !!---------------------------------------------------------------------!
        !! Get ghosts of nn2
        !!---------------------------------------------------------------------!
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,   &
                &    'particles_mapping_ghosts failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF

    !-------------------------------------------------------------------------!
    ! Copy particle positions, field values and D
    !-------------------------------------------------------------------------!
    !IF (ASSOCIATED(Particles%ops)) THEN
        !CALL particles_dcop_deallocate(Particles,info)
        !IF (info .NE. 0) THEN
            !info = ppm_error_error
            !CALL ppm_error(ppm_err_dealloc,caller,   &
                !&    'particles_dcop_deallocate failed',__LINE__,info)
            !GOTO 9999
        !ENDIF
    !ENDIF

   ! haeckic: do the stuff for the new particles for anisotropic

    Particles_old => Particles
    Particles => NULL()

    !Particles is allocated with size Mpart - ie WITH the ghosts particles
    CALL ppm_alloc_particles(Particles,Particles_old%Mpart,&
        ppm_param_alloc_fit,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,   &
            &    'ppm_alloc_particles failed',__LINE__,info)
        GOTO 9999
    ENDIF

    ! Soft copy
    xp => Particles%xp
    Particles = Particles_old
    Particles%xp => xp
    !move DC operators from old to new particles
    ! (only move their definitions
! !     Particles%ops => Particles_old%ops
! !     Particles_old%ops => NULL()
! !     IF (ASSOCIATED(Particles%ops)) THEN
! !         DO i=1,Particles%ops%max_opsid
! !             Particles%ops%desc(i)%is_computed = .FALSE.
! !         ENDDO
! !     ENDIF

    ! Set all arrays to unmapped
    Particles%wpi => NULL()
    Particles%wps => NULL()
    IF (Particles%nwpv.GT.0) THEN
        Particles%wpv => NULL()
    ENDIF
    Particles%nwpi = 0
    Particles%nwps = 0
    Particles%nwpv = 0
    Particles%max_wpiid = 0
    Particles%max_wpsid = 0
    Particles%max_wpvid = 0


    !transfer nvlist and vlist to Particles
    !Particles%neighlists = Particles_old%neighlists
    Particles%nvlist => Particles_old%nvlist
    Particles%vlist => Particles_old%vlist

    Particles_old%nvlist => NULL()
    Particles_old%vlist => NULL()
    Particles_old%neighlists = .FALSE.

    !link the old particles to the new ones in the data structure:
    Particles%Particles_cross => Particles_old

    IF (ppm_dim .EQ. 2) THEN
         tensor_length = 4
    ELSE
         tensor_length = 9
    ENDIF

    CALL particles_allocate_wpv(Particles,Particles%D_id,tensor_length,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_fit)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'particles_allocate_wpv failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    Particles%wpv(Particles%D_id)%name = Particles_old%wpv(Particles_old%D_id)%name
    CALL particles_allocate_wpv(Particles,Particles%G_id,tensor_length,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_fit)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'particles_allocate_wpv failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF
    Particles%wpv(Particles%G_id)%name = Particles_old%wpv(Particles_old%G_id)%name
    Particles%nwpv=2

    xp => Get_xp(Particles,with_ghosts=.TRUE.)
    D => Get_wpv(Particles,Particles%D_id,with_ghosts=.TRUE.)
    inv => Get_wpv(Particles,Particles%G_id,with_ghosts=.TRUE.)
    xp_old => Get_xp(Particles_old,with_ghosts=.TRUE.)
    D_old => Get_wpv(Particles_old,Particles_old%D_id,with_ghosts=.TRUE.)
    inv_old => Get_wpv(Particles_old,Particles_old%G_id,with_ghosts=.TRUE.)
    ! copy the properties 
    DO ip=1,Particles%Mpart

         xp(1:ppm_dim,ip) = xp_old(1:ppm_dim,ip)

         IF (ppm_dim .EQ. 2) THEN
            D(1,ip)   = D_old(1,ip)
            D(2,ip)   = D_old(2,ip)
            D(3,ip)   = D_old(3,ip)
            D(4,ip)   = D_old(4,ip)
            inv(1,ip) = inv_old(1,ip)
            inv(2,ip) = inv_old(2,ip)
            inv(3,ip) = inv_old(3,ip)
            inv(4,ip) = inv_old(4,ip)
         ELSE
            D(1,ip)   = D_old(1,ip)
            D(2,ip)   = D_old(2,ip)
            D(3,ip)   = D_old(3,ip)
            D(4,ip)   = D_old(4,ip)
            D(5,ip)   = D_old(5,ip)
            D(6,ip)   = D_old(6,ip)
            D(7,ip)   = D_old(7,ip)
            D(8,ip)   = D_old(8,ip)
            D(9,ip)   = D_old(9,ip)
            inv(1,ip) = inv_old(1,ip)
            inv(2,ip) = inv_old(2,ip)
            inv(3,ip) = inv_old(3,ip)
            inv(4,ip) = inv_old(4,ip)
            inv(5,ip) = inv_old(5,ip)
            inv(6,ip) = inv_old(6,ip)
            inv(7,ip) = inv_old(7,ip)
            inv(8,ip) = inv_old(8,ip)
            inv(9,ip) = inv_old(9,ip)
         ENDIF

    ENDDO

    xp  => Set_xp(Particles,read_only=.TRUE.)
    D   => Set_wpv(Particles,Particles%D_id,read_only=.TRUE.)
    inv => Set_wpv(Particles,Particles%G_id)
    xp_old  => Set_xp(Particles_old,read_only=.TRUE.)
    D_old   => Set_wpv(Particles_old,Particles_old%D_id,read_only=.TRUE.)
    inv_old => Set_wpv(Particles_old,Particles_old%G_id,read_only=.TRUE.)

!     write(*,*) 'particles',Particles%wpv(1)%name, Particles%wpv(2)%name,Particles%wpv(3)%name
!     write(*,*) Particles%max_wpvid,Particles%nwps, Particles%nwpv
!     write(*,*) Particles%D_id, Particles%Dtilde_id,Particles%G_id
!     
!     write(*,*) 'particles2',Particles_old%wpv(1)%name, Particles_old%wpv(2)%name,Particles_old%wpv(3)%name
!     write(*,*) Particles_old%max_wpvid,Particles_old%nwps, Particles_old%nwpv
!     write(*,*) Particles_old%D_id, Particles_old%Dtilde_id,Particles_old%G_id

    ! Now, work on xp and keep xp_old for interpolation
    ! 
    ! Particles:
    ! - xp
    ! - D
    ! - G_id
    ! 
    ! Particles_old:
    ! - xp
    ! - wp
    ! - D
    ! - Dtilde
    ! - G_id
    ! 


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

    ! haeckic: do the gradient descent
    CALL sop_gradient_descent(Particles_old,Particles, &
      nvlist_cross,vlist_cross,                &
      nneighmin_cross,nneighmax_cross,num_it,opts,info,wp_fun=wp_fun,&
      D_fun=D_fun,wp_grad_fun=wp_grad_fun,threshold=Psi_threshold,&
      need_deriv=need_derivatives,stats=stats)

    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,   &
            &    'sop_gradient_descent failed',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (PRESENT(stats)) &
        stats%nb_grad_desc_steps = stats%nb_grad_desc_steps + num_it

#if debug_verbosity > 2
    CALL sop_dump_debug(Particles%xp,ppm_dim,Particles%Npart,1004,info)
    CALL sop_dump_debug(Particles%wps(Particles%D_id)%vec,Particles%Npart,1005,info)
#endif

    !!-------------------------------------------------------------------------!
    !! Compute field values at new particle locations
    !!-------------------------------------------------------------------------!

#if debug_verbosity > 1
    WRITE(filename,'(A,I0)') 'P_aftergraddesc_',Particles%itime
    CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif

    IF (PRESENT(wp_fun)) THEN
        !!---------------------------------------------------------------------!
        !! If the field is known as a function, do nothing
        !!---------------------------------------------------------------------!
    ELSE
        !!---------------------------------------------------------------------!
        !! Otherwise, use particle-to-particle interpolation
        !!---------------------------------------------------------------------!
        CALL sop_interpolate(Particles_old,Particles,opts,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,   &
                &    'sop_interpolate failed',__LINE__,info)
            GOTO 9999
        ENDIF
        
    ENDIF

#if debug_verbosity > 1
    WRITE(filename,'(A,I0)') 'P_afterinterpolate_',Particles%itime
    CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif

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

    Particles%Particles_cross => NULL()

    CALL ppm_alloc_particles(Particles_old,Particles%Npart,&
        ppm_param_dealloc,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_dealloc,caller,   &
            &    'ppm_alloc_particles failed',__LINE__,info)
        GOTO 9999
    ENDIF

    IF(ASSOCIATED(nvlist_cross)) DEALLOCATE(nvlist_cross)
    IF(ASSOCIATED(vlist_cross)) DEALLOCATE(vlist_cross)

    IF(ASSOCIATED(xp))    CALL ppm_write(ppm_rank,caller,'forgot to Set xp',info)
    IF(ASSOCIATED(xp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set xp_old',info)
    IF(ASSOCIATED(D))     CALL ppm_write(ppm_rank,caller,'forgot to Set D',info)
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
