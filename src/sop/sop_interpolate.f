!!!----------------------------------------------------------------------------!
!!! Interpolate the field variable from the old particles' positions 
!!!      to the new ones 
!!! Assumes that the operators have been computed by&
!!!     correct_diff_operators_interp
!!! 
!!!----------------------------------------------------------------------------!
SUBROUTINE sop_interpolate(Particles_old,Particles,opts,nn2_id,info)

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_inl_xset_vlist
    USE ppm_module_dcops

    IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,       INTENT(IN   )   :: Particles_old
    !!! Old set of particles
    TYPE(ppm_t_particles), POINTER,       INTENT(INOUT)   :: Particles
    !!! New set of particles
    TYPE(sop_t_opts),  POINTER,           INTENT(IN   )   :: opts
    !!! Options
    INTEGER,                              INTENT(IN   )   :: nn2_id
    !!! index where nearest neighbour distances for Particles_old are stored
    INTEGER,                              INTENT(  OUT)   :: info
    !!! Return status, 0 upon success

    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                         :: ip,iq,ineigh,i,prop_id
    REAL(KIND(1.d0))                :: t0
    CHARACTER(LEN=256)              :: filename,cbuf
    CHARACTER(LEN=256)              :: caller = 'sop_interpolate'

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old
    REAL(MK),     DIMENSION(:),   POINTER      :: wp_old
    REAL(MK),     DIMENSION(:),   POINTER      :: D_old
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

    REAL(MK),     DIMENSION(:),   POINTER      :: nn2 => NULL()
    INTEGER,      DIMENSION(:),   POINTER      :: nvlist_cross => NULL()
    INTEGER,      DIMENSION(:,:), POINTER      :: vlist_cross => NULL()
    !REAL(MK),     DIMENSION(:,:), ALLOCATABLE  :: eta,eta_interp
    REAL(MK),     DIMENSION(:,:), POINTER      :: eta
    !INTEGER                                    :: nneighmax_cross
    !INTEGER                                    :: nneighmin_cross
    !should be removed once the argument lists for the inl routines
    !have been updated to inhomogeneous ghostlayers
    REAL(MK),DIMENSION(2*ppm_dim)       :: ghostlayer
    INTEGER,     DIMENSION(ppm_dim)            :: order_d


    !!-------------------------------------------------------------------------!
    !! Initialise
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    !------------------------------------------------------------------
    ! Since particles have moved during gradient descent, 
    !we need to recompute the cross neigbour lists
    !------------------------------------------------------------------
    !!---------------------------------------------------------------------!
    !! FIXME: We need to ensure cross-neighbour lists are large enough 
    !! to compute the interpolation kernels
    !!---------------------------------------------------------------------!
    !!-----------------------------------------------------------------!
    !! Construct cross-set neighbour lists
    !! (note: D_old is used as rcp_old=rcp_over_D*D_old
    !!  on output, D_old may have been changed artificially to increase
    !! rcp_old. Do not use it anymore (except for computing rcp_old))
    !!-----------------------------------------------------------------!
    D_old => Get_wps(Particles_old,Particles_old%D_id)
    D_old = D_old * opts%rcp_over_D
    ghostlayer=Particles%cutoff
    CALL ppm_inl_xset_vlist(Particles%active_topoid,Particles%xp,Particles%Npart,&
        Particles%Mpart,Particles_old%xp,Particles_old%Npart,&
        Particles_old%Npart,D_old,Particles%skin,&
        ghostlayer,info,Particles%vlist_cross,Particles%nvlist_cross)
    Particles%neighlists_cross = .TRUE.
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,&
            'ppm_inl_xset_vlist failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    D_old = D_old / opts%rcp_over_D
    D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)

    Particles%nneighmin_cross = &
        MINVAL(Particles%nvlist_cross(1:Particles%Npart))
    Particles%nneighmax_cross = &
        MAXVAL(Particles%nvlist_cross(1:Particles%Npart))

    write(cbuf,*) 'Nvlist_cross min/max = ',Particles%nneighmin_cross,&
        Particles%nneighmax_cross
    CALL ppm_write(ppm_rank,caller,cbuf,info)

#if debug_verbosity > 0
    IF (Particles%nneighmin_cross .LT. opts%nneigh_critical) THEN
        WRITE(cbuf,*) 'Too few cross-neighbours, something is wrong'
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info= -1
        GOTO 9999
    ENDIF
#endif

    !!---------------------------------------------------------------------!
    !! Compute interpolation kernels
    !!---------------------------------------------------------------------!
    D_old  => Get_wps(Particles_old,Particles_old%D_id)
    nn2    => Get_wps(Particles_old,nn2_id)
    xp_old => Get_xp(Particles_old)
    IF (opts%order_approx .GE. 0) THEN
        order_d = 0
        CALL ppm_part_dcops(Particles,Particles%eta_id,opts%c,info,&
            order_deriv=order_d,order_approx=opts%order_approx,&
            isinterp=.TRUE.,Particles_old=Particles_old,nn_sq_id=nn2_id)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_part_dcops failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        !CALL sop_dcop_order2_interp(xp_old,D_old,&
            !nn2,Particles%xp,Particles%wps(Particles%rcp_id)%vec,&
            !Particles%Npart,Particles%Mpart,Particles%nvlist,Particles%vlist, &
            !nvlist_cross,vlist_cross,eta,eta_interp,c,&
            !Particles%nneighmin,Particles%nneighmax,nneighmin_cross,&
            !nneighmax_cross,info) 
    ELSE
        CALL ppm_write(ppm_rank,caller,&
            'invalid interpolation order',info)
        info = -1
        GOTO 9999
    ENDIF
    nn2   => Set_wps(Particles_old,nn2_id,read_only=.TRUE.)
    D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)
    xp_old => Set_xp(Particles_old,read_only=.TRUE.)

    !eta_interp => Get_wpv(Particles,Particles%eta_id)

    !!---------------------------------------------------------------------!
    !! Interpolate fields onto new positions
    !! (reallocate arrays if number of particles has changed)
    !!---------------------------------------------------------------------!
    ! Loop through all properties i for which wps_m(i) = 1
    DO i=1,Particles_old%max_wpsid
        IF (Particles_old%wps_m(i).EQ.1) THEN
            !skip the properties that do not need to be interpolated
            IF (i.EQ.nn2_id) CYCLE
            IF (i.EQ.Particles_old%level_id) CYCLE
            IF (i.EQ.Particles_old%D_id) CYCLE
            IF (i.EQ.Particles_old%rcp_id) CYCLE
            IF (i.EQ.Particles_old%Dtilde_id) CYCLE
            prop_id = 0
            CALL particles_allocate_wps(Particles,prop_id,info)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'particles_allocate_wps failed.',info)
                info = -1
                GOTO 9999
            ENDIF
            !wp     => Get_wps(Particles,prop_id)
            !wp_old => Get_wps(Particles_old,i,with_ghosts=.TRUE.)
            !xp     => Get_xp(Particles)
            !xp_old => Get_xp(Particles,with_ghosts=.TRUE.)

            CALL particles_apply_dcops(Particles,prop_id,prop_id,Particles%eta_id,0,&
                info,Particles_old,i)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'particles_apply_dcops failed.',info)
                info = -1
                GOTO 9999
            ENDIF
            !CALL sop_interpolate_particles(&
                !xp,wp_old,xp,Particles%Npart,nvlist_cross,&
                !vlist_cross,eta,eta_interp,opts,info,wp=wp)
            !IF (info .NE. 0) THEN
                !CALL ppm_write(ppm_rank,caller,'sop_interp_particles failed.',info)
                !info = -1
                !GOTO 9999
            !ENDIF
            !wp     => Set_wps(Particles,prop_id)
            !wp_old => Set_wps(Particles_old,i,read_only=.TRUE.)
            !xp     => Set_xp(Particles,read_only=.TRUE.)
            !xp_old => Set_xp(Particles,read_only=.TRUE.)
        ENDIF
    ENDDO


    IF (opts%level_set) THEN
            CALL particles_apply_dcops(Particles,Particles%level_id,&
                Particles%level_id,Particles%eta_id,0,info,&
                Particles_old,Particles_old%level_id)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'particles_apply_dcops failed.',info)
                info = -1
                GOTO 9999
            ENDIF

            !MAJOR FIXME!!!
            ! do something better to compute the gradients of level
        CALL ppm_part_dcops(Particles,Particles%eta_id,opts%c,info,&
            islaplacian=.TRUE.,isinterp=.TRUE.,&
            Particles_old=Particles_old,nn_sq_id=nn2_id)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'ppm_part_dcops failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        level          => Get_wps(Particles,Particles%level_id)
        level_grad     => Get_wpv(Particles,Particles%level_grad_id)
        level_old      => Get_wps(Particles_old,Particles_old%level_id)
        xp_old         => Get_xp(Particles_old)
        xp             => Get_xp(Particles)
        eta            => Get_wpv(Particles,Particles%eta_id)
        vlist_cross    => Particles%vlist_cross
        nvlist_cross   => Particles%nvlist_cross
            DO ip = 1,Particles%Npart 
                level_grad(1:ppm_dim,ip) = 0._MK
            ENDDO
            DO ip = 1,Particles%Npart
                DO ineigh = 1,nvlist_cross(ip)
                    iq = vlist_cross(ineigh,ip)
                    level_grad(1:ppm_dim,ip) = level_grad(1:ppm_dim,ip) + &
                        (level(ip)+level_old(iq)) * 0.5_MK*       & 
                        (xp_old(1:ppm_dim,iq)-xp(1:ppm_dim,ip)) * eta(ineigh,ip)
                ENDDO
            ENDDO
        level          => Set_wps(Particles,Particles%level_id)
        level_grad     => Set_wpv(Particles,Particles%level_grad_id)
        level_old      => Set_wps(Particles_old,Particles_old%level_id,&
            read_only=.TRUE.)
        eta            => Set_wpv(Particles,Particles%eta_id,read_only=.TRUE.)
        xp_old         => Set_xp(Particles_old,read_only=.TRUE.)
        xp             => Set_xp(Particles,read_only=.TRUE.)
        nvlist_cross   => NULL()
        vlist_cross    => NULL()

        !level          => Get_wps(Particles,Particles%level_id)
        !level_grad     => Get_wpv(Particles,Particles%level_grad_id)
        !level_old      => Get_wps(Particles_old,Particles_old%level_id)
        !CALL sop_interpolate_particles(&
            !Particles_old%xp,level_old,&
            !Particles%xp,Particles%Npart,nvlist_cross,&
            !vlist_cross,eta,eta_interp,opts,info,wp=level,wp_grad=level_grad)
        !level          => Set_wps(Particles,Particles%level_id)
        !level_grad     => Set_wpv(Particles,Particles%level_grad_id)
        !level_old      => Set_wps(Particles_old,Particles_old%level_id,&
            !read_only=.TRUE.)
        !IF (info .NE. 0) THEN
            !CALL ppm_write(ppm_rank,caller,'sop_interp_particles failed.',info)
            !info = -1
            !GOTO 9999
        !ENDIF
    ENDIF
    !eta_interp => Set_wpv(Particles,Particles%eta_id,read_only=.TRUE.)

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999  CONTINUE ! jump here upon error

END SUBROUTINE sop_interpolate

#undef __KIND
