!!!----------------------------------------------------------------------------!
!!! Interpolate the field variables from the old particles' positions 
!!! to the new ones 
!!!----------------------------------------------------------------------------!
SUBROUTINE sop_interpolate(Particles_old,Particles,opts,info)

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_inl_xset_k_vlist
    USE ppm_module_dcops
    USE ppm_module_io_vtk

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
    REAL(MK),     DIMENSION(:),   POINTER      :: D => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp => NULL()
    !REAL(MK),     DIMENSION(:,:), POINTER      :: wp_grad => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad => NULL()

    INTEGER,      DIMENSION(:),   POINTER      :: nvlist_cross => NULL()
    INTEGER,      DIMENSION(:,:), POINTER      :: vlist_cross => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: eta
    INTEGER,     DIMENSION(:),ALLOCATABLE      :: order
    INTEGER,     DIMENSION(:),ALLOCATABLE      :: degree
    REAL(MK),    DIMENSION(ppm_dim)            :: coeffs
    !should be removed once the argument lists for the inl routines
    !have been updated to inhomogeneous ghostlayers
    REAL(MK),DIMENSION(2*ppm_dim)              :: ghostlayer


    !!-------------------------------------------------------------------------!
    !! Initialise
    !!-------------------------------------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif
    !not necessary, unless this routine is called externally
    Particles%Particles_cross => Particles_old

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
    D_old => Get_wps(Particles_old,Particles_old%D_id,with_ghosts=.TRUE.)
    ghostlayer=Particles%cutoff
    CALL ppm_inl_xset_k_vlist(Particles%active_topoid,Particles%xp,&
        Particles%Npart,Particles%Mpart,Particles_old%xp,Particles_old%Npart,&
        Particles_old%Mpart,D_old,opts%nneigh_critical,&
        ghostlayer,info,Particles%vlist_cross,Particles%nvlist_cross)
    Particles%neighlists_cross = .TRUE.
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'ppm_inl_xset_k_vlist failed.',__LINE__,info)
        GOTO 9999
    ENDIF
    D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)

    Particles%nneighmin_cross = &
        MINVAL(Particles%nvlist_cross(1:Particles%Npart))
    Particles%nneighmax_cross = &
        MAXVAL(Particles%nvlist_cross(1:Particles%Npart))

    !write(cbuf,*) 'Nvlist_cross min/max = ',Particles%nneighmin_cross,&
        !Particles%nneighmax_cross
    !CALL ppm_write(ppm_rank,caller,cbuf,info)

#if debug_verbosity > 0
    IF (Particles%nneighmin_cross .LT. opts%nneigh_critical) THEN
        i=0
        call particles_allocate_wps(Particles,i,info,name='nvlist_cross')
        wp => get_wps(Particles,i)
        FORALL(ip=1:Particles%Npart) wp(ip) = DBLE(Particles%nvlist_cross(ip))
        wp => set_wps(Particles,i,read_only=.true.)
        call ppm_vtk_particle_cloud('P_old_dbg',Particles_old,info)
        call ppm_vtk_particle_cloud('P_new_dbg',Particles,info)

        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Too few cross-neighbours, something wrong',__LINE__,info)
        GOTO 9999
    ENDIF
#endif

    !!---------------------------------------------------------------------!
    !! Compute interpolation kernels
    !!---------------------------------------------------------------------!
    IF (opts%order_approx .GE. 0) THEN
        ALLOCATE(order(1),degree(ppm_dim))
        order = opts%order_approx
        degree = 0 !zeroth-order derivative => interpolation
        CALL particles_dcop_define(Particles,Particles%eta_id,(/1._MK/),degree,&
            order,1,info,name="interp",interp=.TRUE.)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_define failed',__LINE__,info)
            GOTO 9999
        ENDIF
        DEALLOCATE(order,degree)
        CALL particles_dcop_compute(Particles,Particles%eta_id,info,c=opts%c)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_compute failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ELSE
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,'invalid interpolation order',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF


    !!---------------------------------------------------------------------!
    !! Interpolate scalar fields onto new positions
    !! (reallocate arrays if number of particles has changed)
    !!---------------------------------------------------------------------!
    ! Loop through all properties i that are mapped
    DO i=1,Particles_old%max_wpsid
        IF (Particles_old%wps(i)%is_mapped) THEN
            !skip the properties that do not need to be interpolated
            IF (i.EQ.Particles_old%nn_sq_id) CYCLE
            IF (i.EQ.Particles_old%level_id) CYCLE
            IF (i.EQ.Particles_old%D_id) CYCLE
            IF (i.EQ.Particles_old%rcp_id) CYCLE
            IF (i.EQ.Particles_old%Dtilde_id) CYCLE
            IF (.NOT.Particles_old%wps(i)%map_parts) CYCLE
            prop_id = i
            CALL particles_allocate_wps(Particles,prop_id,info,&
                name=Particles_old%wps(i)%name)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_alloc,caller,&
                    'particles_allocate_wps failed',__LINE__,info)
                GOTO 9999
            ENDIF

            CALL particles_dcop_apply(Particles,i,prop_id,Particles%eta_id,info)
            IF (info .NE. 0) THEN
                info = ppm_error_error
                CALL ppm_error(ppm_err_sub_failed,caller,&
                    'particles_dcop_apply failed',__LINE__,info)
                GOTO 9999
            ENDIF

        ENDIF
    ENDDO


    IF (opts%level_set) THEN
        CALL particles_dcop_apply(Particles,Particles_old%level_id,&
            Particles%level_id,Particles%eta_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_apply failed',__LINE__,info)
            GOTO 9999
        ENDIF

        !MAJOR FIXME!!!
        ! do something better to compute the gradients of level
        CALL particles_dcop_free(Particles,Particles%eta_id,info)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_free failed',__LINE__,info)
            GOTO 9999
        ENDIF
        ALLOCATE(order(ppm_dim),degree(ppm_dim**2))
        order = 2
        !define Laplacian operator
        degree = 0
        FORALL(i=1:ppm_dim) degree((i-1)*ppm_dim+i)=2 !Laplacian
        coeffs = 1._MK
        CALL particles_dcop_define(Particles,Particles%eta_id,coeffs,degree,&
            order,ppm_dim,info,name="interp",interp=.TRUE.)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_define failed',__LINE__,info)
            GOTO 9999
        ENDIF
        DEALLOCATE(order,degree)
        CALL particles_dcop_compute(Particles,Particles%eta_id,info,c=opts%c)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_compute failed',__LINE__,info)
            GOTO 9999
        ENDIF
        level          => Get_wps(Particles,Particles%level_id)
        level_grad     => Get_wpv(Particles,Particles%level_grad_id)
        level_old      => Get_wps(Particles_old,Particles_old%level_id)
        xp_old         => Get_xp(Particles_old)
        xp             => Get_xp(Particles)
        eta            => Get_dcop(Particles,Particles%eta_id)
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
        eta            => Set_dcop(Particles,Particles%eta_id)
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

    !-------------------------------------------------------------------------!
    ! Free DC operator
    !-------------------------------------------------------------------------!
    CALL particles_dcop_free(Particles,Particles%eta_id,info)
    IF (info.NE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,'particles_dcop_free failed',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999  CONTINUE ! jump here upon error

END SUBROUTINE sop_interpolate

#undef __KIND
