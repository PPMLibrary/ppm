SUBROUTINE DTYPE(sop_interpolate)(Particles_old,Particles,opts,info)
    !!!------------------------------------------------------------------------!
    !!! Interpolate the field variables from the old particles' positions 
    !!! to the new ones 
    !!!------------------------------------------------------------------------!

    !-------------------------------------------------------------------------
    !  Modules
    !-------------------------------------------------------------------------
    USE ppm_module_inl_xset_k_vlist
    USE ppm_module_dcops
    USE ppm_module_io_vtk

    IMPLICIT NONE

    DEFINE_MK()
    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(DTYPE(ppm_t_particles)), POINTER,INTENT(IN   )   :: Particles_old
    !!! Old set of particles
    TYPE(DTYPE(ppm_t_particles)), POINTER,INTENT(INOUT)   :: Particles
    !!! New set of particles
    TYPE(DTYPE(sop_t_opts)),  POINTER,    INTENT(IN   )   :: opts
    !!! Options
    INTEGER,                              INTENT(  OUT)   :: info
    !!! Return status, 0 upon success

    !-------------------------------------------------------------------------
    !  Local variables
    !-------------------------------------------------------------------------
    INTEGER                         :: ip,iq,ineigh,i,prop_id,eta_id
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
    INTEGER                                    :: memory_used


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
    !---------------------------------------------------------------------!
    ! We need to ensure cross-neighbour lists are large enough 
    ! to compute the interpolation kernels
    !   We use kd-trees to do that.
    !---------------------------------------------------------------------!

    CALL particles_neighlists_xset(Particles,Particles_old,&
        Particles%active_topoid,info)!,knn=opts%nneigh_critical)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_sub_failed,caller,&
            'particles_neighlists_xset failed.',__LINE__,info)
        GOTO 9999
    ENDIF

    !!---------------------------------------------------------------------!
    !! Compute interpolation kernels
    !!---------------------------------------------------------------------!
    IF (opts%order_approx .GE. 0) THEN
        ALLOCATE(order(1),degree(ppm_dim))
        order = opts%order_approx
        degree = 0 !zeroth-order derivative => interpolation
        eta_id = 0
        CALL particles_dcop_define(Particles,eta_id,(/1._MK/),degree,&
            order,1,info,name="interp",interp=.TRUE.)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_define failed',__LINE__,info)
            GOTO 9999
        ENDIF
        DEALLOCATE(order,degree)
        CALL particles_dcop_compute(Particles,eta_id,info,c=opts%c)
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

    !---------------------------------------------------------------------!
    ! Interpolate scalar fields onto new positions
    ! (reallocate arrays if number of particles has changed)
    !---------------------------------------------------------------------!
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

            CALL particles_dcop_apply(Particles,i,prop_id,eta_id,info)
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
            Particles%level_id,eta_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_apply failed',__LINE__,info)
            GOTO 9999
        ENDIF

        !MAJOR FIXME!!!
        ! do something better to compute the gradients of level
        CALL particles_dcop_free(Particles,eta_id,info)
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
        CALL particles_dcop_define(Particles,eta_id,coeffs,degree,&
            order,ppm_dim,info,name="interp",interp=.TRUE.)
        IF (info.NE.0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_dcop_define failed',__LINE__,info)
            GOTO 9999
        ENDIF
        DEALLOCATE(order,degree)
        CALL particles_dcop_compute(Particles,eta_id,info,c=opts%c)
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
        eta            => Get_dcop(Particles,eta_id)
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
        eta            => Set_dcop(Particles,eta_id)
        xp_old         => Set_xp(Particles_old,read_only=.TRUE.)
        xp             => Set_xp(Particles,read_only=.TRUE.)
        nvlist_cross   => NULL()
        vlist_cross    => NULL()
    ENDIF

    memory_used = 0
    DO i=1,Particles%ops%max_opsid
        IF (ASSOCIATED(Particles%ops%ker(i)%vec)) THEN
            memory_used = memory_used + 8*SIZE(Particles%ops%ker(i)%vec)
        ENDIF
    ENDDO
    !WRITE(cbuf,*) '      dc_ops: ', memory_used
    memory_used_total = memory_used_total+memory_used

    !-------------------------------------------------------------------------!
    ! Free DC operator
    !-------------------------------------------------------------------------!
    CALL particles_dcop_free(Particles,eta_id,info)
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

END SUBROUTINE DTYPE(sop_interpolate)
