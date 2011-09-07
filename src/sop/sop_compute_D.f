!!!----------------------------------------------------------------------------!
!!! Computes anisotropic requirements for given Particles
!!! 
!!! This functions computes Dtilde and D for anisotropic particles.
!!! This is done in the following way. The gradient direction is the
!!! shorter axis of the ellipse, which is scaled with D_fun. The orthogonal
!!! direction is scaled s.t. the value difference is considered.
!!! 
!!!----------------------------------------------------------------------------!

    pure function f0_grad_fun(pos)

        use ppm_module_data, ONLY: ppm_dim
    integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
    real(mk),parameter              :: pi = 3.1415926535897931_mk
        real(mk), dimension(ppm_dim)             :: f0_grad_fun
        real(mk), dimension(ppm_dim), intent(in) :: pos
        real(mk)                                 :: beta

        beta = 500.0_mk
        f0_grad_fun(1) = -beta*2*(pos(1)-0.5)*exp(-beta*((pos(1)-0.5)**2))
        f0_grad_fun(2) = 0

    end function f0_grad_fun

pure function f0_fun(pos)

        use ppm_module_data, ONLY: ppm_dim
    integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
    real(mk),parameter              :: pi = 3.1415926535897931_mk
        real(mk)                                 :: f0_fun
        real(mk), dimension(ppm_dim), intent(in) :: pos
        real(mk)                                 :: beta
         
        beta = 500.0_mk
        f0_fun = exp(-beta*((pos(1)-0.5)**2))

    end function f0_fun

SUBROUTINE sop_compute_D(Particles,D_fun,opts,info,     &
                            wp_fun,wp_grad_fun,stats)

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
    CHARACTER(LEN = 64)                        :: myformat
    CHARACTER(LEN = 256)                       :: cbuf
    CHARACTER(LEN = 256)                       :: caller='sop_compute_req'
    REAL(KIND(1.D0))                           :: t0

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: D_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: rcp_old => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level_old => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad_old => NULL()


    ! HAECKIC: some variables adapted

    REAL(MK),     DIMENSION(:,:), POINTER      :: xp => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: inv => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: rcp => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: D => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: Dtilde => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: wp => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: wp_grad => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: level => NULL()
    REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_A => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_B => NULL()

    REAL(MK),     DIMENSION(:,:), POINTER      :: eta => NULL()
    REAL(MK)                                   :: min_D,orth_len,old_scale,new_scale,new_scale_long,temp_scale,max_g,max_w
    LOGICAL                                    :: need_derivatives
    REAL(MK),     DIMENSION(ppm_dim)           :: dummy_grad
    INTEGER                                    :: topo_id,eta_id
    REAL(MK),     DIMENSION(ppm_dim)           :: coeffs
    INTEGER,      DIMENSION(ppm_dim)           :: order
    INTEGER,      DIMENSION(ppm_dim*ppm_dim)   :: degree
    REAL(MK),     DIMENSION(ppm_dim)           :: wp_grad_fun0,wp_grad_fun_proj,wp_dir, wp_dir_temp, vec
    REAL(MK)                                   :: dummy_wp

    !-------------------------------------------------------------------------!
    ! Initialize
    !-------------------------------------------------------------------------!
    info = 0

#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif

    dummy_grad = 0._MK
    dummy_wp = 0._MK

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

    !check that we are dealing with anisotropic particles 
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


    
    ! HAECKIC: here a wpv is allocated
    !          Isn't that done in adapt_particles.f??
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
                    ppm_dim,info,with_ghosts=.TRUE.,name='adapt_wpgrad')
                IF (info.NE.0) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_alloc,caller,&
                        'particles_allocate_wpv failed', __LINE__,info)
                    GOTO 9999
                ENDIF
            ELSE
                IF (.NOT.Particles%wpv(adapt_wpgradid)%is_mapped) THEN
                    CALL particles_allocate_wpv(Particles,adapt_wpgradid,&
                        ppm_dim,info,with_ghosts=.TRUE.,&
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


    
    ! HAECKIC: nothing adapted here
    ! crash if not enough neighbours
    ! adapt the plotting
    IF (need_derivatives .AND. Particles%nneighmin.LT.opts%nneigh_critical) THEN
!         xp => Get_xp(Particles)
!         rcp => Get_wps(Particles,Particles%rcp_id)
!         WRITE(cbuf,*) 'Not enough neighbours'
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         WRITE(cbuf,'(2(A,I5,2X))') 'nneigh_critical ',opts%nneigh_critical,&
!             'nneigh_toobig ',opts%nneigh_toobig
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         WRITE(cbuf,'(A,I5)') 'We have nneighmin = ',Particles%nneighmin
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         WRITE(cbuf,*) 'Writing debug data to file fort.1230+rank'
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         WRITE(myformat,'(A,I1,A)') '(',ppm_dim+1,'(E30.16,2X),I4)'
!         DO ip=1,Particles%Npart
!             WRITE(1230+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
!         ENDDO
!         DO ip=Particles%Npart+1,Particles%Mpart
!             WRITE(1250+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),0
!         ENDDO
!         CALL ppm_write(ppm_rank,caller,&
!             'Calling neighlists anyway, but crashing just after that',info)
!         CALL particles_neighlists(Particles,topo_id,info)
!         IF (info .NE. 0) THEN
!             CALL ppm_write(ppm_rank,caller,'particles_neighlists failed.',info)
!             info = -1
!             GOTO 9999
!         ENDIF
!         WRITE(cbuf,'(A,I5,1X,I5)') 'nneighmin is now: ',Particles%nneighmin
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         WRITE(cbuf,*) 'Writing debug data to file fort.1240+rank'
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         DO ip=1,Particles%Npart
!             WRITE(1240+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),Particles%nvlist(ip)
!         ENDDO
!         DO ip=Particles%Npart+1,Particles%Mpart
!             WRITE(1260+ppm_rank,myformat) xp(1:ppm_dim,ip),rcp(ip),0
!         ENDDO
!         xp=>NULL()
!         rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
!         WRITE(cbuf,'(2(A,I6),2(A,E30.20))') 'Npart=',&
!             Particles%Npart,' Mpart=',Particles%Mpart,&
!             ' cutoff=',Particles%cutoff,' skin=',Particles%skin
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!         info = -1
!         GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Compute D (desired resolution)
    !!-------------------------------------------------------------------------!
    ! (re)allocate Dtilde, which is a tensor

    ! HAECKIC: here a wpv is allocated
    CALL particles_allocate_wpv(Particles,Particles%Dtilde_id,Particles%tensor_length,info,&
        with_ghosts=.TRUE.,iopt=ppm_param_alloc_grow,name='D_tilde')
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'particles_allocate_wpv failed', __LINE__,info)
        GOTO 9999
    ENDIF


    IF (.NOT.PRESENT(wp_fun)) THEN
      
    ! HAECKIC: what if D_fun depends on wp_fun and we don not have it

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

         ! HAECKIC: anisotropic gradient
         CALL particles_get_grad_aniso(Particles,adapt_wpgradid,info,with_ghosts=.TRUE.)
         IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                  'particles_get_grad_aniso failed', __LINE__,info)
            GOTO 9999
         ENDIF

         
         ! HAECKIC: this is a temporary test for plane example -----------
         wp_grad => Get_wpv(Particles,adapt_wpgradid)
         wp => Get_wps(Particles,Particles%adapt_wpid)
         xp => Get_xp(Particles)
         old_scale = 0.0_mk
         new_scale = 0.0_mk
         max_g = 0.0_mk
         max_w = 0.0_mk
         DO ip=1,Particles%Npart
            wp_grad_fun0 = f0_grad_fun(xp(1:ppm_dim,ip))
            old_scale = old_scale+SUM((wp_grad(1:ppm_dim,ip)-wp_grad_fun0)**2)
            dummy_wp = f0_fun(xp(1:ppm_dim,ip))
            new_scale = new_scale+((wp(ip)-dummy_wp)**2)
            IF (max_g .LT. SUM((wp_grad(1:ppm_dim,ip)-wp_grad_fun0)**2)) THEN
               max_g = SUM((wp_grad(1:ppm_dim,ip)-wp_grad_fun0)**2)
            ENDIF
            IF (max_w .LT. (wp(ip)-dummy_wp)**2) THEN
               max_w = ((wp(ip)-dummy_wp)**2)
            ENDIF
         ENDDO
         write(*,*) '----Errors----'
         write(*,*) 'Squared Error per particle in gradient: ', old_scale/Particles%Npart
         write(*,*) 'Max Error per particle in gradient: ', max_g
         write(*,*) 'Squared Error per particle in field: ', new_scale/Particles%Npart
         write(*,*) 'Max Error per particle in field: ', max_w
         wp_grad => Set_wpv(Particles,adapt_wpgradid,read_only=.TRUE.)
         wp => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
         xp => Set_xp(Particles,read_only=.TRUE.)
         ! ------------------------------------------------------------------

    ENDIF if_needs_derivatives
     
    ! HAECKIC: Entire next part is changed

    !-------------------------------------------------------------------------!
    ! Compute D_tilde
    !-------------------------------------------------------------------------!
    if_D_needs_grad: IF (opts%D_needs_gradients) THEN

        xp => Get_xp(Particles)
        Dtilde => Get_wpv(Particles,Particles%Dtilde_id)

        IF (PRESENT(wp_grad_fun)) THEN

        ELSE
            wp_grad => Get_wpv(Particles,adapt_wpgradid,with_ghosts=.TRUE.)

             IF (PRESENT(wp_fun)) THEN
             
             ELSE
               ! get wp approx
             ENDIF

        ENDIF

        ! For each particle calc the gradient using given functions wp_grad_fun
        DO ip=1,Particles%Npart
            ! gradient
            IF (PRESENT(wp_grad_fun)) THEN
                  wp_grad_fun0 = wp_grad_fun(xp(1:ppm_dim,ip))
                  old_scale = SQRT(SUM(wp_grad_fun0**2))
                  new_scale = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun0,opts)/old_scale
            ELSE
                  IF (PRESENT(wp_fun)) THEN
                     wp_grad_fun0 = wp_grad(1:ppm_dim,ip)
                     old_scale = SQRT(SUM(wp_grad_fun0**2))
                     new_scale = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun0,opts)/old_scale
                  ELSE
                     wp_grad_fun0 = wp_grad(1:ppm_dim,ip)
                     old_scale = SQRT(SUM(wp_grad_fun0**2))
                     new_scale = D_fun(dummy_wp,wp_grad_fun0,opts)/old_scale
                  ENDIF
            ENDIF

            ! HAECKIC: TODO: determine orthogonal direction by projection of neighbors in this direction
            ! 1. project gradient of neighbors and apply d_fun to make of it
            ! 2. project smaller axes: 0 -> max, else: relative take antiproportional to smaller axes of neigh

            IF (.NOT. Particles%neighlists) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,caller,&
                  'need neighbour lists to be uptodate', __LINE__,info)
               GOTO 9999
            ENDIF
            
            ! use a dummy wp
            ! get the maximum gradient towards orthogonal direction
            wp_grad_fun_proj = (/0.0_mk,0.0_mk/)
            wp_dir = (/-wp_grad_fun0(2),wp_grad_fun0(1)/)
            DO ineigh=1,Particles%nvlist(ip)
                  iq = Particles%vlist(ineigh,ip)
                  IF (PRESENT(wp_grad_fun)) THEN
                     wp_dir_temp = (wp_dir*wp_grad_fun(xp(1:ppm_dim,iq)))/SQRT(SUM(wp_dir**2))
                  ELSE
                     wp_dir_temp = (wp_dir*wp_grad(1:ppm_dim,iq))/SQRT(SUM(wp_dir**2))
                  ENDIF
                  IF (SQRT(SUM(wp_grad_fun_proj**2)) .GT. SQRT(SUM(wp_dir_temp**2))) THEN
                     wp_grad_fun_proj = wp_dir_temp
                  ENDIF
            ENDDO
      
            IF (PRESENT(wp_grad_fun)) THEN
                  orth_len = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun_proj,opts)/old_scale
            ELSE
                  IF (PRESENT(wp_fun)) THEN
                     orth_len = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun_proj,opts)/old_scale
                  ELSE
                     orth_len = D_fun(dummy_wp,wp_grad_fun_proj,opts)/old_scale
                  ENDIF
            ENDIF
           
            ! set the tensor
            CALL ppm_alloc(Matrix_A,(/ Particles%tensor_length /),ppm_param_alloc_fit,info)

            IF (ppm_dim .EQ. 2) THEN
            
               Matrix_A(1) = -orth_len*wp_grad_fun0(2)
               Matrix_A(3) =  orth_len*wp_grad_fun0(1)
               Matrix_A(2) = new_scale*wp_grad_fun0(1)
               Matrix_A(4) = new_scale*wp_grad_fun0(2)

               ! Check minimum length of axes
               vec = (/Matrix_A(1) , Matrix_A(3)/)
               IF (SQRT(SUM(vec**2)) .LT. opts%minimum_D) THEN
                  Matrix_A(1) = opts%minimum_D*Matrix_A(1)/SQRT(SUM(vec**2))
                  Matrix_A(3) = opts%minimum_D*Matrix_A(3)/SQRT(SUM(vec**2))
               ENDIF
               
               vec = (/Matrix_A(2) , Matrix_A(4)/)
               IF (SQRT(SUM(vec**2)) .LT. opts%minimum_D) THEN
                  Matrix_A(2) = opts%minimum_D*Matrix_A(2)/SQRT(SUM(vec**2))
                  Matrix_A(4) = opts%minimum_D*Matrix_A(4)/SQRT(SUM(vec**2))
               ENDIF
               
               ! Check maximum length of axes
               vec = (/Matrix_A(1) , Matrix_A(3)/)
               IF (SQRT(SUM(vec**2)) .GT. opts%maximum_D) THEN
                  Matrix_A(1) = opts%maximum_D*Matrix_A(1)/SQRT(SUM(vec**2))
                  Matrix_A(3) = opts%maximum_D*Matrix_A(3)/SQRT(SUM(vec**2))
               ENDIF
               
               vec = (/Matrix_A(2) , Matrix_A(4)/)
               IF (SQRT(SUM(vec**2)) .GT. opts%maximum_D) THEN
                  Matrix_A(2) = opts%maximum_D*Matrix_A(2)/SQRT(SUM(vec**2))
                  Matrix_A(4) = opts%maximum_D*Matrix_A(4)/SQRT(SUM(vec**2))
               ENDIF

               CALL particles_inverse_matrix(Matrix_A,Matrix_B,info)

               Dtilde(1,ip) = Matrix_B(1)
               Dtilde(2,ip) = Matrix_B(2)
               Dtilde(3,ip) = Matrix_B(3)
               Dtilde(4,ip) = Matrix_B(4)

            ELSE
               
               ! HAECKIC: TODO: 3D Dtilde computation
               ! Take any direction in plane and use projection to scale them
   
            ENDIF
               
        ENDDO
        
        IF (PRESENT(wp_grad_fun)) THEN

        ELSE
             wp_grad => Set_wpv(Particles,adapt_wpgradid,read_only=.TRUE.)

             IF (PRESENT(wp_fun)) THEN
             
             ELSE
                  ! set wp
             ENDIF

        ENDIF
        
        Dtilde => Set_wpv(Particles,Particles%Dtilde_id)
        xp => Set_xp(Particles,read_only=.TRUE.)

        ! Get ghosts for D_tilde
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_mapping_ghosts failed',info)
            info = -1
            GOTO 9999
        ENDIF

    ELSE ! .NOT. D_needs_grad
         
         ! HAECKIC: not treated case in anisotropic set up
         IF (PRESENT(wp_fun)) THEN

         ELSE
         
         ENDIF

    ENDIF if_D_needs_grad

    !-------------------------------------------------------------------------!
    ! Rescale D_tilde
    !-------------------------------------------------------------------------!
    ! minimum length of gradient axes, maximum length of gradient axis
    ! HAECKIC: DROPPED
!     Dtilde => Get_wpv(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
!     DO ip = 1,Particles%Mpart
!         CALL particles_inverse_matrix(Dtilde(:,ip),Matrix_B,info)
!         old_scale = sqrt(Matrix_B(2)*Matrix_B(2) + Matrix_B(4)*Matrix_B(4))
!         new_scale = MIN(MAX(old_scale,opts%minimum_D),opts%maximum_D)
!         IF (.NOT.(new_scale .EQ. old_scale)) THEN
!             ! inverse scaling
!             new_scale = old_scale/new_scale
!             Dtilde(1:Particles%tensor_length,ip) = new_scale * Dtilde(1:Particles%tensor_length,ip)
!         ENDIF
!     ENDDO
!     Dtilde => Set_wpv(Particles,Particles%Dtilde_id,&
!         read_only=.FALSE.,ghosts_ok=.TRUE.)

    !---------------------------------------------------------------------!
    ! Update the real tensors of the particles using 1/rcp_over_D*Dtilde
    !---------------------------------------------------------------------!
    inv => Get_wpv(Particles,Particles%G_id)
    Dtilde => Get_wpv(Particles,Particles%Dtilde_id)
    DO ip=1,Particles%Npart
        ! inverse scaling
        new_scale = 1/(opts%rcp_over_D)
        inv(1:Particles%tensor_length,ip) = new_scale * Dtilde(1:Particles%tensor_length,ip)

    ENDDO
    inv => Set_wpv(Particles,Particles%G_id)
    Dtilde => Set_wpv(Particles,Particles%Dtilde_id,read_only=.TRUE.)

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

    ! HAECKIC: here a wpv is allocated
    IF (Particles%D_id .EQ. 0 ) THEN
        CALL particles_allocate_wpv(Particles,Particles%D_id,Particles%tensor_length,&
            info,name='D')
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_alloc,caller,&
                'particles_allocate_wps failed',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF


    ! HAECKIC: completely different
    !------------------------------------------------------------------------------!
    ! D is Dtilde scaled with minium length of shorter axis over all neighbours iq
    !------------------------------------------------------------------------------!
    D      => Get_wpv(Particles,Particles%D_id)
    Dtilde => Get_wpv(Particles,Particles%Dtilde_id,with_ghosts=.TRUE.)
    DO ip=1,Particles%Npart
        ! Copy the tensor and get smallest in neighborhood
        ! HAECKIC: Scaling of axes can be differently
        ! we acutally should project parts of it on different axes
        D(1:Particles%tensor_length,ip) = Dtilde(1:Particles%tensor_length,ip)
        CALL particles_shorter_axis(Particles,ip,old_scale)
        new_scale = old_scale
        CALL particles_longer_axis(Particles,ip,new_scale_long)
        DO ineigh=1,Particles%nvlist(ip)
            iq = Particles%vlist(ineigh,ip)
            CALL particles_shorter_axis(Particles,iq,temp_scale) 
            IF (temp_scale.LT.new_scale) THEN
                new_scale = temp_scale
                CALL particles_longer_axis(Particles,iq,new_scale_long)
            ENDIF
        ENDDO
        ! scale the axes of the ellipse
        ! inverse scaling
        new_scale = old_scale/new_scale
        CALL particles_longer_axis(Particles,ip,old_scale)
        new_scale_long = old_scale/new_scale_long
        IF (ppm_dim.eq.2) THEN
            ! HAECKIC: inverse scaling with vector
            D(1,ip) = new_scale_long*D(1,ip)
            D(2,ip) = new_scale_long*D(2,ip)
            D(3,ip) = new_scale*D(3,ip)
            D(4,ip) = new_scale*D(4,ip)
        ELSE
            !HAECKIC: TODO 3D case
        ENDIF
    ENDDO
    D      => Set_wpv(Particles,Particles%D_id)
    Dtilde => Set_wpv(Particles,Particles%Dtilde_id,read_only=.TRUE.)


! #if debug_verbosity > 0
!     D => Get_wpv(Particles,Particles%D_id)
! #ifdef __MPI
!     CALL MPI_Allreduce(MINVAL(D(1:Particles%Npart)),min_D,1,&
!         ppm_mpi_kind,MPI_MIN,ppm_comm,info)
! #else
!     min_D =MINVAL(D(1:Particles%Npart))
! #endif
!     IF (ppm_rank .EQ.0) THEN
!         WRITE(cbuf,'(A,E12.4)') 'Min D = ',min_D
!         CALL ppm_write(ppm_rank,caller,cbuf,info)
!     ENDIF
!     D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
! #endif

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error


END SUBROUTINE sop_compute_D

#undef __KIND
