!!!----------------------------------------------------------------------------!
!!!
!!! Do the gradient descent on xp, using interpolation from {xp_old} to
!!! {xp} at each step (requires cross--neighbour-list)
!!! All particle properties (D, rcp, etc...) are assumed to be referring
!!! to xp
!!! 
!!! On input: vlist is up-to-date
!!! On output: vlist and vlist_cross have enough elements to 
!!!            allow for computing corrected kernels)
!!! 
!!!----------------------------------------------------------------------------!

SUBROUTINE sop_gradient_descent(Particles_old,Particles, &
        nvlist_cross,vlist_cross,    &
        nneighmin_cross,nneighmax_cross,num_it,opts,info, &
        wp_fun,D_fun,wp_grad_fun,threshold,need_deriv,stats)

    !USE ppm_module_inl_xset_vlist
    USE ppm_module_io_vtk

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
    TYPE(ppm_t_particles), POINTER,        INTENT(IN   )   :: Particles_old
    TYPE(ppm_t_particles), POINTER,        INTENT(INOUT)   :: Particles
    INTEGER,      DIMENSION(:),  POINTER,  INTENT(INOUT)   :: nvlist_cross
    INTEGER,      DIMENSION(:,:),POINTER,  INTENT(INOUT)   :: vlist_cross
    INTEGER,                               INTENT(INOUT)   :: nneighmax_cross
    INTEGER,                               INTENT(INOUT)   :: nneighmin_cross
    INTEGER,                               INTENT(  OUT)   :: num_it
    TYPE(sop_t_opts), POINTER,             INTENT(IN   )   :: opts
    INTEGER,                               INTENT(  OUT)   :: info

    !optional arguments
    REAL(MK), OPTIONAL,                    INTENT(IN)      :: threshold
    LOGICAL,  OPTIONAL,                    INTENT(IN)      :: need_deriv
    !Monitor function
    OPTIONAL                                               :: D_fun
    !Field function (usually known only during initialisation)
    OPTIONAL                                               :: wp_fun
    !Gradient of the field function (usually known only during initialisation)
    OPTIONAL                                               :: wp_grad_fun
    TYPE(sop_t_stats),  POINTER,OPTIONAL,  INTENT(  OUT)  :: stats
    !!! statistics on output
    ! argument-functions need an interface
    INTERFACE
        FUNCTION D_fun(f,dfdx,opts,lap_f)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
            USE ppm_module_sop_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK)                                     :: D_fun
            REAL(MK),                         INTENT(IN) :: f
            REAL(MK),DIMENSION(ppm_dim),      INTENT(IN) :: dfdx
            TYPE(sop_t_opts),POINTER,         INTENT(IN) :: opts
            REAL(MK),OPTIONAL,                INTENT(IN) :: lap_f
        END FUNCTION D_fun

        FUNCTION wp_grad_fun(pos)
            USE ppm_module_data, ONLY: ppm_dim
            USE ppm_module_typedef
#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
            REAL(MK),DIMENSION(ppm_dim)             :: wp_grad_fun
        END FUNCTION wp_grad_fun
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
    END INTERFACE

    ! local variables
    INTEGER                             :: ip, it_adapt,iq,ineigh,di,j,k
    INTEGER                             :: iunit
    CHARACTER(LEN = 256)                :: filename,cbuf
    CHARACTER(LEN = 256)                :: caller='sop_gradient_descent'
    REAL(KIND(1.D0))                    :: t0

    REAL(MK)                            :: step,step_previous,alpha1,alpha2
    REAL(MK)                            :: step_min,step_max,step_stall
    INTEGER                             :: it_adapt_max
    REAL(MK)                            :: Psi_max,Psi_global,Psi_global_old
    REAL(MK)                            :: Psi_1,Psi_2,Psi_thresh
    REAL(MK),DIMENSION(ppm_dim)         :: dist,dist2,dummy_grad
    REAL(MK),DIMENSION(:,:),POINTER     :: Gradient_Psi => NULL()
    INTEGER                             :: nneigh_adapt

    REAL(MK),DIMENSION(:,:),POINTER     :: xp => NULL()
    REAL(MK),DIMENSION(:,:),POINTER     :: D => NULL()
    REAL(MK),DIMENSION(:,:),POINTER     :: inv => NULL()
    REAL(MK),DIMENSION(:,:),POINTER     :: Dtilde => NULL()
    REAL(MK),DIMENSION(:,:),POINTER     :: xp_old => NULL()
    REAL(MK),DIMENSION(:),POINTER       :: D_old => NULL()

    REAL(MK)                            :: tmpvar1,tmpvar2,minDold, dist1s, dist2s, dists, min_dist
    REAL(MK)                            :: weight,weight_sum
    REAL(MK)                            :: almostzero
    INTEGER                             :: tmpvari1,tmpvari2
    LOGICAL                             :: need_derivatives
    INTEGER                             :: topo_id
    LOGICAL                             :: adding_particles
    INTEGER                             :: nb_spawn
    
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_A => NULL()
    REAL(MK),     DIMENSION(:),   POINTER      :: Matrix_B => NULL()

    !should be removed once the argument lists for the inl routines
    !have been updated to inhomogeneous ghostlayers
    REAL(MK),DIMENSION(2*ppm_dim)       :: ghostlayer


    !!-------------------------------------------------------------------------!
    !! Initialize
    !!-----------------------------2nd--------------------------------------------!
    info = 0
#if debug_verbosity > 0
    CALL substart(caller,t0,info)
#endif
    dummy_grad=0._MK
    almostzero=EPSILON(1._MK)
    topo_id = Particles%active_topoid
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

    IF (.NOT.Particles%neighlists) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'compute neighbour lists before',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

    IF (.NOT.Particles%has_ghosts) THEN
        info = ppm_error_fatal
        CALL ppm_error(ppm_err_alloc,caller,   &
            &  'Need ghosts particles to be updated on entry',&
            &  __LINE__,info)
        GOTO 9999
    ENDIF

#if debug_verbosity > 0
    IF(PRESENT(wp_grad_fun)) THEN
        WRITE(cbuf,'(A)') 'Using analytical expressions to compute D'
        IF (.NOT. PRESENT(D_fun)) THEN
            WRITE(cbuf,'(A)') 'Incompatible options. D_fun must be present'
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF
    ELSE
        WRITE(cbuf,'(A)') 'Using interpolation routines to compute D'
    ENDIF
    IF (ppm_rank .EQ.0) THEN
        CALL ppm_write(ppm_rank,caller,cbuf,info)
    ENDIF
#endif
    IF (PRESENT(need_deriv)) THEN
        need_derivatives=need_deriv
    ELSE
        need_derivatives=.TRUE.
    ENDIF

    ! number of steps
    it_adapt = 0
    ! New particles added?
    adding_particles=.TRUE.

    ! some potential stuff needed
    Psi_max = HUGE(1._MK)
    Psi_global = HUGE(1._MK)
    Psi_global_old = HUGE(1._MK)

    ! Adaptivity stopping criterion
    IF (PRESENT(threshold)) THEN
        Psi_thresh = threshold
    ELSE
        Psi_thresh = 2.5_MK
    ENDIF

    IF (MK .EQ. KIND(1.D0)) THEN
        step_min = 1D-8; step_stall = 1D-14
        !step_min = 1D-3; step_stall = 1D-14
    ELSE
        step_min = 1E-8; step_stall = 1E-14
        !step_min = 1E-3; step_stall = 1E-14
    ENDIF
    step_max = 0.1_MK ! 0.1_MK
    it_adapt_max = 1000

    step = 1._MK
    nneigh_adapt = opts%nneigh_theo

    7099 CONTINUE

    IF (ASSOCIATED(Gradient_Psi)) DEALLOCATE(Gradient_Psi)
    ALLOCATE(Gradient_Psi(ppm_dim,Particles%Mpart),STAT=info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'allocation failed',info)
        info = -1
        GOTO 9999
    ENDIF

    !!-------------------------------------------------------------------------!
    !! Gradient descent loop until stopping criterion is met
    !!-------------------------------------------------------------------------!
    !! until one of the following is violated:
    !! 1. step > step_stall, i.e step is so small we do not change anything
    !! 2. it_adapt < adapt_max, i.e. until max number of iteration reached
    !! 3. Psi_max > Psi_thresh OR new particles added, i.e. max potential is below threshold
    it_adapt_loop: DO WHILE (step .GT. step_stall .AND. &
            &          it_adapt .LT. it_adapt_max .AND. &
            &     ((Psi_max .GT. Psi_thresh)) .OR. adding_particles)
    
        it_adapt = it_adapt + 1

        write(*,*) 'step: ', it_adapt, it_adapt_max

        !!---------------------------------------------------------------------!
        !! Update number of particles by fusion/insertion
        !! of particles (need vlist to be up-to-date)
        !! (on output, Mpart becomes meaningless.)
        !!---------------------------------------------------------------------!
        !NOTE: if we know D_tilde analytically, we should not supply
        ! new particles with a value for D (i.e. not allocate the array yet)
        ! do fuse/spawn first, then apply_bc, get_ghosts, etc, then alloc and 
        ! compute D

        !FIXME: in theory, we should do apply_bc+remap+get_ghosts here
        ! (Fusing particles requires knowing the ghosts and particles
        ! have moved after the linesearch in the previous iteration)
        ! There MAY be a better way of doing this...
        !!---------------------------------------------------------------------!

        write(*,*) Particles%Npart, Particles%Mpart

        CALL particles_apply_bc(Particles,topo_id,info)
        CALL particles_mapping_partial(Particles,topo_id,info)
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        CALL particles_neighlists(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_neighlists failed',__LINE__,info)
            GOTO 9999
        ENDIF


!         !call check_duplicates(Particles)


        !Delete (fuse) particles that are too close to each other
        !(needs ghost particles to be up-to-date)

        CALL sop_fuse_particles(Particles,opts,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'sop_fuse_particles failed',__LINE__,info)
            GOTO 9999
        ENDIF

!          min_dist = 1.0_mk
!          DO ip=1,Particles%Npart
!             DO j=1,Particles%nvlist(ip)
!                iq = Particles%vlist(j,ip)
!                CALL particles_anisotropic_distance(Particles,ip,iq,dist1s,info)
!                CALL particles_anisotropic_distance(Particles,iq,ip,dist2s,info)
!                dists = MIN(dist1s,dist2s)
!                IF (dists .LT. min_dist) THEN
!                   min_dist = dists
!                   k = ip
!                   di = iq
!                ENDIF
! 
!             ENDDO
!          ENDDO
!          write(*,*) 'before 1', min_dist, k, di

        !we only removed particles, but they didnt move.
        Particles%areinside=.TRUE.
        Particles%ontopology=.TRUE.
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'mapping ghost failed',__LINE__,info)
            GOTO 9999
        ENDIF
        CALL particles_neighlists(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'neighlists failed here',__LINE__,info)
            GOTO 9999
        ENDIF

!          min_dist = 1.0_mk
!          DO ip=1,Particles%Npart
!             DO j=1,Particles%nvlist(ip)
!                iq = Particles%vlist(j,ip)
!                CALL particles_anisotropic_distance(Particles,ip,iq,dist1s,info)
!                CALL particles_anisotropic_distance(Particles,iq,ip,dist2s,info)
!                dists = MIN(dist1s,dist2s)
!                IF (dists .LT. min_dist) THEN
!                   min_dist = dists
!                ENDIF
! 
!             ENDDO
!          ENDDO
!          write(*,*) 'after ', min_dist


        !call check_duplicates(Particles)

        !if(it_adapt.LE.50)then
        !Insert (spawn) new particles where needed

        !write(*,*) 'before: ', Particles%Npart, Particles%Mpart
         min_dist = 1.0_mk
         DO ip=1,Particles%Npart
            DO j=1,Particles%nvlist(ip)
               iq = Particles%vlist(j,ip)

                     

               CALL particles_anisotropic_distance(Particles,ip,iq,dist1s,info)
               CALL particles_anisotropic_distance(Particles,iq,ip,dist2s,info)
               
               IF (ip.eq.216 .and. iq .eq. 243 ) THEN
                  write(*,*) 'disti ',ip, dist1s-dist2s, Particles%xp(:,ip)
               ENDIF
               IF (ip.eq.243 .and. iq .eq. 216 ) THEN
                  write(*,*) 'disti ',ip, dist1s-dist2s, Particles%xp(:,ip)
               ENDIF
               
               dists = MIN(dist1s,dist2s)
               IF (dists .LT. min_dist) THEN
                  min_dist = dists
                  k = ip
                  di = iq
               ENDIF

            ENDDO
         ENDDO
         write(*,*) 'before ', min_dist, k, di


        IF (min_dist+1d-12 .LT. opts%fuse_radius) THEN
            write(*,*) 'EEEEEERRRRRRRRRROOOOOOOOOOOORRRRRRRRRRRR'
         ENDIF

        CALL  sop_spawn_particles(Particles,opts,info,&
            nb_part_added=nb_spawn,wp_fun=wp_fun)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'sop_spawn_particles failed',__LINE__,info)
            GOTO 9999
        ENDIF

        !write(*,*) 'after: ', Particles%Npart, Particles%Mpart
        ! haeckic: to be changed
        adding_particles = nb_spawn .GT. 0
        !endif

        !call check_duplicates(Particles)
        CALL particles_updated_positions(Particles,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_updated_positions failed',__LINE__,info)
            GOTO 9999
        ENDIF


!          DO ip=1,Particles%Npart
!             write(*,*) ip,Particles%wpv(Particles%G_id)%vec(1:4,ip)
!          ENDDO

        !!---------------------------------------------------------------------!
        !! Ensure that particles satisfy the boundary conditions
        !! (Here because partial remapping fails otherwise).
        !!---------------------------------------------------------------------!
        CALL particles_apply_bc(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_apply_bc failed',__LINE__,info)
            GOTO 9999
        ENDIF

        CALL particles_mapping_partial(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_mapping_partial failed',__LINE__,info)
            GOTO 9999
        ENDIF

        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'particles_mapping_ghosts failed',__LINE__,info)
            GOTO 9999
        ENDIF


!         Compute_D: IF (PRESENT(wp_grad_fun).OR. &
!             (.NOT.need_derivatives.AND.PRESENT(wp_fun))) THEN
!             !!-----------------------------------------------------------------!
!             !! Get anisotropic requirements directly using the given funcs
!             !!-----------------------------------------------------------------!
!             CALL sop_compute_D(Particles,D_fun,opts,info,     &
!                 wp_fun=wp_fun,wp_grad_fun=wp_grad_fun)
!             IF (info .NE. 0) THEN
!                 info = ppm_error_error
!                 CALL ppm_error(ppm_err_sub_failed,caller,&
!                     'sop_compute_D failed',__LINE__,info)
!                 GOTO 9999
!             ENDIF
!         ELSE
!             ! do not update D, we make that using the xset and the D of the old particles
! 
!         ENDIF Compute_D
! 
! 
!         IF (.NOT. PRESENT(wp_fun)) THEN
!             !------------------------------------------------------------------!
!             ! Update D and cutoff radii
!             !------------------------------------------------------------------!
! 
!             ! The following is used to prevent particles with small D from
!             ! drifting away inside the domain, eventually generating plenty
!             ! of other small particles that can potentially fill the whole
!             ! box (defeating the point of having an adaptive scheme...)
!             ! Roughly, it means that a particle can have a small D only
!             ! if it has neighbours from the older generation (D_old) that
!             ! also have a small D.
! !             D_old => Get_wps(Particles_old,Particles_old%Dtilde_id)
! !             D_old = D_old * opts%rcp_over_D
! !             ghostlayer=Particles%cutoff
! !             CALL ppm_inl_xset_vlist(topo_id,Particles%xp,Particles%Npart,&
! !                 Particles%Mpart,Particles_old%xp,Particles_old%Npart,&
! !                 Particles_old%Npart,D_old,Particles%skin,&
! !                 ghostlayer,info,vlist_cross,nvlist_cross)
! !             D_old = D_old / opts%rcp_over_D
! !             D_old => Set_wps(Particles_old,Particles_old%Dtilde_id,&
! !                 read_only=.TRUE.)
! !             IF (info .NE. 0) THEN
! !                 CALL ppm_write(ppm_rank,caller,&
! !                     'ppm_inl_xset_vlist failed.',info)
! !                 info = -1
! !                 GOTO 9999
! !             ENDIF
! ! 
! ! #if debug_verbosity > 0
! !             IF (MINVAL(nvlist_cross(1:Particles%Npart)).LE.0) THEN
! !                 CALL ppm_write(ppm_rank,caller,&
! !                     'Insufficient number of xset neighbours to compute D',info)
! !                 info = -1
! !                 GOTO 9999
! !             ENDIF
! ! #endif
! ! 
! !             D     => Get_wps(Particles,    Particles%D_id)
! !             D_old => Get_wps(Particles_old,Particles_old%Dtilde_id)
! !             rcp => Get_wps(Particles,Particles%rcp_id)
! !             DO ip=1,Particles%Npart
! !                 minDold=HUGE(1._MK)
! !                 DO ineigh=1,nvlist_cross(ip)
! !                     iq=vlist_cross(ineigh,ip)
! !                     minDold=MIN(minDold,D_old(iq))
! !                 ENDDO
! !                 !D(ip) = MAX(D(ip),minDold)
! !                 !when Dtilde/D is large, increase rcp
! !                 rcp(ip) = opts%rcp_over_D * MIN(D(ip),2._MK*minDold)
! !                 D(ip) = minDold
! !             ENDDO
! !             D     => Set_wps(Particles,    Particles%D_id)
! !             D_old => Set_wps(Particles_old,Particles_old%Dtilde_id,&
! !                 read_only=.TRUE.)
! !             rcp => Set_wps(Particles,Particles%rcp_id)
!         ELSE
!             !------------------------------------------------------------------!
!             ! Update cutoff radii
!             !------------------------------------------------------------------!
!             D => Get_wpv(Particles,Particles%D_id)
!             Dtilde => Get_wpv(Particles,Particles%Dtilde_id)
!             inv => Get_wpv(Particles,Particles%G_id)
! 
!             DO ip=1,Particles%Npart
!                 ! haeckic: do the right thing here 
!                 !Update inv using Dtilde computed before, do the min?
!                 !rcp(ip) = opts%rcp_over_D * MIN(Dtilde(ip),2._MK*D(ip))
!                 IF (ppm_dim .EQ. 2) THEN
! !                      Matrix_A(1) = 0
! !                      Matrix_A(3) = 0
! !                      Matrix_A(2) = 0
! !                      Matrix_A(4) = 0
! ! 
! !                      CALL particles_inverse_matrix(Matrix_A,Matrix_B,info)
! 
! !                      Dtilde(1,ip) = Matrix_B(1)
! !                      Dtilde(2,ip) = Matrix_B(2)
! !                      Dtilde(3,ip) = Matrix_B(3)
! !                      Dtilde(4,ip) = Matrix_B(4)
!                        inv(1,ip) = 1/(opts%rcp_over_D)*Dtilde(1,ip)
!                        inv(2,ip) = 1/(opts%rcp_over_D)*Dtilde(2,ip)
!                        inv(3,ip) = 1/(opts%rcp_over_D)*Dtilde(3,ip)
!                        inv(4,ip) = 1/(opts%rcp_over_D)*Dtilde(4,ip)
! 
!                        !haeckic: consider the neighbor seach for minimum
! 
! 
!                 ELSE
!                      ! haeckic: determine 3d orthogonal axes
!                      ! largest
! !                      Matrix_A(1) = -factor*new_scale*wp_grad_fun0(2)
! !                      Matrix_A(4) =  factor*new_scale*wp_grad_fun0(1)
! !                      Matrix_A(7) =  factor*new_scale*wp_grad_fun0(1)
! ! 
! !                      ! middle
! !                      Matrix_A(2) = -factor*new_scale*wp_grad_fun0(2)
! !                      Matrix_A(5) =  factor*new_scale*wp_grad_fun0(1)
! !                      Matrix_A(8) =  factor*new_scale*wp_grad_fun0(1)
! !                      
! !                      !smallest
! !                      Matrix_A(3) = new_scale*wp_grad_fun0(1)
! !                      Matrix_A(6) = new_scale*wp_grad_fun0(2)
! !                      Matrix_A(9) = new_scale*wp_grad_fun0(3)
! ! 
! !                      CALL particles_inverse_matrix(Matrix_A,Matrix_B,info)
! ! 
! !                      Dtilde(1,ip) = Matrix_B(1)
! !                      Dtilde(2,ip) = Matrix_B(2)
! !                      Dtilde(3,ip) = Matrix_B(3)
! !                      Dtilde(4,ip) = Matrix_B(4)
! !                      Dtilde(5,ip) = Matrix_B(5)
! !                      Dtilde(6,ip) = Matrix_B(6)
! !                      Dtilde(7,ip) = Matrix_B(7)
! !                      Dtilde(8,ip) = Matrix_B(8)
! !                      Dtilde(9,ip) = Matrix_B(9)
! !                      
!                 ENDIF
!             ENDDO
!             Dtilde => Set_wpv(Particles,Particles%Dtilde_id,read_only=.TRUE.)
!             D => Set_wpv(Particles,Particles%D_id,read_only=.TRUE.)
!             inv => Set_wpv(Particles,Particles%G_id)
!         ENDIF

#if debug_verbosity > 1
        WRITE(filename,'(A,I0,A,I0)') 'P_duringgraddesc_',&
            Particles%itime,'_',it_adapt
        CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif

        CALL particles_updated_cutoff(Particles,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,&
                'particles_updated_cutoff failed.',info)
            info = -1
            GOTO 9999
        ENDIF
        !---------------------------------------------------------------------!
        ! Update ghosts
        !---------------------------------------------------------------------!
        CALL particles_mapping_ghosts(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,&
                'particles_mapping_ghosts failed.',info)
            info = -1
            GOTO 9999
        ENDIF

        !---------------------------------------------------------------------!
        ! Update neighbour lists
        !---------------------------------------------------------------------!
        CALL particles_neighlists(Particles,topo_id,info)
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'particles_neighlists failed.',info)
            info = -1
            GOTO 9999
        ENDIF

        min_dist = 1.0_mk
         DO ip=1,Particles%Npart
            DO j=1,Particles%nvlist(ip)
               iq = Particles%vlist(j,ip)
               CALL particles_anisotropic_distance(Particles,ip,iq,dist1s,info)
               CALL particles_anisotropic_distance(Particles,iq,ip,dist2s,info)
               dists = MIN(dist1s,dist2s)
               IF (dists .LT. min_dist) THEN
                  min_dist = dists
                  k = ip
                  di = iq
               ENDIF

            ENDDO
         ENDDO
!          write(*,*) 'after ', min_dist, k, di, Particles%xp(:,k), Particles%xp(:,di)

#if debug_verbosity > 1
        CALL sop_dump_debug(Particles%xp,ppm_dim,Particles%Npart,&
            20000+it_adapt,info)
        CALL sop_dump_debug(Particles%wps(Particles%rcp_id)%vec,&
            Particles%Npart,30000+it_adapt,info)
        CALL sop_dump_debug(Particles%nvlist,Particles%Npart,40000+it_adapt,info)
        CALL sop_dump_debug(Particles%wps(Particles%D_id)%vec,&
            Particles%Npart,50000+it_adapt,info)
#endif
        !!---------------------------------------------------------------------!
        !! /begin Line search **
        !!---------------------------------------------------------------------!

        !!---------------------------------------------------------------------!
        !! Reallocate arrays whose sizes have changed
        !!---------------------------------------------------------------------!
        IF (SIZE(Gradient_Psi,2).LT.Particles%Mpart) THEN
            DEALLOCATE(Gradient_Psi)
            ALLOCATE(Gradient_Psi(ppm_dim,Particles%Mpart),STAT=info)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'allocation failed',info)
                info = -1
                GOTO 9999
            ENDIF
        ENDIF

        !!---------------------------------------------------------------------!
        !! Compute gradient of the potential
        !! (need ghosts for xp and D)
        !! (on output, the ghost values for Gradient_Psi have been updated)
        !!---------------------------------------------------------------------!
        CALL sop_gradient_psi(Particles,topo_id,Gradient_Psi,Psi_global,&
            Psi_max,opts,info) 
        IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,caller,'sop_gradient_psi failed.',info)
            info = -1
            GOTO 9999
        ENDIF

        write(*,*) 'PSI MAX', Psi_max, opts%fuse_radius
        !!---------------------------------------------------------------------!
        !! Writeout potential-vs-time to file
        !!---------------------------------------------------------------------!
        !IF (ppm_rank .EQ. 0) THEN
            !WRITE(filename,'(A,A)') TRIM(debugdir),'Psi_global.dat'
            !iunit=20
            !OPEN(UNIT=iunit,FILE=filename,FORM='FORMATTED',&
                !ACCESS='sequential',POSITION='APPEND',IOSTAT=info)
            !WRITE(iunit,'(3(E16.6,2X),I8,2X)') Psi_global, Psi_max,step,&
                !Particles%Npart
            !CLOSE(iunit)
        !ENDIF

        Psi_global_old = Psi_global
        Psi_1 = HUGE(1._MK)
        alpha1 = -1._MK
        alpha2 = -1._MK
        step_previous = 0._MK

        !!---------------------------------------------------------------------!
        !! Evaluate potential after different step sizes
        !!---------------------------------------------------------------------!
        linesearch_loop: DO WHILE (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK)
            xp => Get_xp(Particles,with_ghosts=.TRUE.)
            DO ip=1,Particles%Mpart
                xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + &
                    (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
            ENDDO
            xp => Set_xp(Particles,ghosts_ok=.TRUE.)
            step_previous = step

            !haeckic: why get D?
            D => Get_wpv(Particles,Particles%D_id,with_ghosts=.TRUE.)
            CALL sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)
            IF (info .NE. 0) THEN
                CALL ppm_write(ppm_rank,caller,'sop_potential_psi failed.',info)
                info = -1
                GOTO 9999
            ENDIF
            D => Set_wpv(Particles,Particles%D_id,read_only=.TRUE.)

            IF (Psi_global .LT. Psi_global_old) THEN 
                IF (Psi_global .LT. Psi_1) THEN
                    Psi_1 = Psi_global
                    alpha1 = step
                ENDIF
                step = 2._MK * step
                write(*,*) '2*step'
                IF (step .GT. step_max) THEN
                  write(*,*) 'EXIT: step > step_max'
                  EXIT linesearch_loop
                ENDIF
            ELSE IF (Psi_global .GT. Psi_global_old) THEN 
                Psi_2 = Psi_global
                alpha2 = step
                step = 0.5_MK * step
                write(*,*) '0.5*step'
                IF (step .LT. step_min) THEN
                  write(*,*) 'EXIT: step < step_min'
                  EXIT linesearch_loop
                ENDIF
            ELSE
                step = 10._MK * step
                IF (step .GT. step_max) THEN
                  write(*,*) 'EXIT: step > step_max'
                  EXIT linesearch_loop
                ENDIF
            ENDIF

            !FIXME
!             IF (ABS(Psi_global_old - Psi_global)/Psi_global_old .LT. 1E-4)THEN
!                write(*,*) 'EXIT: ABS'
!                EXIT linesearch_loop
!             ENDIF

        ENDDO linesearch_loop

        1000 CONTINUE

        IF (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK) THEN
            step = MAX(step_min,MIN(step,step_max))
            write(*,*) 'alpha < 0'
        ELSE
            !Quadratic fit
            Psi_1 = Psi_1 - Psi_global_old
            Psi_2 = Psi_2 - Psi_global_old
            step = 0.5_MK * (alpha1**2 * Psi_2 - alpha2**2 * Psi_1) / &
                (alpha1*Psi_2 - alpha2*Psi_1)
            write(*,*) 'quadratic fit'
        ENDIF

        !!---------------------------------------------------------------------!
        !! Choose best step size (from quadratic fit)
        !! Move particles (no need to move the ghosts, since we will have
        !! to get them through a local mapping anyway...)
        !!---------------------------------------------------------------------!
        !Move particles (including ghosts)
        write(*,*) 'Do a step', step
        write(*,*) 'Pot_now: ', Psi_global,' Pot_before: ', Psi_global_old
        xp => Get_xp(Particles)
        DO ip=1,Particles%Npart
            xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + &
            & (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
        ENDDO
        xp => Set_xp(Particles)

#if debug_verbosity > 2
#ifdef __MPI
        CALL MPI_Allreduce(step*MAXVAL(&
            ABS(Gradient_Psi(1:ppm_dim,1:Particles%Npart))),&
            tmpvar2,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
        IF (ppm_rank.EQ.0) THEN
            WRITE(cbuf,*) 'Moved particles: max displacement= ',tmpvar2
            CALL ppm_write(ppm_rank,caller,cbuf,info)
        ENDIF
#endif
#endif
        !!---------------------------------------------------------------------!
        !! /end Line search **
        !!---------------------------------------------------------------------!

#if debug_verbosity > 0
#ifdef __MPI
        CALL MPI_Allreduce(Psi_max,tmpvar2,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
        CALL MPI_Allreduce(Particles%Npart,tmpvari1,1,&
            MPI_INTEGER,MPI_SUM,ppm_comm,info)
        CALL MPI_Allreduce(Particles%Mpart,tmpvari2,1,&
            MPI_INTEGER,MPI_SUM,ppm_comm,info)
        IF (ppm_rank.EQ.0) THEN
            WRITE(cbuf,'(A,I3,2(A,E11.4),A,I4,1X,I4,1X,A,I6,A,I6,A,E7.2)') &
                'it_adapt= ',it_adapt,&
                ' Psi_mean= ',Psi_global/REAL(tmpvari1,MK),' Psi_max= ',tmpvar2, &
                ' Nneigh= ', Particles%nneighmin, Particles%nneighmax, &
                'Np=',tmpvari1,' Mp=',tmpvari2,' step=',step
            CALL ppm_write(ppm_rank,caller,cbuf,info)
        ENDIF
#endif
#endif
        CALL MPI_Allreduce(adding_particles,adding_particles,1,&
            MPI_LOGICAL,MPI_LOR,ppm_comm,info)

    ENDDO it_adapt_loop

    write(*,*) 'Iteration stopped: '
    write(*,*) 'stepsize: ', (step .GT. step_stall)
    write(*,*) 'step: ', (it_adapt .LT. it_adapt_max)
    write(*,*) 'Psi_max: ', (Psi_max .GT. Psi_thresh)
    write(*,*) 'adding: ', adding_particles

    !------------------------------------------------------------------
    ! Since particles have moved, we need to remap them
    !------------------------------------------------------------------
    CALL particles_apply_bc(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,&
            'particles_apply_bc failed',info)
        info = -1
        GOTO 9999
    ENDIF
    CALL particles_mapping_partial(Particles,topo_id,info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,&
            'particles_apply_bc failed',info)
        info = -1
        GOTO 9999
    ENDIF


#if debug_verbosity > 0
    WRITE(cbuf,'(2(A,E11.4))') 'Finished adapt loop. Psi_mean = ',&
        Psi_global,' Psi_max = ',Psi_max
    IF (ppm_rank.EQ.0) & 
        CALL ppm_write(ppm_rank,caller,cbuf,info)
#endif

    !returns number of iterations
    num_it = it_adapt 

    !!-------------------------------------------------------------------------!
    !! Finalize
    !!-------------------------------------------------------------------------!
    DEALLOCATE(Gradient_Psi)

#if debug_verbosity > 0
    CALL substop(caller,t0,info)
#endif

    9999 CONTINUE ! jump here upon error

END SUBROUTINE sop_gradient_descent
