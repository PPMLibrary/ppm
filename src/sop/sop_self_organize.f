      SUBROUTINE DTYPE(sop_self_organize)(this,D_fun,opts,info,Adapt_Field,    &
              wp_fun,wp_grad_fun,wplap_id,return_grad,level_fun,level_grad_fun,&
              smallest_sv,nb_fun,stats)
            !-------------------------------------------------------------------------
            !  Modules
            !-------------------------------------------------------------------------
            !!!----------------------------------------------------------------------!
            !!! This is the core routine of the adaptive particle method
            !!!
            !!! It adapts the particles using a monitor function D_fun that depends on
            !!! the current field wp and/or its gradient wp_grad
            !!! If wp and its gradient are not provided analytically,
            !!! they are approximated
            !!! using particle-particle interpolation with DC operators.
            !!! After adaptation, the properties carried by particles are interpolated
            !!! on the new positions.
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
            !!!----------------------------------------------------------------------!


            !-------------------------------------------------------------------------
            !  Modules
            !-------------------------------------------------------------------------
            USE ppm_module_map_part
            USE ppm_module_io_vtk

          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          CLASS(DTYPE(ppm_t_sop))                               :: this
          !!! particles
          TYPE(DTYPE(ppm_t_options_sop)),       INTENT(INOUT)   :: opts
          !!! options
          INTEGER,                              INTENT(  OUT)   :: info

          !optional arguments

          CLASS(ppm_t_main_abstr), TARGET, OPTIONAL               :: Adapt_Field
          !!! Field, or property, on which which the monitor function is applied.
          !!! This is thus the field that the particles will try to best discretize.
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
          TYPE(DTYPE(sop_t_stats)),POINTER,OPTIONAL,INTENT(OUT) :: stats
          !!! statistics on output

          ! argument-functions need an interface
          INTERFACE
              !Monitor function
              FUNCTION D_fun(f1,dfdx,opts,f2)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT DTYPE(ppm_t_options_sop)
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK)                                   :: D_fun
                  REAL(MK),                   INTENT(IN)     :: f1
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)     :: dfdx
                  TYPE(DTYPE(ppm_t_options_sop)),INTENT(IN)  :: opts
                  REAL(MK),OPTIONAL,          INTENT(IN)     :: f2
              END FUNCTION D_fun

              !Function that returns the width of the narrow band
              FUNCTION nb_fun(kappa,scale_D)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK)                                   :: nb_fun
                  REAL(MK),                INTENT(IN)        :: kappa
                  REAL(MK),                INTENT(IN)        :: scale_D
              END FUNCTION nb_fun

              !Field function (usually known only during initialisation)
              FUNCTION wp_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)     :: pos
                  REAL(MK)                                   :: wp_fun
              END FUNCTION wp_fun

              !Level function (usually known only during initialisation)
              FUNCTION level_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)     :: pos
                  REAL(MK)                                   :: level_fun
              END FUNCTION level_fun

              !Gradient of the field func. (usually known only during initialisation)
              FUNCTION wp_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK),DIMENSION(ppm_dim)                :: wp_grad_fun
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)     :: pos
              END FUNCTION wp_grad_fun

              !Gradient of the level func. (usually known only during initialisation)
              FUNCTION level_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  IMPORT                                     :: MK
                  IMPLICIT NONE
                  REAL(MK),DIMENSION(ppm_dim)                :: level_grad_fun
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)     :: pos
              END FUNCTION level_grad_fun
          END INTERFACE

          ! local variables
          INTEGER                                    :: i,ip,ineigh,iq
          INTEGER                                    :: iunit
          CHARACTER(LEN=256)                         :: filename
          REAL(MK)                                   :: dist2
          REAL(MK)                                   :: Psi_threshold
          INTEGER                                    :: num_it

          REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: D_old => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: Dtilde_old => NULL()
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

          INTEGER,      DIMENSION(:),   POINTER      :: nvlist => NULL()
          INTEGER,      DIMENSION(:,:), POINTER      :: vlist  => NULL()

          INTEGER,      DIMENSION(:),   POINTER      :: nvlist_cross => NULL()
          INTEGER,      DIMENSION(:,:), POINTER      :: vlist_cross => NULL()
          INTEGER                                    :: nneighmax_cross
          INTEGER                                    :: nneighmin_cross
          LOGICAL                                    :: need_derivatives
          INTEGER                                    :: memory_used

          CLASS(ppm_t_operator_discr_),POINTER       :: op => NULL()

          TYPE(DTYPE(ppm_t_sop))                     :: Particles_old
          TYPE(ppm_t_options_op)                     :: opts_op
          TYPE(ppm_t_operator)                       :: Interp
          CLASS(ppm_t_operator_discr),POINTER        :: InterpDiscr => NULL()
          CLASS(DTYPE(ppm_t_part_prop)_),POINTER     :: nn_sq => NULL()
          CLASS(ppm_t_discr_data),POINTER            :: discr_data => NULL()

          INTEGER, DIMENSION(ppm_dim)                :: degree
          REAL(ppm_kind_double),DIMENSION(1)         :: coeffs

          REAL(MK), PARAMETER :: big=HUGE(1.0_MK)

          start_subroutine("sop_adapt")

          !-------------------------------------------------------------------------!
          ! Checks consistency of parameters
          !-------------------------------------------------------------------------!
          check_false(<#(PRESENT(wp_grad_fun) .AND. .NOT.PRESENT(wp_fun))#>,&
                  &  "provided analytical gradients but not analytical function values. This case is not yet implemented")
          IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
              IF (opts%D_needs_gradients) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,caller,   &
                      & 'provided analytical function values but no &
                      & function gradients. This case is not implemented',&
                      &  __LINE__,info)
                  GOTO 9999
              ENDIF
              IF (this%level_set) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,caller,   &
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
              stdout("Calling adapt_particles with conflicting options",&
                     "Warning: on exit, wp_lap will have nonsense values")
          ENDIF

          check_true(<#(PRESENT(wp_fun).OR.PRESENT(Adapt_Field))#>,&
              "Need to either provide an analytical function wp_fun, or a discretized Field, Adapt_Field, on which particles should adapt")

          !check data structure exists
          check_true(<#ASSOCIATED(this%xp)#>,&
                  "Particles structure had not been defined. Call allocate first")

          !check that we are dealing with adaptive particles
          check_true(<#this%adaptive#>,&
                  &  "These particles have not been declared as adaptive")

          !check all particles are inside the computational domain
          check_true(<#this%flags(ppm_part_areinside)#>,&
                  "Some particles may be outside the domain. Apply BC first")

          !check that neighbour lists have been computed
          check_true(<#this%has_neighlist()#>,&
                  "Neighbour lists are not up-to-date. Call neighlists first")

          IF (opts%level_set) THEN
              IF (.NOT. this%level_set) THEN
                  fail("Need to enable level-set for Particles first")
              ENDIF
              IF (PRESENT(level_fun) .NEQV. PRESENT(wp_fun)) THEN
                  fail("Incompatible optional arguments")
              ENDIF
          ENDIF

          !-------------------------------------------------------------------------!
          ! Determines whether we will need to approximate derivatives
          !-------------------------------------------------------------------------!
          need_derivatives=.FALSE.
          IF (.NOT.PRESENT(wp_grad_fun)) THEN
              IF (this%level_set .OR. opts%D_needs_gradients) &
                  need_derivatives=.TRUE.
          ENDIF

          IF (PRESENT(stats)) stats%min_sv=big

          ! Check that the scalar field on which particles are supposed to adapt
          ! has been defined or is provided by an analytical function
          IF (.NOT.PRESENT(wp_fun)) THEN
              !Points this%adapt_wp to the array that stores the discretized
              !resolution function.
              SELECT TYPE(Adapt_Field)
              CLASS IS (ppm_t_field_)
              CALL this%get_discr(Adapt_Field,discr_data,info)
                  SELECT TYPE(discr_data)
                  CLASS IS (DTYPE(ppm_t_part_prop)_)
                      this%adapt_wp => discr_data
                  END SELECT
              CLASS IS (DTYPE(ppm_t_part_prop)_)
                  this%adapt_wp => Adapt_Field
              CLASS DEFAULT
                  fail("Wrong type for argument Adapt_Field")
              END SELECT
              or_fail("Could not access data array for Adapt_Field")
              check_associated(<#this%adapt_wp#>)
              check_true(<#this%adapt_wp%flags(ppm_ppt_ghosts)#>,&
                      &  "need to get the ghosts for adapt_wp")
          ENDIF

          IF (opts%level_set) THEN
              ! Check that the level set has been defined or
              ! is provided by an analytical function
              IF (PRESENT(level_fun)) THEN
                  check_true(<#PRESENT(level_grad_fun)#>,&
                          &  "need analytical gradients for level function")
              ELSE
                  !check if a scalar property has already been specified as the argument
                  !for the resolution function (i.e. the particles will adapt to resolve
                  !this property well)
                  check_associated(<#this%level#>,&
                          &  "need to define level set property first")
                  check_associated(<#this%level_grad#>,&
                          &  "need to define level set gradient property first")
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
                  IF (.NOT. ASSOCIATED(this%adapt_wp_grad)) THEN
                      !no array has already been specified for wp_grad
                      !need to allocate one
                      CALL this%create_prop(info,part_prop=this%adapt_wp_grad,&
                      &    lda=ppm_dim,name='adapt_wp_grad')
                      or_fail("could not create property for adapt_wp_grad")
                  ELSE
                      IF (.NOT.this%adapt_wp_grad%flags(ppm_ppt_partial)) THEN
                          CALL this%realloc_prop(this%adapt_wp_grad,&
                          &    info,lda=ppm_dim)
                          or_fail("failed to reallocate property adapt_wp_grad")
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

          !TODO
          CALL this%compute_D(D_fun,opts,info,     &
              wp_fun=wp_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
              level_grad_fun=level_grad_fun,nb_fun=nb_fun,stats=stats)
          or_fail("sop_compute_D failed")

          !REMOVME???
          CALL this%get(this%rcp,rcp,info)
          or_fail("could not get pointer to rcp")

          CALL this%get(this%D,D,info,read_only=.TRUE.)
          or_fail("could not get pointer to D")

          DO ip=1,this%Npart
              rcp(ip) = opts%rcp_over_D * D(ip)
          ENDDO
                  IF (opts%level_set) THEN
                      IF (PRESENT(wp_fun)) THEN
                          CALL this%get_xp(xp,info)
                          or_fail("get_xp")
                          DO ip=1,this%Npart
                              IF (ABS(level_fun(xp(1:ppm_dim,ip))) .GT. &
                                  &   opts%nb_width2*nb_fun(wp_fun(xp(1:ppm_dim,ip)),&
                                  &                                       opts%scale_D)) THEN
                                  rcp(ip) = 1.3_MK * rcp(ip)
                              ENDIF
                          ENDDO
                          CALL this%set_xp(xp,info,read_only=.TRUE.)
                          or_fail("set_xp")
                      ELSE
                          CALL this%get(this%level,level,info,read_only=.TRUE.)
                          or_fail("get level")
                          CALL this%get(this%adapt_wp,wp,info,read_only=.TRUE.)
                          or_fail("get adapt_wp")
                          DO ip=1,this%Npart
                              IF (ABS(level(ip)) .GT. &
                                  &   opts%nb_width2*nb_fun(wp(ip),opts%scale_D)) THEN
                                  rcp(ip) = 1.3_MK * rcp(ip)
                              ENDIF
                          ENDDO
                      ENDIF
                  ENDIF

#if debug_verbosity > 1
          WRITE(filename,'(A,I0)') 'P_SOP_Init_',this%itime
          CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif

          CALL this%updated_cutoff(1.3_MK*this%ghostlayer,info)
          or_fail("particles_updated_cutoff failed")


          !need the ghosts for D_old to be correct
          CALL this%map_ghosts(info)
          or_fail("map_ghosts failed")

          CALL this%comp_neighlist(info)
          or_fail("Failed to compute Verlet lists")
          !REMOVME???

          !!-------------------------------------------------------------------------!
          !! Find the nearest neighbour for each real particle
          !! (used in the definition of the primitive functions, for
          !! interpolation. Not needed if the functions are known analytically)
          !! We do it now, since we already have neighbour lists for that particle set
          !!-------------------------------------------------------------------------!
          IF (.NOT.PRESENT(wp_fun)) THEN
              degree =  0
              coeffs = 1.0_MK

              CALL Interp%create(ppm_dim,coeffs,degree,info,name="interpolant")
              or_fail("could not create interpolating operator")

              CALL opts_op%create(ppm_param_op_dcpse,info,&
              &    interp=.TRUE.,order=opts%order_approx,c=REAL(opts%c,ppm_kind_double))
              or_fail("failed to initialize option object for operator")

              !compute nearest-neighbour distances
              CALL this%create_prop(info,part_prop=nn_sq,name='nn_sq')
              or_fail("could not create property for nn_sq")
          ENDIF

          !-------------------------------------------------------------------------!
          ! Copy particle positions, field values and D
          !-------------------------------------------------------------------------!

          !Move data from the original particle set to a temporary particle set.
          !(this is done by doing a shallow copy and nullifying the pointers one
          !by one...)

          !TODO
          !Particles_old = this

          !These are kept on Particles_old
          this%xp => NULL()
          this%props => NULL()
          this%D   => NULL()
          this%Dtilde   => NULL()
          this%rcp => NULL()
          this%adapt_wp => NULL()
          this%level => NULL()
          this%level_grad => NULL()
          this%gi => NULL()
          this%field_ptr => NULL()

          !These are kept on "this"
          Particles_old%ops => NULL()
          Particles_old%maps => NULL()
          Particles_old%pcost => NULL()
          Particles_old%stats => NULL()
          Particles_old%neighs => NULL()
          Particles_old%flags(ppm_part_neighlists) = .FALSE.

          !move DC operators from old to new particles
          ! (only move their definitions
          !IF (ASSOCIATED(Particles_old%ops)) THEN
              !this%ops => Particles_old%ops
              !Particles_old%ops => NULL()
              !op => this%ops%begin()
              !DO WHILE (ASSOCIATED(op))
                  !op%flags(ppm_ops_iscomputed) = .FALSE.
                  !op => this%ops%next()
              !ENDDO
          !ENDIF

          ! All the particle properties are kept on the old set Particles_old
          ! and will be interpolated back onto the new set after self-organization.
          ! In the meantime, the new set of particles has no property to carry around.

          !link the old particles to the new ones in the data structure:
          !this%Particles_cross => Particles_old

          !Allocate a new array for particle positions
          ldc(1) = ppm_dim
          ldc(2) = this%Mpart
          CALL ppm_alloc(this%xp,ldc,ppm_param_alloc_fit,info)
          or_fail_alloc("xp")

          !Define some properties on the new set
          CALL this%create_prop(info,part_prop=this%D,&
          &    name=Particles_old%D%name,with_ghosts=.TRUE.)
          or_fail("could not create property for D")

          CALL this%create_prop(info,part_prop=this%Dtilde,&
          &    name=Particles_old%Dtilde%name,with_ghosts=.TRUE.)
          or_fail("could not create property for Dtilde")

          CALL this%create_prop(info,part_prop=this%rcp,&
          &    name=Particles_old%rcp%name,with_ghosts=.TRUE.)
          or_fail("could not create property for rcp")


          ! We ask for read_only arrays because we know that we dont need
          ! a ghost mapping after modifying them.
          CALL this%get_xp(xp,info,with_ghosts=.TRUE.)
          CALL this%get(this%D,D,info,with_ghosts=.TRUE.,read_only=.TRUE.)
          CALL this%get(this%Dtilde,Dtilde,info,with_ghosts=.TRUE.,read_only=.TRUE.)
          CALL this%get(this%rcp,rcp,info,with_ghosts=.TRUE.,read_only=.TRUE.)
          CALL Particles_old%get_xp(xp_old,info,with_ghosts=.TRUE.)
          CALL Particles_old%get(Particles_old%D,D_old,info,&
              with_ghosts=.TRUE.,read_only=.TRUE.)
          CALL Particles_old%get(Particles_old%Dtilde,Dtilde_old,info,&
              with_ghosts=.TRUE.,read_only=.TRUE.)
          CALL Particles_old%get(Particles_old%rcp,rcp_old,info,&
              with_ghosts=.TRUE.,read_only=.TRUE.)

          DO ip=1,this%Mpart
              xp(1:ppm_dim,ip) = xp_old(1:ppm_dim,ip)
              D(ip) = D_old(ip)
              Dtilde(ip) = Dtilde_old(ip)
              rcp(ip) = rcp_old(ip)
          ENDDO

          CALL this%set_xp(xp,info,read_only=.TRUE.)
          CALL Particles_old%set_xp(xp_old,info,read_only=.TRUE.)

          !If we use a level-set method with narrow band
          ! we have to keep track of the distance function
          ! DURING adaptation. We thus have to copy the
          ! level function (and its gradient)
          IF (Particles_old%level_set .AND. .NOT.PRESENT(level_fun)) THEN
              CALL this%create_prop(info,part_prop=this%level,&
                  name=Particles_old%level%name,zero=.TRUE.)
              or_fail("could not create property for level")
              CALL this%create_prop(info,part_prop=this%level_grad,&
                  name=Particles_old%level_grad%name,lda=ppm_dim,zero=.TRUE.)
              or_fail("could not create property for level_grad")

              CALL this%get(this%level,level,info,read_only=.FALSE.)
              CALL this%get(this%level_grad,level_grad,info,read_only=.FALSE.)
              CALL Particles_old%get(Particles_old%level,level_old,info,read_only=.TRUE.)
              CALL Particles_old%get(Particles_old%level_grad,level_grad_old,&
                  info,read_only=.TRUE.)
              DO ip=1,this%Npart
                  level(ip) = level_old(ip)
                  level_grad(1:ppm_dim,ip) = level_grad_old(1:ppm_dim,ip)
              ENDDO
          ENDIF

          IF (opts%level_set .AND. .NOT.PRESENT(wp_fun)) THEN
              CALL this%create_prop(info,part_prop=this%adapt_wp,&
                  name=this%adapt_wp%name)
              or_fail("could not create property for adapt_wp")


              CALL this%get(this%adapt_wp,wp,info,read_only=.FALSE.)
              CALL Particles_old%get(Particles_old%adapt_wp,wp_old,info,read_only=.TRUE.)
              DO ip=1,this%Npart
                  wp(ip)= wp_old(ip)
              ENDDO
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

          !TODO
          !IF (this%level_set) THEN
              !CALL sop_gradient_descent_ls(Particles_old,Particles, &
                  !nvlist_cross,vlist_cross,                &
                  !nneighmin_cross,nneighmax_cross,num_it,opts,info,wp_fun=wp_fun,&
                  !D_fun=D_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
                  !level_grad_fun=level_grad_fun,threshold=Psi_threshold,&
                  !need_deriv=need_derivatives,nb_fun=nb_fun,stats=stats)
          !ELSE
              !CALL sop_gradient_descent(Particles_old,Particles, &
                  !nvlist_cross,vlist_cross,                &
                  !nneighmin_cross,nneighmax_cross,num_it,opts,info,wp_fun=wp_fun,&
                  !D_fun=D_fun,wp_grad_fun=wp_grad_fun,threshold=Psi_threshold,&
                  !need_deriv=need_derivatives,stats=stats)
          !ENDIF

          !TODO
          !num_it is uninitialized
          !yaser
          IF (PRESENT(stats)) THEN
             stats%nb_grad_desc_steps = stats%nb_grad_desc_steps !+ num_it
          ENDIF

          !!-------------------------------------------------------------------------!
          !! Compute field values at new particle locations
          !!-------------------------------------------------------------------------!

          IF (PRESENT(wp_fun)) THEN
              !!---------------------------------------------------------------------!
              !! If the field is known as a function, do nothing
              !!---------------------------------------------------------------------!
          ELSE
              !!---------------------------------------------------------------------!
              !! Otherwise, use particle-to-particle interpolation
              !!---------------------------------------------------------------------!

              !TODO
      !        CALL sop_interpolate(Particles_old,Particles,opts,info)

              CALL Interp%discretize_on(this,InterpDiscr,opts_op,info,&
                  Discr_src=Particles_old,nn_sq=nn_sq)
              or_fail("could not discretize interpolating operator")

              CALL opts_op%destroy(info)
              or_fail("failed to destroy opts_op")

          ENDIF

#if debug_verbosity > 1
          WRITE(filename,'(A,I0)') 'P_SOP_Done_',this%itime
          CALL ppm_vtk_particle_cloud(filename,this,info)
#endif

          !IF (ppm_rank .EQ. 0) THEN
              !WRITE(filename,'(A,A)') TRIM(debugdir),'Stats.dat'
              !iunit=20
              !OPEN(UNIT=iunit,FILE=filename,FORM='FORMATTED',&
                  !ACCESS='sequential',POSITION='APPEND',IOSTAT=info)
              !WRITE(iunit,'(I8,2X,I8,2X)') this%Npart,nb_grad_desc_steps
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
#if debug_verbosity > 2
          !!evaluate memory usage

          !stdout("Memory used for Particles_old:")

          !memory_used = 8*SIZE(Particles_old%xp)
          !WRITE(*,*) '          xp: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 2*(SIZE(Particles_old%vlist)+SIZE(Particles_old%nvlist))
          !WRITE(*,*) '       vlist: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 2*(SIZE(Particles_old%vlist_cross)+SIZE(Particles_old%nvlist_cross))
          !WRITE(*,*) '      xvlist: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !DO i=1,Particles_old%max_wpiid
              !IF (ASSOCIATED(Particles_old%wpi(i)%vec)) THEN
                  !memory_used = memory_used + 2*SIZE(Particles_old%wpi(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wpi: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 0
          !DO i=1,Particles_old%max_wpsid
              !IF (ASSOCIATED(Particles_old%wps(i)%vec)) THEN
                  !memory_used = memory_used + 8*SIZE(Particles_old%wps(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wps: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 0
          !DO i=1,Particles_old%max_wpvid
              !IF (ASSOCIATED(Particles_old%wpv(i)%vec)) THEN
                  !memory_used = memory_used + 8*SIZE(Particles_old%wpv(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wpv: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !WRITE(*,*) 'Memory used for Particles:'

          !memory_used = 8*SIZE(this%xp)
          !WRITE(*,*) '          xp: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 2*(SIZE(this%vlist)+SIZE(this%nvlist))
          !WRITE(*,*) '       vlist: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 2*(SIZE(this%vlist_cross)+SIZE(this%nvlist_cross))
          !WRITE(*,*) '      xvlist: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !DO i=1,this%max_wpiid
              !IF (ASSOCIATED(this%wpi(i)%vec)) THEN
                  !memory_used = memory_used + 2*SIZE(this%wpi(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wpi: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 0
          !DO i=1,this%max_wpsid
              !IF (ASSOCIATED(this%wps(i)%vec)) THEN
                  !memory_used = memory_used + 8*SIZE(this%wps(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wps: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !memory_used = 0
          !DO i=1,this%max_wpvid
              !IF (ASSOCIATED(this%wpv(i)%vec)) THEN
                  !memory_used = memory_used + 8*SIZE(this%wpv(i)%vec)
              !ENDIF
          !ENDDO
          !WRITE(*,*) '         wpv: ', memory_used
          !memory_used_total = memory_used_total+memory_used

          !WRITE(*,*) ' Grand TOTAL: ', memory_used_total,' i.e.  ',&
              !REAL(memory_used_total,MK)/1E6,' MB'
#endif


          CALL Particles_old%destroy(info)
          or_fail("failed to destroy old set of particles")

          CALL InterpDiscr%destroy(info)
          or_fail("failed to destroy discretized interpolation operator")

          CALL Interp%destroy(info)
          or_fail("failed to destroy interpolation operator")

          IF (ASSOCIATED(nvlist_cross)) DEALLOCATE(nvlist_cross)
          IF (ASSOCIATED(vlist_cross)) DEALLOCATE(vlist_cross)

          IF (ASSOCIATED(xp))    CALL ppm_write(ppm_rank,caller,'forgot to Set xp',info)
          IF (ASSOCIATED(xp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set xp_old',info)
          IF (ASSOCIATED(D))     CALL ppm_write(ppm_rank,caller,'forgot to Set D',info)
          IF (ASSOCIATED(D_old)) CALL ppm_write(ppm_rank,caller,'forgot to Set D_old',info)
          IF (ASSOCIATED(xp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set xp_old',info)
          IF (ASSOCIATED(wp_old))CALL ppm_write(ppm_rank,caller,'forgot to Set wp_old',info)
          IF (ASSOCIATED(level)) CALL ppm_write(ppm_rank,caller,'forgot to Set level',info)
          IF (ASSOCIATED(level_old))CALL ppm_write(ppm_rank,caller,'forgot to Set level_old',info)
          IF (ASSOCIATED(level_grad))CALL ppm_write(ppm_rank,caller,'forgot to Set level_grad',info)
          IF (ASSOCIATED(level_grad_old)) &
            CALL ppm_write(ppm_rank,caller,'forgot to Set level_grad_old',info)


          end_subroutine()
      END SUBROUTINE DTYPE(sop_self_organize)
