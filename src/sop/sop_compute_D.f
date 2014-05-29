      !!!----------------------------------------------------------------------------!
      !!! Computes resolution field D
      !!!----------------------------------------------------------------------------!

      SUBROUTINE DTYPE(sop_compute_D)(this,D_fun,opts,info, &
      &          wp_fun,wp_grad_fun,level_fun,level_grad_fun,nb_fun, &
      &          only_D_tilde,stats)

          IMPLICIT NONE
#ifdef __MPI
          INCLUDE 'mpif.h'
#endif

          DEFINE_MK()
          ! arguments
          CLASS(DTYPE(ppm_t_sop)),              INTENT(INOUT)   :: this
          !!! particles
          TYPE(DTYPE(ppm_t_options_sop)),       INTENT(IN   )   :: opts
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
          TYPE(DTYPE(sop_t_stats)),POINTER,OPTIONAL,INTENT(OUT) :: stats
          !!! statistics on output

          ! argument-functions need an interface
          INTERFACE
              !Monitor function
              FUNCTION D_fun(f1,dfdx,opts,f2)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  import DTYPE(ppm_t_options_sop)
                  DEFINE_MK()
                  REAL(MK)                               :: D_fun
                  REAL(MK),                   INTENT(IN) :: f1
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN) :: dfdx
                  TYPE(DTYPE(ppm_t_options_sop)),   INTENT(IN) :: opts
                  REAL(MK),OPTIONAL,          INTENT(IN) :: f2
              END FUNCTION D_fun

              !Function that returns the width of the narrow band
              FUNCTION nb_fun(kappa,scale_D)
                  USE ppm_module_interfaces
                  DEFINE_MK()
                  REAL(MK)                             :: nb_fun
                  REAL(MK),                INTENT(IN)  :: kappa
                  REAL(MK),                INTENT(IN)  :: scale_D
              END FUNCTION nb_fun

              !Field function (usually known only during initialisation)
              FUNCTION wp_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)        :: pos
                  REAL(MK)                                      :: wp_fun
              END FUNCTION wp_fun

              !Level function (usually known only during initialisation)
              FUNCTION level_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)        :: pos
                  REAL(MK)                                      :: level_fun
              END FUNCTION level_fun

              !Gradient of the field func. (usually known only during initialisation)
              FUNCTION wp_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim)                      :: wp_grad_fun
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
              END FUNCTION wp_grad_fun

              !Gradient of the level func. (usually known only during initialisation)
              FUNCTION level_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_interfaces
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim)                      :: level_grad_fun
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: pos
              END FUNCTION level_grad_fun
          END INTERFACE

          ! local variables
          INTEGER                                    :: i,ip,ineigh,iq
          CHARACTER(LEN=64)                          :: myformat

!           REAL(MK),     DIMENSION(:,:), POINTER      :: xp_old => NULL()
!           REAL(MK),     DIMENSION(:),   POINTER      :: wp_old => NULL()
!           REAL(MK),     DIMENSION(:),   POINTER      :: D_old => NULL()
!           REAL(MK),     DIMENSION(:),   POINTER      :: rcp_old => NULL()
!           REAL(MK),     DIMENSION(:),   POINTER      :: level_old => NULL()
!           REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad_old => NULL()

          REAL(MK),     DIMENSION(:,:), POINTER      :: xp => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: rcp => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: D => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: Dtilde => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: wp => NULL()
          REAL(MK),     DIMENSION(:,:), POINTER      :: wp_grad => NULL()
          REAL(MK),     DIMENSION(:),   POINTER      :: level => NULL()
          REAL(MK),     DIMENSION(:,:), POINTER      :: level_grad => NULL()

!           REAL(MK),     DIMENSION(:,:), POINTER      :: eta => NULL()
          REAL(MK)                                   :: min_D
          LOGICAL                                    :: need_derivatives
          REAL(MK),     DIMENSION(ppm_dim)           :: dummy_grad
          INTEGER                                    :: topo_id !,eta_id
          REAL(ppm_kind_double),     DIMENSION(ppm_dim)           :: coeffs
          INTEGER,      DIMENSION(ppm_dim)           :: order
          INTEGER,      DIMENSION(ppm_dim*ppm_dim)   :: degree
          REAL(MK),     DIMENSION(ppm_dim)           :: wp_grad_fun0
          REAL(MK)                                   :: alpha
          CLASS(DTYPE(ppm_t_neighlist)_),POINTER:: NList
          TYPE(ppm_t_operator)                  :: Op
          TYPE(ppm_t_options_op)                :: opts_op
          CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
          INTEGER, DIMENSION(:),    POINTER     :: nvlist => NULL()
          INTEGER, DIMENSION(:,:),  POINTER     :: vlist => NULL()

          start_subroutine("sop_compute_D")

          dummy_grad=0._MK
          topo_id = this%active_topoid

          !-------------------------------------------------------------------------!
          ! Checks consistency of parameters
          !-------------------------------------------------------------------------!
          check_false(<#(PRESENT(wp_grad_fun) .AND. .NOT.PRESENT(wp_fun))#>,&
          &  "provided analytical gradients but not analytical function values. This case is not yet implemented")
          IF (.NOT.PRESENT(wp_grad_fun) .AND. PRESENT(wp_fun)) THEN
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
          !-------------------------------------------------------------------------!
          ! Perform consistency checks
          !-------------------------------------------------------------------------!
          !check data structure exists
          IF (.NOT.ASSOCIATED(this%xp)) THEN
              fail("Particles structure had not been defined. Call allocate first")
          ENDIF

          !check that we are dealing with adaptive particles
          IF (.NOT.this%adaptive) THEN
              fail("These particles have not been declared as adaptive")
          ENDIF

          !check all particles are inside the computational domain
          check_true(<#this%flags(ppm_part_areinside)#>,&
                  "Some particles may be outside the domain. Apply BC first")

          IF (opts%level_set) THEN
              IF (.NOT. this%level_set) THEN
                  fail("Need to enable level-set for Particles first")
              ENDIF
              IF (PRESENT(wp_fun) .AND..NOT.PRESENT(level_fun)) THEN
                  fail("Need to provide analytical level function")
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

          ! Check that the scalar field on which particles are supposed to adapt
          ! has been defined or is provided by an analytical function
          IF (.NOT.PRESENT(wp_fun)) THEN
              !check if a scalar property has already been specified as the argument
              !for the resolution function (i.e. the particles will adapt to resolve
              !this property well)
              check_associated(<#this%adapt_wp#>)
          ENDIF

          IF (opts%level_set) THEN
              ! Check that the level set has been defined or
              ! is provided by an analytical function
              IF (PRESENT(level_fun)) THEN
                  IF (.NOT.PRESENT(level_grad_fun)) THEN
                      fail("need to analytical gradients for level function")
                  ENDIF
              ELSE
                  !check if a scalar property has already been specified as the argument
                  !for the resolution function (i.e. the particles will adapt to resolve
                  !this property well)
                  IF (.NOT.ASSOCIATED(this%level)) THEN
                      fail("need to define level_id first")
                  ENDIF
                  IF (.NOT.ASSOCIATED(this%level_grad)) THEN
                      fail("need to define level_grad first")
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
                  check_associated(<#this%adapt_wp_grad#>)
                  IF (.NOT.this%adapt_wp_grad%flags(ppm_ppt_partial)) THEN
                     CALL this%realloc_prop(this%adapt_wp_grad,info,lda=ppm_dim)
                     or_fail("failed to reallocate property adapt_wp_grad")
                  ENDIF
              ENDIF
          ENDIF

          NList => this%get_neighlist()
          check_associated(NList,"Could not find neighbour list for that particle set")
          ! crash if not enough neighbours
          IF (need_derivatives .AND. NList%nneighmin.LT.opts%nneigh_critical) THEN
              WRITE(cbuf,'(2(A,I5,2X))') 'nneigh_critical ',opts%nneigh_critical,&
                  'nneigh_toobig ',opts%nneigh_toobig
              CALL ppm_write(ppm_rank,caller,cbuf,info)
              WRITE(cbuf,'(A,I5)') 'We have nneighmin = ',NList%nneighmin
              CALL ppm_write(ppm_rank,caller,cbuf,info)
              fail("Not enough neighbours")
          ENDIF

          !!-------------------------------------------------------------------------!
          !! Compute D (desired resolution)
          !!-------------------------------------------------------------------------!
          ! (re)allocate Dtilde
          IF (.NOT. ASSOCIATED(this%Dtilde)) THEN
              CALL this%create_prop(info,part_prop=this%Dtilde,lda=1,&
                  name='Dtilde',with_ghosts=.TRUE.)
              or_fail("could not create property for Dtilde")
          ELSE
              IF (.NOT.this%Dtilde%flags(ppm_ppt_partial)) THEN
                  CALL this%realloc_prop(this%Dtilde,info,lda=1,with_ghosts=.TRUE.)
                  or_fail("failed to reallocate property Dtilde")
              ENDIF
          ENDIF

          !!-------------------------------------------------------------------------!
          !! Case where we need to approximate derivatives
          !!-------------------------------------------------------------------------!
          if_needs_derivatives: IF (need_derivatives) THEN
              ! Compute DC operator for gradients
              coeffs=1.D0; degree = 0
              FORALL(i=1:ppm_dim) degree((i-1)*ppm_dim+i)=1 !Gradient

              CALL Op%create(ppm_dim,coeffs,degree,info,name="Gradient")
              or_fail("Failed to create Gradient operator")
              CALL opts_op%create(ppm_param_op_dcpse,info,order=3,&
                  c=REAL(opts%c,ppm_kind_double),vector=.true.)
              or_fail("failed to initialize option object for operator")
              CALL Op%discretize_on(this,DCop,opts_op,info)
              or_fail("Failed to discretize Gradient operator on particle set")

          ENDIF if_needs_derivatives

          !-------------------------------------------------------------------------!
          ! Compute D_tilde
          !-------------------------------------------------------------------------!

          if_D_needs_grad: IF (opts%D_needs_gradients) THEN
              check_false(<#opts%level_set#>,&
                      "D_needs_gradients with level_sets: not supported")

              IF (.NOT.PRESENT(wp_grad_fun)) THEN
                  !Compute gradients
                  CALL DCop%compute(this%adapt_wp,this%adapt_wp_grad,info)
                  or_fail("Failed to compute gradient of adapt field")

                  CALL DCop%destroy(info)
                  or_fail("Failed to destroy DC operator")
              ENDIF

              !Compute Dtilde on real particles

              CALL this%get_xp(xp,info)
              or_fail("Failed to access xp")

              CALL this%get(this%Dtilde,Dtilde,info)
              or_fail("Failed to access Dtilde")

              IF (PRESENT(wp_fun)) THEN
                  DO ip=1,this%Npart
                      wp_grad_fun0= wp_grad_fun(xp(1:ppm_dim,ip))
                      Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),wp_grad_fun0,opts)
                  ENDDO
              ELSE
                  CALL this%get(this%adapt_wp,wp,info,read_only=.TRUE.)
                  or_fail("Failed to access wp")

                  CALL this%get(this%adapt_wp_grad,wp_grad,info,read_only=.TRUE.)
                  or_fail("Failed to access wp_grad")

                  DO ip=1,this%Npart
                      Dtilde(ip) = D_fun(wp(ip),wp_grad(1:ppm_dim,ip),opts)
                  ENDDO
              ENDIF
              CALL this%set_xp(xp,info,read_only=.TRUE.)

              ! Get ghosts for D_tilde
              CALL this%map_ghosts(info)
              or_fail("ghost mappings failed")

          ELSE ! .NOT. D_needs_grad

              CALL this%get(this%Dtilde,Dtilde,info,with_ghosts=.TRUE.,read_only=.TRUE.)
              or_fail("Failed to access Dtilde")
              IF (PRESENT(wp_fun)) THEN
                  CALL this%get_xp(xp,info,with_ghosts=.TRUE.)
                  or_fail("Failed to access xp")
                  IF (opts%level_set) THEN
                      DO ip=1,this%Mpart
                          Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),dummy_grad,&
                              opts,level_fun(xp(1:ppm_dim,ip)))
                      ENDDO
                  ELSE
                      DO ip=1,this%Mpart
                          Dtilde(ip) = D_fun(wp_fun(xp(1:ppm_dim,ip)),dummy_grad,opts)
                      ENDDO
                  ENDIF
                  CALL this%set_xp(xp,info,read_only=.TRUE.)
              ELSE
                  CALL this%get(this%adapt_wp,wp,info,with_ghosts=.TRUE.,read_only=.TRUE.)
                  or_fail("Failed to access wp")
                  IF (opts%level_set) THEN
                      CALL this%get(this%level,level,info,with_ghosts=.TRUE.,read_only=.TRUE.)
                      or_fail("Failed to access level")
                      DO ip=1,this%Mpart
                          Dtilde(ip) = D_fun(wp(ip),dummy_grad,opts,level(ip))
                      ENDDO
                  ELSE
                      DO ip=1,this%Mpart
                          Dtilde(ip) = D_fun(wp(ip),dummy_grad,opts)
                      ENDDO
                  ENDIF
              ENDIF

          ENDIF if_D_needs_grad

          !-------------------------------------------------------------------------!
          ! Rescale D_tilde
          !-------------------------------------------------------------------------!
          CALL this%get(this%Dtilde,Dtilde,info,with_ghosts=.TRUE.,read_only=.TRUE.)
          or_fail("failed to access Dtilde")
          DO ip = 1,this%Mpart
              Dtilde(ip) = MIN(MAX(Dtilde(ip),opts%minimum_D),opts%maximum_D)
          ENDDO

          IF (PRESENT(only_D_tilde)) THEN
              IF (only_D_tilde) GOTO 8000
          ENDIF

          !---------------------------------------------------------------------!
          ! Update cutoff radii
          !---------------------------------------------------------------------!
          CALL this%get(this%rcp,rcp,info)
          or_fail("failed to access rcp")

          CALL this%get(this%Dtilde,Dtilde,info,read_only=.TRUE.)
          or_fail("failed to access Dtilde")

          DO ip=1,this%Npart
              !rcp(ip) = opts%rcp_over_D * Dtilde(ip)
              !TESTING THIS:
              rcp(ip) = Dtilde(ip)
          ENDDO

          CALL this%updated_cutoff(MAXVAL(rcp),info)
          or_fail("Failed to notify PPM that the cutoff has changed")

          !---------------------------------------------------------------------!
          ! Update ghosts
          !---------------------------------------------------------------------!
          CALL this%map_ghosts(info)
          or_fail("ghost mappings failed")
          !---------------------------------------------------------------------!
          ! Update neighbour lists
          !---------------------------------------------------------------------!
          CALL this%comp_neighlist(info)
          or_fail("failed to compute neighbour lists")

          CALL this%get_vlist(nvlist,vlist,info)
          or_fail("failed to access neighbour lists")

          IF (.NOT. ASSOCIATED(this%D)) THEN
              CALL this%create_prop(info,part_prop=this%D,lda=1,name='D')
              or_fail("could not create property for D")
          ENDIF
          !---------------------------------------------------------------------!
          ! D^(n+1) = min(D_tilde^(n+1)(iq)) over all neighbours iq
          !---------------------------------------------------------------------!
          CALL this%get(this%D,D,info)
          or_fail("failed to access D")
          CALL this%get(this%Dtilde,Dtilde,info,read_only=.TRUE.,with_ghosts=.TRUE.)
          or_fail("failed to access Dtilde")
          CALL this%get_xp(xp,info,with_ghosts=.TRUE.)
          or_fail("failed to access xp")
          DO ip=1,this%Npart
              D(ip)=Dtilde(ip)
              DO ineigh=1,nvlist(ip)
                  iq=vlist(ineigh,ip)

              !either this ....
                  IF (Dtilde(iq).GE.Dtilde(ip)) CYCLE

                  !alpha = (sqrt(sum((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2))/&
                      !Dtilde(iq)-1._MK) / (opts%rcp_over_D-1._MK)
                  !if(alpha.le.0) then
                      !D(ip)=MIN(D(ip),Dtilde(iq))
                  !else
                      !D(ip)=MIN(D(ip),sqrt(alpha)*Dtilde(ip)+(1._MK-sqrt(alpha))*Dtilde(iq))
                  !endif

                  alpha = (SQRT(SUM((xp(1:ppm_dim,ip)-xp(1:ppm_dim,iq))**2))/&
                      Dtilde(iq)-1._MK)
                  IF (alpha.LE.0) then
                     D(ip)=MIN(D(ip),Dtilde(iq))
                  ELSE
                     D(ip)=MIN(D(ip),2._MK**(alpha)*Dtilde(iq))
                  ENDIF

              !.... or this:
                  !IF (Dtilde(iq).LT.D(ip)) &
                      !D(ip)=Dtilde(iq)
              ENDDO
          ENDDO
          CALL this%set_xp(xp,info,read_only=.TRUE.)

          8000 CONTINUE


          end_subroutine()
      END SUBROUTINE DTYPE(sop_compute_D)
