      SUBROUTINE DTYPE(sop_gradient_descent)(Particles_old,Particles, &
      &          nvlist_cross,vlist_cross,nneighmin_cross,            &
      &          nneighmax_cross,num_it,opts,info,wp_fun,D_fun,       &
      &          wp_grad_fun,threshold,need_deriv,stats)
          !!!---------------------------------------------------------------------!
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
          !!!---------------------------------------------------------------------!

          USE ppm_module_mpi
          USE ppm_module_util_time
          USE ppm_module_inl_xset_vlist
          USE ppm_module_io_vtk
#ifdef __USE_LBFGS
          USE ppm_module_lbfgs
#endif
          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          TYPE(DTYPE(ppm_t_particles)), POINTER, INTENT(IN   )   :: Particles_old
          TYPE(DTYPE(ppm_t_particles)), POINTER, INTENT(INOUT)   :: Particles
          INTEGER,      DIMENSION(:),  POINTER,  INTENT(INOUT)   :: nvlist_cross
          INTEGER,      DIMENSION(:,:),POINTER,  INTENT(INOUT)   :: vlist_cross
          INTEGER,                               INTENT(INOUT)   :: nneighmax_cross
          INTEGER,                               INTENT(INOUT)   :: nneighmin_cross
          INTEGER,                               INTENT(  OUT)   :: num_it
          TYPE(DTYPE(sop_t_opts)), POINTER,      INTENT(IN   )   :: opts
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
          TYPE(DTYPE(sop_t_stats)),POINTER,OPTIONAL,INTENT(OUT)  :: stats
          !!! statistics on output
          ! argument-functions need an interface
          INTERFACE
              FUNCTION D_fun(f,dfdx,opts,lap_f)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  USE ppm_module_sop_typedef
                  DEFINE_MK()
                  REAL(MK)                                     :: D_fun
                  REAL(MK),                         INTENT(IN) :: f
                  REAL(MK),DIMENSION(ppm_dim),      INTENT(IN) :: dfdx
                  TYPE(DTYPE(sop_t_opts)),POINTER,         INTENT(IN) :: opts
                  REAL(MK),OPTIONAL,                INTENT(IN) :: lap_f
              END FUNCTION D_fun

              FUNCTION wp_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
                  REAL(MK),DIMENSION(ppm_dim)             :: wp_grad_fun
              END FUNCTION wp_grad_fun
              FUNCTION wp_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)        :: pos
                  REAL(MK)                                      :: wp_fun
              END FUNCTION wp_fun

          END INTERFACE

          ! local variables
          INTEGER                             :: ip, it_adapt,iq,ineigh,di,i
          INTEGER                             :: iunit
          CHARACTER(LEN = 256)                :: filename,cbuf
          CHARACTER(LEN = 256)                :: caller='sop_gradient_descent'
          REAL(ppm_kind_double) :: t0

          REAL(MK)                            :: step,step_previous,alpha1,alpha2
          REAL(MK)                            :: step_min,step_max,step_stall
          INTEGER                             :: it_adapt_max
          REAL(MK)                            :: Psi_max,Psi_global,Psi_global_old
          REAL(MK)                            :: gradPsi_max, gradPsi_thresh
          REAL(MK)                            :: Psi_1,Psi_2,Psi_thresh
          REAL(MK),DIMENSION(ppm_dim)         :: dist,dist2,dummy_grad
          REAL(MK),DIMENSION(:,:),POINTER     :: Gradient_Psi => NULL()
          INTEGER                             :: nneigh_adapt

          REAL(MK),DIMENSION(:,:),POINTER     :: xp => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: D => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: rcp => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: Dtilde => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: Dtilde_old => NULL()
          REAL(MK),DIMENSION(:,:),POINTER     :: xp_old => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: D_old => NULL()

          REAL(MK)                            :: tmpvar1,tmpvar2
          REAL(MK)                            :: minDold,minDtildeold
          REAL(MK)                            :: weight,weight_sum
          REAL(MK)                            :: almostzero
          INTEGER                             :: tmpvari1,tmpvari2
          LOGICAL                             :: need_derivatives
          INTEGER                             :: topo_id
          LOGICAL                             :: adding_particles
          INTEGER                             :: nb_spawn,nb_fuse

          !should be removed once the argument lists for the inl routines
          !have been updated to inhomogeneous ghostlayers
          REAL(MK),DIMENSION(2*ppm_dim)       :: ghostlayer

          !testing this
          INTEGER,DIMENSION(:),POINTER        :: move_part => NULL()
          INTEGER,DIMENSION(:),POINTER        :: fuse_part => NULL()
          REAL(MK)                            :: rcp_over_D_save
#ifdef __USE_LBFGS
          LOGICAL                             :: lbfgs_continue
          REAL(MK),DIMENSION(:),ALLOCATABLE   :: Work
          REAL(MK),DIMENSION(:),ALLOCATABLE   :: DIAG
          REAL(MK),DIMENSION(:),ALLOCATABLE   :: xp_lbfgs,grad_lbfgs
#endif

          REAL(MK), PARAMETER :: big=HUGE(1._MK)

          REAL(KIND(1.D0))                    :: t1,t2
          REAL(KIND(1.D0))                    :: ls_t1,ls_t2
          REAL(KIND(1.D0))                    :: add_t1,add_t2
          REAL(KIND(1.D0))                    :: del_t1,del_t2
          REAL(KIND(1.D0))                    :: compD_t1,compD_t2


#ifdef __USE_DEL_METHOD2
          opts%attractive_radius0 = 0.0_MK
          opts%fuse_radius = 0.00_MK
#endif
          !!-------------------------------------------------------------------------!
          !! Initialize
          !!-------------------------------------------------------------------------!
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

          it_adapt = 0
          adding_particles=.TRUE.
          Psi_max = big
          gradPsi_max = big
          Psi_global = big
          Psi_global_old = big

          ! Adaptivity stopping criterion
          IF (PRESENT(threshold)) THEN
              Psi_thresh = threshold
          ELSE
              Psi_thresh = 2.5_MK
          ENDIF
          gradPsi_thresh = 3E2

          IF (MK .EQ. KIND(1.D0)) THEN
              !step_min = 1D-8; step_stall = 1D-14
              step_min = 1D-3; step_stall = 1D-14
          ELSE
              !step_min = 1E-8; step_stall = 1E-14
              step_min = 1E-3; step_stall = 1E-14
          ENDIF
          step_max = 1.4_MK ! 0.1_MK
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

          !fuse_id = 0
          CALL particles_allocate_wpi(Particles,fuse_id,info,&
              iopt=ppm_param_alloc_fit,name='fuse_part')
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,caller,&
                  'particles_allocate_wpi failed', __LINE__,info)
              GOTO 9999
          ENDIF
          fuse_part => get_wpi(Particles,fuse_id)
          DO ip=1,Particles%Npart
              fuse_part(ip) = 0
          ENDDO
          fuse_part => set_wpi(Particles,fuse_id)

          !nb_neigh_id = 0
          CALL particles_allocate_wpi(Particles,nb_neigh_id,info,&
              iopt=ppm_param_alloc_fit,name='nb_neigh',zero=.TRUE.)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_alloc,caller,&
                  'particles_allocate_wpi failed', __LINE__,info)
              GOTO 9999
          ENDIF


              !FIXME: in theory, we should do apply_bc+remap+get_ghosts here
              ! (Fusing particles requires knowing the ghosts and particles
              ! have moved after the linesearch in the previous iteration)
              ! There MAY be a better way of doing this...
              !!---------------------------------------------------------------------!
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

          !!-------------------------------------------------------------------------!
          !! Gradient descent loop until stopping criterion is met
          !!-------------------------------------------------------------------------!
          adaptation_ok = .FALSE.
#ifdef __USE_LBFGS
          lbfgs_continue = .FALSE.
#endif
          !it_adapt_loop: DO WHILE (step .GT. step_stall .AND. &
                  !&          it_adapt .LT. it_adapt_max .AND. &
                  !&     ((Psi_max .GT. Psi_thresh)) .OR. adding_particles .OR.&
                  !&       gradPsi_max .GT. gradPsi_thresh)

          it_adapt_loop: DO WHILE (.NOT.adaptation_ok .OR. Psi_max .GT. Psi_thresh)

              it_adapt = it_adapt + 1

              !!---------------------------------------------------------------------!
              !! Update number of particles by fusion/insertion
              !! of particles (need vlist to be up-to-date)
              !! (on output, Mpart becomes meaningless.)
              !!---------------------------------------------------------------------!
              !NOTE: if we know D_tilde analytically, we should not supply
              ! new particles with a value for D (i.e. not allocate the array yet)
              ! do fuse/spawn first, then apply_bc, get_ghosts, etc, then alloc and
              ! compute D

              ! Small hack to speed up ghost_get and neighlists
              ! (for fusion/insertion of particles, we only need neighbours that are
              ! very close, so we can safely reduce the cutoff radii)
              ! The neighbour lists will anyway have to be recomputed later on
              ! Note that we cannot reduce the cutoffs of particles for which
              ! the ratio Dtilde/D is large.
              ! The factor 1.2 is the same as the one in sop_spawn (check_nn)
              !
              !rcp => Get_wps(Particles,Particles%rcp_id)
              !D      => Get_wps(Particles,Particles%D_id)
              !Dtilde => Get_wps(Particles,Particles%Dtilde_id)
              !DO ip=1,Particles%Npart
                  !rcp(ip) = rcp(ip) * &
                      !MIN(1._MK, 1.2_MK/opts%rcp_over_D * Dtilde(ip)/D(ip))
              !ENDDO
              !rcp => Set_wps(Particles,Particles%rcp_id)
              !D      => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
              !Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)
              !CALL particles_updated_cutoff(Particles,info)
              !IF (info .NE. 0) THEN
                  !info = ppm_error_error
                  !CALL ppm_error(ppm_err_sub_failed,caller,&
                      !'particles_updated_cutoff failed',__LINE__,info)
                  !GOTO 9999
              !ENDIF

              CALL particles_mapping_ghosts(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'particles_mapping_ghosts failed',__LINE__,info)
                  GOTO 9999
              ENDIF

              CALL particles_neighlists(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'particles_neighlists failed',__LINE__,info)
                  GOTO 9999
              ENDIF

              !call check_duplicates(Particles)

              CALL ppm_util_time(add_t1)

              !count neighbours at a distance < D and decide
              ! whether we need to add new particles
              adaptation_ok = .TRUE.
              CALL sop_close_neighbours(Particles,opts,info)

              !Insert (spawn) new particles where needed

              !if (lbfgs_continue .AND. gradPsi_max.LT.1e-1 .OR. &
                  !&  Particles%nneighmin.LT.opts%nneigh_critical) then

              Dtilde => Particles%wps(Particles%Dtilde_id)%vec
              DO ip=1,Particles%Npart
                  write(400+it_adapt,'(3(E20.12,2X))') Particles%xp(1:2,ip),Dtilde(ip)
              ENDDO
              Dtilde => NULL()

              IF (.not.adaptation_ok) then
#ifdef __USE_DEL_METHOD2
                  CALL  sop_spawn2_particles(Particles,opts,info,&
                      nb_part_added=nb_spawn,wp_fun=wp_fun)
#else
                  CALL  sop_spawn_particles(Particles,opts,info,&
                      nb_part_added=nb_spawn,wp_fun=wp_fun)
#endif
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,caller,&
                          'sop_spawn_particles failed',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  adding_particles = (nb_spawn .GT. 0)
                  Particles%neighlists=.FALSE.
#ifdef __USE_LBFGS
                  lbfgs_continue = .FALSE.
#endif
              else
                  adding_particles = .FALSE.
              ENDIF
              Dtilde => Particles%wps(Particles%Dtilde_id)%vec
              DO ip=1,Particles%Npart
                  write(500+it_adapt,'(3(E20.12,2X))') Particles%xp(1:2,ip),Dtilde(ip)
              ENDDO
              Dtilde => NULL()
              CALL check_duplicates(Particles)

#ifdef __USE_LBFGS
              IF (adding_particles) lbfgs_continue = .FALSE.
#endif

              IF (adding_particles) THEN
                  CALL particles_updated_positions(Particles,info)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,caller,&
                          'particles_updated_positions failed',__LINE__,info)
                      GOTO 9999
                  ENDIF

                  !---------------------------------------------------------------------!
                  ! Ensure that particles satisfy the boundary conditions
                  ! (Here because partial remapping fails otherwise).
                  !---------------------------------------------------------------------!
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

              ENDIF

              CALL particles_mapping_ghosts(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'particles_mapping_ghosts failed',__LINE__,info)
                  GOTO 9999
              ENDIF
              CALL particles_neighlists(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'particles_neighlists failed',__LINE__,info)
                  GOTO 9999
              ENDIF


              CALL ppm_util_time(add_t2)
              Particles%stats%t_add = Particles%stats%t_add + (add_t2-add_t1)


#if debug_verbosity > 2
              WRITE(filename,'(A,I0,A,I0)') 'P_SOP_afterSpawn_',&
                  Particles%itime,'_',it_adapt
              CALL ppm_vtk_particle_cloud(filename,Particles,info,&
                  wpi_list=(/nb_neigh_id,fuse_id/),&
                  wps_list=(/(i,i=1,0)/),wpv_list=(/(i,i=1,0)/))
#endif

              !Delete (fuse) particles that are too close to each other
              !(needs ghost particles to be up-to-date)

              CALL ppm_util_time(del_t1)

#ifdef __USE_DEL_METHOD2
              CALL sop_fuse2_particles(Particles,opts,info,nb_part_del=nb_fuse)
#else
              CALL sop_fuse_particles(Particles,opts,info,nb_part_del=nb_fuse)
#endif
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'sop_fuse_particles failed',__LINE__,info)
                  GOTO 9999
              ENDIF
              !we only removed particles - they didnt move.
              Particles%areinside=.TRUE.
              Particles%ontopology=.TRUE.

              CALL ppm_util_time(del_t2)
              Particles%stats%t_del = Particles%stats%t_del + (del_t2-del_t1)

#ifdef __USE_LBFGS
              IF (nb_fuse .GT. 0) lbfgs_continue = .FALSE.
#endif

              CALL particles_mapping_ghosts(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_sub_failed,caller,&
                      'particles_mapping_ghosts failed',__LINE__,info)
                  GOTO 9999
              ENDIF

              CALL ppm_util_time(compD_t1)

              Compute_D: IF (PRESENT(wp_grad_fun).OR. &
                  (.NOT.need_derivatives.AND.PRESENT(wp_fun))) THEN
                  !!-----------------------------------------------------------------!
                  !! Get D directly from a given function
                  !!-----------------------------------------------------------------!
                  CALL sop_compute_D(Particles,D_fun,opts,info,     &
                      wp_fun=wp_fun,wp_grad_fun=wp_grad_fun)
                  IF (info .NE. 0) THEN
                      info = ppm_error_error
                      CALL ppm_error(ppm_err_sub_failed,caller,&
                          'sop_compute_D failed',__LINE__,info)
                      GOTO 9999
                  ENDIF
              ELSE
                  ! do not update D

              ENDIF Compute_D

              CALL ppm_util_time(compD_t2)
              Particles%stats%t_compD = Particles%stats%t_compD + (compD_t2-compD_t1)

              IF (.NOT. PRESENT(wp_fun)) THEN
                  !------------------------------------------------------------------!
                  ! Update D and cutoff radii
                  !------------------------------------------------------------------!

                  !The following is used to prevent particles with small D from
                  ! drifting away inside the domain, eventually generating plenty
                  ! of other small particles that can potentially fill the whole
                  ! box (defeating the point of having an adaptive scheme...)
                  ! It means that a particle can have a small D only
                  ! if it has neighbours from the initial generation (D_old) that
                  ! also have a small D.

              CALL ppm_util_time(t1)

              !WRITE(filename,'(A,I0,A,I0)') 'debug_old_',Particles%itime,'_',it_adapt
              !CALL ppm_vtk_particle_cloud(filename,Particles_old,info)
              !WRITE(filename,'(A,I0,A,I0)') 'debug_new_',Particles%itime,'_',it_adapt
              !CALL ppm_vtk_particle_cloud(filename,Particles,info)

                  D_old => Get_wps(Particles_old,Particles_old%D_id)
                  D_old = D_old * opts%rcp_over_D
                  ghostlayer=Particles%cutoff
                  CALL ppm_inl_xset_vlist(topo_id,Particles%xp,Particles%Npart,&
                      Particles%Mpart,Particles_old%xp,Particles_old%Npart,&
                      Particles_old%Npart,D_old,Particles%skin,&
                      ghostlayer,info,vlist_cross,nvlist_cross)
                  D_old = D_old / opts%rcp_over_D
                  D_old => Set_wps(Particles_old,Particles_old%D_id,&
                      read_only=.TRUE.)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'ppm_inl_xset_vlist failed.',info)
                      info = -1
                      GOTO 9999
                  ENDIF
                  !write(*,*) 'MIN MAX XSET INL : ', &
                  ! MINVAL(nvlist_cross(1:Particles%Npart)), &
                      !MAXVAL(nvlist_cross(1:Particles%Npart))

#if debug_verbosity > 0
                  IF (it_adapt.EQ.1 .AND. &
                      MINVAL(nvlist_cross(1:Particles%Npart)).LE.0) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'Insufficient number of xset neighbours to compute D',info)
                      info = -1
                      GOTO 9999
                  ENDIF
#endif

                  Dtilde  => Get_wps(Particles,Particles%Dtilde_id)
                  D     => Get_wps(Particles,    Particles%D_id)
                  Dtilde_old => Get_wps(Particles_old,Particles_old%Dtilde_id)
                  D_old => Get_wps(Particles_old,Particles_old%D_id)
                  rcp => Get_wps(Particles,Particles%rcp_id)
                  xp => Get_xp(Particles)
                  xp_old => Get_xp(Particles_old,with_ghosts=.TRUE.)

              CALL ppm_util_time(t2)
              Particles%stats%t_xset_inl = Particles%stats%t_xset_inl + (t2-t1)
              Particles%stats%nb_xset_inl = Particles%stats%nb_xset_inl+1


              !Linear interpolation of D_tilde
                  DO ip=1,Particles%Npart
                      IF (nvlist_cross(ip).EQ.0) then
                          Dtilde(ip) = opts%maximum_D
                          D(ip) = opts%maximum_D
                      else
                          minDold=big
                          tmpvar1=big
                          Dtilde(ip) = 0.0_MK
                          weight_sum = 0.0_MK

                          DO ineigh=1,nvlist_cross(ip)
                              iq=vlist_cross(ineigh,ip)
                              minDold=MIN(minDold,Dtilde_old(iq))
                          ENDDO

                          n_loop: DO ineigh=1,nvlist_cross(ip)
                              iq=vlist_cross(ineigh,ip)

                              weight = SUM(((xp(1:ppm_dim,ip)-&
                                  xp_old(1:ppm_dim,iq))/D(ip))**2)! ** 2
                                  !xp_old(1:ppm_dim,iq))/Dtilde_old(iq))**2)! ** 2
                              IF (weight .LT. almostzero) THEN
                                  Dtilde(ip) = Dtilde_old(iq)
                                  weight_sum = 1.0_MK
                                  EXIT n_loop
                              ELSE
                                  weight_sum = weight_sum + 1._MK / weight
                                  Dtilde(ip) = Dtilde(ip) + Dtilde_old(iq) / weight
                              ENDIF

                          ENDDO n_loop

                          Dtilde(ip) = Dtilde(ip) / weight_sum
                          !D(ip) = MAX(D(ip),minDold)
                          !D(ip) = minDold

                          D(ip) = Dtilde(ip)

                          !n_loop2: DO ineigh=1,nvlist_cross(ip)
                              !iq=vlist_cross(ineigh,ip)

                              !weight = SQRT(SUM((xp(1:ppm_dim,ip)-&
                                  !xp_old(1:ppm_dim,iq)**2)))/Dtilde_old(iq) -1._MK
                              !IF (weight .LT. 0._MK) THEN
                                  !D(ip) = MIN(D(ip),Dtilde_old(iq))
                              !ELSE
                                  !D(ip) = MIN(D(ip),2._MK**(weight)*Dtilde_old(iq))
                              !ENDIF

                          !ENDDO n_loop2
                          !D(ip) = Dtilde(ip)

                      ENDIF

                      !when Dtilde/D is large, increase rcp
                      rcp(ip) = opts%rcp_over_D * MIN(Dtilde(ip),2._MK*D(ip))
                      !rcp(ip) = opts%rcp_over_D * MIN(D(ip),2._MK*minDold)
                      !rcp(ip) = opts%rcp_over_D * D(ip)
                  ENDDO
                  D     => Set_wps(Particles,    Particles%D_id)
                  D_old => Set_wps(Particles_old,Particles_old%D_id,&
                      read_only=.TRUE.)
                  Dtilde_old => Set_wps(Particles_old,Particles_old%Dtilde_id,&
                      read_only=.TRUE.)
                  rcp => Set_wps(Particles,Particles%rcp_id)
                  Dtilde  => Set_wps(Particles,Particles%Dtilde_id)
                  xp => Set_xp(Particles,read_only=.TRUE.)
                  xp_old => Set_xp(Particles_old,read_only=.TRUE.)
              ELSE
                  !------------------------------------------------------------------!
                  ! Update cutoff radii
                  !------------------------------------------------------------------!
                  D => Get_wps(Particles,Particles%D_id)
                  Dtilde => Get_wps(Particles,Particles%Dtilde_id)
                  rcp => Get_wps(Particles,Particles%rcp_id)
                  DO ip=1,Particles%Npart
                      !so that particles near the edge between regions of high
                      ! and low resolutions get the right number of neighbours
                      rcp(ip) = opts%rcp_over_D * MIN(Dtilde(ip),2._MK*D(ip))
                      !rcp(ip) = opts%rcp_over_D * D(ip)
                  ENDDO
                  Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)
                  D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
                  rcp => Set_wps(Particles,Particles%rcp_id)
              ENDIF

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

#if debug_verbosity>1
              CALL check_duplicates(Particles)
#endif

              !!---------------------------------------------------------------------!
              !! /begin Line search **
              !!---------------------------------------------------------------------!
              Particles%stats%nb_ls = Particles%stats%nb_ls + 1

              CALL ppm_util_time(ls_t1)

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

#if debug_verbosity > 1
              CALL particles_allocate_wps(Particles,potential_before_id,info,&
                  iopt=ppm_param_alloc_fit,name='potential_before')
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'allocation failed',info)
                  info = -1
                  GOTO 9999
              ENDIF
#endif

#ifdef __USE_LBFGS
              !L-BFGS

              CALL sop_gradient_psi(Particles,topo_id,Gradient_Psi,Psi_global,&
                  Psi_max,opts,info,gradPsi_max=gradPsi_max)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'sop_gradient_psi failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF

              !call particles_check_arrays(Particles,info)
              !IF (info .NE. 0) THEN
                  !CALL ppm_write(ppm_rank,caller,'particles_check_arrays failed.',info)
                  !info = -1
                  !GOTO 9999
              !ENDIF
              !call check_duplicates(Particles)
              CALL particles_check_arrays(Particles,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_check_arrays before failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF

              IF ( lbfgs_continue ) THEN
                  info = 1
              ELSE
                  IF (ALLOCATED(Work)) DEALLOCATE(Work)
                  IF (ALLOCATED(DIAG)) DEALLOCATE(DIAG)
                  ALLOCATE(Work(3*Particles%Mpart*(2*3+1)+2*3))
                  ALLOCATE(DIAG(3*Particles%Mpart))
                  Work=0._MK
                  DIAG=1._MK
                  info = 0
              ENDIF
              rcp => Get_wps(Particles,potential_before_id)
              Dtilde => Get_wps(Particles,Particles%Dtilde_id)
              DO ip = 1,Particles%Npart
                  rcp(ip) = SQRT(SUM(Gradient_Psi(1:ppm_dim,ip)**2))/Dtilde(ip)
              ENDDO
              rcp => Set_wps(Particles,potential_before_id)
              Dtilde => Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)


              !IF (lbfgs_continue .AND. gradPsi_max .LE. 5E-2) THEN
              IF (lbfgs_continue .AND. gradPsi_max .LE. 5E-12) THEN
                  adaptation_ok = .TRUE.
              ELSE

                  write(*,*) 'before LBFGS: info = ',info, Particles%Npart
                  Dtilde=>Get_wps(Particles,Particles%Dtilde_id)
                  CALL LBFGS(ppm_dim*Particles%Npart,3,&
                      Particles%xp(1:ppm_dim,1:Particles%Npart),&
                      Psi_global,Gradient_Psi(1:ppm_dim,1:Particles%Npart),&
                      .FALSE.,DIAG,(/1,0/),1D-3,ppm_myepsd,Work,info,scaling=Dtilde)
                  Dtilde=>Set_wps(Particles,Particles%Dtilde_id,read_only=.TRUE.)
                  write(*,*) 'after LBFGS: info = ',info
                  IF (info.LT.0) THEN
                      lbfgs_continue = .FALSE.

                      adaptation_ok = .TRUE.
                      CALL sop_close_neighbours(Particles,opts,info)

                      IF (adaptation_ok) THEN
                          CALL ppm_write(ppm_rank,caller,&
                              'LBFGS failed, but adaptation is ok',info)
                      ELSE
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_sub_failed,caller,&
                              'LBFGS failed', __LINE__,info)
                          GOTO 9999
                      ENDIF
                  ELSE IF (info.GT.0) THEN
                      lbfgs_continue = .TRUE.
                      adaptation_ok = .FALSE.
                  ELSE
                      IF (lbfgs_continue) THEN
                          write(*,*) 'found suitable minimum in LBFGS'
                          adaptation_ok = .TRUE.
                      ENDIF
                  ENDIF
              ENDIF

              CALL particles_check_arrays(Particles,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_check_arrays after failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF


              !CALL sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)
              !IF (info .NE. 0) THEN
                  !CALL ppm_write(ppm_rank,caller,'sop_potential_psi failed.',info)
                  !info = -1
                  !GOTO 9999
              !ENDIF


#elif defined __USE_SD

              !STEEPEST DESCENT

              !!---------------------------------------------------------------------!
              !! Compute gradient of the potential
              !! (need ghosts for xp and D)
              !! (on output, the ghost values for Gradient_Psi have been updated)
              !!---------------------------------------------------------------------!
              CALL sop_gradient_psi(Particles,topo_id,Gradient_Psi,Psi_global,&
                  Psi_max,opts,info,gradPsi_max=gradPsi_max)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'sop_gradient_psi failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF


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
#if debug_verbosity > 1
              CALL sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)
              IF (Psi_global .NE. Psi_global_old) THEN
                  write(*,*) 'global potential different in ', &
                      'sop_gradient_psi and sop_potential_psi'
                  write(*,*) 'Psi_global(sop_potential) = ',Psi_global
                  write(*,*) 'Psi_global(sop_gradient) = ',Psi_global_old
                  write(*,*) 'difference = ', Psi_global_old - Psi_global
                  info = -1
                  GOTO 9999
              ENDIF

              CALL sop_voronoi_MC(Particles,opts,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'sop_voronoi_MC failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF
              WRITE(filename,'(A,I0,A,I0)') 'P_Voronoi_It_',&
                  Particles%itime,'_',it_adapt
              CALL ppm_vtk_particle_cloud(filename,Particles,info)

#endif


              Psi_1 = big
              alpha1 = -1._MK
              alpha2 = -1._MK
              step_previous = 0._MK

              !IF (gradPsi_max .LT. 1e-3) then
                  !step_max = 10.0_MK
              !else
                  !step_max = 1.5_MK
              !endif
              step_max = MIN(0.3_MK / MIN(gradPsi_max,3.0_MK), 10000.0_MK)
              !step_max = MIN(0.3_MK / MIN(gradPsi_max,3.0_MK), 0.5_MK)

              !!---------------------------------------------------------------------!
              !! Evaluate potential after different step sizes
              !!---------------------------------------------------------------------!
              linesearch_loop: DO WHILE (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK)
                  xp => Get_xp(Particles,with_ghosts=.TRUE.)
                  DO ip=1,Particles%Mpart
                      xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) - &
                          (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
                  ENDDO
                  xp => Set_xp(Particles,ghosts_ok=.TRUE.)
                  step_previous = step

                  CALL sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,'sop_potential_psi failed.',info)
                      info = -1
                      GOTO 9999
                  ENDIF

                  IF (Psi_global .LT. Psi_global_old) THEN
                      IF (Psi_global .LT. Psi_1) THEN
                          Psi_1 = Psi_global
                          alpha1 = step
                      ENDIF
                      step = 2._MK * step
                      IF (step .GT. step_max) EXIT linesearch_loop
                  ELSE IF (Psi_global .GT. Psi_global_old) THEN
                      Psi_2 = Psi_global
                      alpha2 = step
                      step = 0.5_MK * step
                      IF (step .LT. step_min) EXIT linesearch_loop
                  ELSE
                      step = 10._MK * step
                      IF (step .GT. step_max) EXIT linesearch_loop
                  ENDIF

                  !FIXME
                  IF (ABS(Psi_global_old - Psi_global)/Psi_global_old .LT. 1E-4) &
                      EXIT linesearch_loop

              ENDDO linesearch_loop

              1000 CONTINUE

              IF (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK) THEN
                  step = MAX(step_min,MIN(step,step_max))
              ELSE
                  !Quadratic fit
                  Psi_1 = Psi_1 - Psi_global_old
                  Psi_2 = Psi_2 - Psi_global_old
                  step = 0.5_MK * (alpha1**2 * Psi_2 - alpha2**2 * Psi_1) / &
                      (alpha1*Psi_2 - alpha2*Psi_1)
              ENDIF

              !!---------------------------------------------------------------------!
              !! Choose best step size (from quadratic fit)
              !! Move particles (no need to move the ghosts, since we will have
              !! to get them through a local mapping anyway...)
              !!---------------------------------------------------------------------!
              !Move particles (including ghosts)
              xp => Get_xp(Particles)
              DO ip=1,Particles%Npart
                  xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) - &
                      (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
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

#else
              WRITE(*,*) 'This routine needs to be compiled with either '
              WRITE(*,*) '__USE_LBFGS or __USE_SD precompiler flags'
              info = -1
              GOTO 9999
#endif
      !end ifdef between LBFGS and  SD  algorithms

              CALL ppm_util_time(ls_t2)
              Particles%stats%t_ls = Particles%stats%t_ls + (ls_t2-ls_t1)


#if debug_verbosity > 0
#ifdef __MPI
              CALL MPI_Allreduce(Psi_max,tmpvar2,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
              CALL MPI_Allreduce(Particles%Npart,tmpvari1,1,&
                  MPI_INTEGER,MPI_SUM,ppm_comm,info)
              CALL MPI_Allreduce(Particles%Mpart,tmpvari2,1,&
                  MPI_INTEGER,MPI_SUM,ppm_comm,info)
              IF (ppm_rank.EQ.0) THEN
                  WRITE(cbuf,'(A,I0,A,E19.10,A,E11.4,A,I4,1X,I4,1X,A,I6,A,I6,A,E9.2)') &
                      'it_adapt= ',it_adapt,&
                      ' V= ',Psi_global, & !/REAL(tmpvari1,MK),
                      ' Psi_max= ',tmpvar2, &
                      ' Nneigh= ', Particles%nneighmin, Particles%nneighmax, &
                      'Np=',tmpvari1,' Mp=',tmpvari2,' step=',step
                  CALL ppm_write(ppm_rank,caller,cbuf,info)
#if debug_verbosity > 2
                  CALL particles_print_stats(Particles,info)
#endif
              ENDIF
#endif
#endif
#ifdef __MPI
              CALL MPI_Allreduce(adding_particles,adding_particles,1,&
                  MPI_LOGICAL,MPI_LOR,ppm_comm,info)
#endif

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

#if debug_verbosity > 1
              WRITE(filename,'(A,I0,A,I0)') 'P_SOP_It_',&
                  Particles%itime,'_',it_adapt
              CALL ppm_vtk_particle_cloud(filename,Particles,info)
#endif


          ENDDO it_adapt_loop


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
#ifdef __USE_LBFGS
          IF (ALLOCATED(Work)) DEALLOCATE(Work)
          IF (ALLOCATED(DIAG)) DEALLOCATE(DIAG)
#endif
#if debug_verbosity > 1
              CALL particles_allocate_wps(Particles,potential_before_id,info,&
                  iopt=ppm_param_dealloc)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'deallocation failed',info)
                  info = -1
                  GOTO 9999
              ENDIF
#endif

#if debug_verbosity < 1
          CALL particles_allocate_wpi(Particles,fuse_id,info,&
              iopt=ppm_param_dealloc)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,&
                  'particles_allocate_wpi (dealloc) failed', __LINE__,info)
              GOTO 9999
          ENDIF
          CALL particles_allocate_wpi(Particles,nb_neigh_id,info,&
              iopt=ppm_param_dealloc)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,caller,&
                  'particles_allocate_wpi (dealloc) failed', __LINE__,info)
              GOTO 9999
          ENDIF
#endif

#if debug_verbosity > 0
          CALL substop(caller,t0,info)
#endif

          9999 CONTINUE ! jump here upon error

      END SUBROUTINE DTYPE(sop_gradient_descent)




      !!!----------------------------------------------------------------------------!
      !!! LEVEL SET implementation of:
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

      SUBROUTINE DTYPE(sop_gradient_descent_ls)(Particles_old,Particles, &
              nvlist_cross,vlist_cross,   &
              nneighmin_cross,nneighmax_cross,num_it,opts,info, &
              wp_fun,D_fun,wp_grad_fun,level_fun,level_grad_fun,&
              threshold,need_deriv,nb_fun,stats)

          USE ppm_module_mpi
          USE ppm_module_inl_xset_vlist
          USE ppm_module_sop_typedef
          IMPLICIT NONE

          DEFINE_MK()
          ! arguments
          TYPE(DTYPE(ppm_t_particles)), POINTER, INTENT(IN   )   :: Particles_old
          TYPE(DTYPE(ppm_t_particles)), POINTER, INTENT(INOUT)   :: Particles
          INTEGER,      DIMENSION(:),  POINTER,  INTENT(INOUT)   :: nvlist_cross
          INTEGER,      DIMENSION(:,:),POINTER,  INTENT(INOUT)   :: vlist_cross
          INTEGER,                               INTENT(INOUT)   :: nneighmax_cross
          INTEGER,                               INTENT(INOUT)   :: nneighmin_cross
          INTEGER,                               INTENT(  OUT)   :: num_it
          TYPE(DTYPE(sop_t_opts)), POINTER,      INTENT(IN   )   :: opts
          INTEGER,                               INTENT(  OUT)   :: info

          !optional arguments
          REAL(MK), OPTIONAL,                    INTENT(IN)      :: threshold
          LOGICAL,  OPTIONAL,                    INTENT(IN)      :: need_deriv
          OPTIONAL                                               :: D_fun
          !Monitor function
          OPTIONAL                                               :: wp_fun
          !Field function (usually known only during initialisation)
          OPTIONAL                                               :: wp_grad_fun
          !Gradient of the field function (usually known only during initialisation)
          OPTIONAL                                               :: level_fun
          !Level function (usually known only during initialisation)
          OPTIONAL                                               :: level_grad_fun
          !Gradient of the level function (usually known only during initialisation)
          ! argument-functions need an interface
          TYPE(DTYPE(sop_t_stats)),POINTER,OPTIONAL,INTENT(OUT)  :: stats
          !!! statistics on output
          INTERFACE
              FUNCTION D_fun(f1,dfdx,opts,f2)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  USE ppm_module_sop_typedef
                  DEFINE_MK()
                  REAL(MK)                               :: D_fun
                  REAL(MK),                   INTENT(IN) :: f1
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN) :: dfdx
                  TYPE(DTYPE(sop_t_opts)),POINTER,INTENT(IN) :: opts
                  REAL(MK),OPTIONAL,          INTENT(IN) :: f2
              END FUNCTION D_fun

              !Function that returns the width of the narrow band
              FUNCTION nb_fun(kappa,scale_D)
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK)                             :: nb_fun
                  REAL(MK),                INTENT(IN)  :: kappa
                  REAL(MK),                INTENT(IN)  :: scale_D
              END FUNCTION nb_fun

              FUNCTION wp_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
                  REAL(MK)                                :: wp_fun
              END FUNCTION wp_fun

              FUNCTION wp_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
                  REAL(MK),DIMENSION(ppm_dim)             :: wp_grad_fun
              END FUNCTION wp_grad_fun

              FUNCTION level_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
                  REAL(MK)                                :: level_fun
              END FUNCTION level_fun

              FUNCTION level_grad_fun(pos)
                  USE ppm_module_data, ONLY: ppm_dim
                  USE ppm_module_typedef
                  DEFINE_MK()
                  REAL(MK),DIMENSION(ppm_dim),INTENT(IN)  :: pos
                  REAL(MK),DIMENSION(ppm_dim)             :: level_grad_fun
              END FUNCTION level_grad_fun

          END INTERFACE

          ! local variables
          INTEGER                             :: ip, it_adapt,iq,ineigh,di
          INTEGER                             :: iunit
          CHARACTER(LEN = 256)                :: filename,cbuf
          CHARACTER(LEN = 256)                :: caller='sop_gradient_descent_ls'
          REAL(ppm_kind_double) :: t0

          REAL(MK)                            :: step,step_previous,alpha1,alpha2
          REAL(MK)                            :: step_min,step_max,step_stall
          INTEGER                             :: it_adapt_max
          REAL(MK)                            :: Psi_max,Psi_global,Psi_global_old
          REAL(MK)                            :: Psi_1,Psi_2,Psi_thresh
          REAL(MK),DIMENSION(ppm_dim)         :: dist,dist2,dummy_grad
          REAL(MK),DIMENSION(:,:),POINTER     :: Gradient_Psi => NULL()
          INTEGER                             :: nneigh_adapt

          REAL(MK),DIMENSION(:,:),POINTER     :: xp => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: D => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: rcp => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: Dtilde => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: level => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: wp => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: wp_old => NULL()
          REAL(MK),DIMENSION(:,:),POINTER     :: level_grad => NULL()
          REAL(MK),DIMENSION(:,:),POINTER     :: xp_old => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: D_old => NULL()
          REAL(MK),DIMENSION(:),POINTER       :: level_old => NULL()
          REAL(MK),DIMENSION(:,:),POINTER     :: level_grad_old => NULL()

          REAL(MK)                            :: tmpvar1,tmpvar2,minDold
          REAL(MK)                            :: weight,weight_sum,lev
          REAL(MK)                            :: almostzero,nb
          INTEGER                             :: tmpvari1,tmpvari2
          LOGICAL                             :: need_derivatives
          INTEGER                             :: topo_id

          !should be removed once the argument lists for the inl routines
          !have been updated to inhomogeneous ghostlayers
          REAL(MK),DIMENSION(2*ppm_dim)       :: ghostlayer

          !!-------------------------------------------------------------------------!
          !! Initialize
          !!-------------------------------------------------------------------------!
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

          it_adapt = 0
          Psi_max = big
          Psi_global = big
          Psi_global_old = big

          ! Adaptivity stopping criterion
          IF (PRESENT(threshold)) THEN
              Psi_thresh = threshold
          ELSE
              Psi_thresh = 2.5_MK
          ENDIF

          IF (MK .EQ. KIND(1.D0)) THEN
              !step_min = 1D-8; step_stall = 1D-14
              step_min = 1D-3; step_stall = 1D-14
          ELSE
              !step_min = 1E-8; step_stall = 1E-14
              step_min = 1E-3; step_stall = 1E-14
          ENDIF

          step_max = 0.1_MK
          it_adapt_max = 9999

          step = 1._MK
          nneigh_adapt = 20 !opts%nneigh_theo

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
          it_adapt_loop: DO WHILE (step .GT. step_stall .AND. &
                  &          it_adapt .LT. it_adapt_max .AND. &
                  &     ((Psi_max .GT. Psi_thresh) .OR.  &
                  Particles%nneighmin .LT. nneigh_adapt))

              it_adapt = it_adapt + 1

              !!---------------------------------------------------------------------!
              !! Update number of particles by fusion/insertion
              !! of particles (need vlist to be up-to-date)
              !! (on output, Mpart becomes meaningless.)
              !!---------------------------------------------------------------------!
              !NOTE: if we know D_tilde analytically, we should not supply
              ! new particles with a value for D (i.e. not allocate the array yet)
              ! do fuse/spawn first, then apply_bc, get_ghosts, etc, then alloc and
              ! compute D
              !
              !FIXME: in theory, we should do apply_bc+remap+get_ghosts here
              ! (Fusing particles requires knowing the ghosts and particles
              ! have moved after the linesearch in the previous iteration)
              ! There MAY be a better way of doing this...
              !!---------------------------------------------------------------------!
              CALL particles_apply_bc(Particles,topo_id,info)
              CALL particles_mapping_partial(Particles,topo_id,info)
              CALL particles_mapping_ghosts(Particles,topo_id,info)
              CALL particles_neighlists(Particles,topo_id,info)


              !Delete (fuse) particles that are too close to each other
              !(needs ghost particles to be up-to-date)
              CALL sop_fuse_particles(Particles,opts,info,wp_fun=wp_fun,&
                  level_fun=level_fun,nb_fun=nb_fun,printp=it_adapt)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'sop_fuse_particles failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF

              !Insert (spawn) new particles where needed
              CALL  sop_spawn_particles(Particles,opts,info,wp_fun=wp_fun,&
                  level_fun=level_fun,nb_fun=nb_fun)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,'sop_spawn_particles failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF

              CALL particles_updated_positions(Particles,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_updated_positions failed',info)
                  info = -1
                  GOTO 9999
              ENDIF

              !!---------------------------------------------------------------------!
              !! Ensure that particles satisfy the boundary conditions
              !! (Here because partial remapping fails otherwise).
              !!---------------------------------------------------------------------!
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

              CALL particles_mapping_ghosts(Particles,topo_id,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_apply_bc failed',info)
                  info = -1
                  GOTO 9999
              ENDIF


              Compute_D: IF (PRESENT(wp_grad_fun).OR. &
                  (.NOT.need_derivatives.AND.PRESENT(wp_fun))) THEN
                  !!-----------------------------------------------------------------!
                  !! Get D directly from a given function
                  !!-----------------------------------------------------------------!
                  CALL sop_compute_D(Particles,D_fun,opts,info,     &
                      wp_fun=wp_fun,wp_grad_fun=wp_grad_fun,level_fun=level_fun,&
                      level_grad_fun=level_grad_fun,nb_fun=nb_fun)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,'sop_compute_D failed',info)
                      info = -1
                      GOTO 9999
                  ENDIF
              ELSE
                  ! do not update D
              ENDIF Compute_D

              !interpolate some quantities on new particles
              ! TODO: not very efficient if only few particles have been added
              IF (.NOT. PRESENT(wp_fun)) THEN

                  !The following is used to prevent particles with small D from
                  ! drifting away inside the domain, eventually generating plenty
                  ! of other small particles that can potentially fill the whole
                  ! box (defeating the point of having an adaptive scheme...)
                  ! Roughly, it means that a particle can have a small D only
                  ! if it has neighbours from the older generation (D_old) that
                  ! also have a small D.
                  D_old => Get_wps(Particles_old,Particles_old%D_id)
                  D_old = D_old * opts%rcp_over_D
                  ghostlayer=Particles%cutoff
                  CALL ppm_inl_xset_vlist(topo_id,Particles%xp,Particles%Npart,&
                      Particles%Mpart,Particles_old%xp,Particles_old%Npart,&
                      Particles_old%Npart,D_old,Particles%skin,&
                      ghostlayer,info,vlist_cross,nvlist_cross)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'ppm_inl_xset_vlist failed.',info)
                      info = -1
                      GOTO 9999
                  ENDIF

#if debug_verbosity > 0
                  IF (MINVAL(nvlist_cross(1:Particles%Npart)).EQ.0) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'xlist empty dumping some debug data',info)
                      DO ip=1,Particles%Npart
                          write(801,*) Particles%xp(1:ppm_dim,ip), nvlist_cross(ip)
                      ENDDO
                      DO ip=1,Particles_old%Npart
                          write(802,*) Particles_old%xp(1:ppm_dim,ip), &
                              opts%rcp_over_D*D_old(ip)
                      ENDDO
                      info = -1
                      GOTO 9999
                  ENDIF
#endif

                  D_old = D_old / opts%rcp_over_D
                  D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)

                  D     => Get_wps(Particles,    Particles%D_id)
                  D_old => Get_wps(Particles_old,Particles_old%D_id)
                  DO ip=1,Particles%Npart
                      minDold=big
                      DO ineigh=1,nvlist_cross(ip)
                          iq=vlist_cross(ineigh,ip)
                          minDold=MIN(minDold,D_old(iq))
                      ENDDO
                      D(ip) = MAX(D(ip),minDold)
                  ENDDO
                  D     => Set_wps(Particles,    Particles%D_id)
                  D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)

                  wp     => Get_wps(Particles,Particles%adapt_wpid)
                  wp_old => Get_wps(Particles_old,Particles_old%adapt_wpid)
                  D_old => Get_wps(Particles_old,Particles_old%D_id)
                  CALL sop_approx_wp(Particles_old%xp,wp_old,D_old,Particles%xp,&
                      wp,Particles%Npart,Particles%Mpart,nvlist_cross,vlist_cross,info)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'sop_approx_wp failed.',info)
                      info = -1
                      GOTO 9999
                  ENDIF
                  D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)
                  wp_old=> Set_wps(Particles_old,Particles_old%adapt_wpid,&
                      read_only=.TRUE.)
                  wp    => Set_wps(Particles,Particles%adapt_wpid)

                  IF (Particles%level_id.EQ.0) THEN
                          !removme
                          CALL ppm_write(ppm_rank,caller,&
                              'NOT sure this can ever happen. &
                              &level_id should not be zero here',info)
                          info = -1
                          GOTO 9999
                          !removme
                  ELSE
                      level => Get_wps(Particles,Particles%level_id)
                      level_grad => Get_wpv(Particles,Particles%level_grad_id)
                      level_old => Get_wps(Particles_old,Particles_old%level_id)
                      level_grad_old => Get_wpv(Particles_old,&
                          Particles_old%level_grad_id)

                      !! FIXME: these routines also get the ghosts
                      !! (either the state of Particles should be updated
                      !! accordingly, or this should be moved outside of the
                      !! routine)
                      D     => Get_wps(Particles,    Particles%D_id)
                      D_old => Get_wps(Particles_old,Particles_old%D_id)
                      CALL sop_approx_wp(Particles_old%xp,D_old,D_old,Particles%xp,&
                          D,Particles%Npart,Particles%Mpart,&
                          nvlist_cross,vlist_cross,info)
                      IF (info .NE. 0) THEN
                          CALL ppm_write(ppm_rank,caller,&
                              'sop_approx_wp failed.',info)
                          info = -1
                          GOTO 9999
                      ENDIF

                      !  CALL sop_approx_wp_1d(Particles_old%xp,level_old,Particles%xp,&
                      !      level,Particles%Npart,Particles%Mpart,&
                      !      nvlist_cross,vlist_cross,info)

                      CALL sop_approx_wp(Particles_old%xp,level_grad_old,ppm_dim,D_old,&
                          Particles%xp,level_grad,Particles%Npart,Particles%Mpart,&
                          nvlist_cross,vlist_cross,info)
                      IF (info .NE. 0) THEN
                          CALL ppm_write(ppm_rank,caller,&
                              'sop_approx_wp failed.',info)
                          info = -1
                          GOTO 9999
                      ENDIF
                      D => Set_wps(Particles,    Particles%D_id)
                      D_old => Set_wps(Particles_old,Particles_old%D_id,read_only=.TRUE.)

                      xp => Get_xp(Particles,with_ghosts=.TRUE.)
                      xp_old => Get_xp(Particles_old)
                      ploop: DO ip=1,Particles%Npart
                          level(ip)  = 0._MK
                          weight_sum = 0._MK
                          DO ineigh = 1,nvlist_cross(ip)
                              iq = vlist_cross(ineigh,ip)
                              dist = xp(1:ppm_dim,ip) - xp_old(1:ppm_dim,iq)
                              weight = SUM(dist**2) ** 2
                              IF (weight .LT. almostzero) THEN
                                  level(ip) = level_old(iq)
                                  CYCLE ploop
                              ELSE
                                  weight_sum = weight_sum + 1._MK / weight
                                  level(ip) = level(ip) + (level_old(iq) + &
                                      SUM(dist*level_grad_old(1:ppm_dim,iq))) / weight
                              ENDIF
                          ENDDO
                          level(ip) = level(ip) / weight_sum
                      ENDDO ploop
                      xp => Set_xp(Particles,read_only=.TRUE.)
                      xp_old => Set_xp(Particles_old,read_only=.TRUE.)

                      level => Set_wps(Particles,Particles%level_id)
                      level_grad => Set_wpv(Particles,Particles%level_grad_id)
                      level_old => Set_wps(Particles_old,&
                          Particles_old%level_id,read_only=.TRUE.)
                      level_grad_old => Set_wpv(Particles_old,&
                          Particles_old%level_grad_id,read_only=.TRUE.)
                  ENDIF


              ENDIF

              !---------------------------------------------------------------------!
              ! Update cutoff radii
              !---------------------------------------------------------------------!
              D => Get_wps(Particles,Particles%D_id)
              rcp => Get_wps(Particles,Particles%rcp_id)
              DO ip=1,Particles%Npart
                  rcp(ip) = opts%rcp_over_D * D(ip)
              ENDDO

                  !HACK (trying this to make sure particles near the edge of the
                  !narrow band get enough neighbours)
                  IF (Particles%level_id.EQ.0) THEN
                      IF (PRESENT(level_fun)) THEN
                          xp => Get_xp(Particles)
                          DO ip=1,Particles%Npart
                              IF (ABS(level_fun(xp(1:ppm_dim,ip))).GT. &
                      opts%nb_width2*nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)) THEN
                                  rcp(ip) = 1.3_MK * rcp(ip)
                              ENDIF
                          ENDDO
                          xp => Set_xp(Particles,read_only=.TRUE.)
                      ELSE
                          write(*,*) 'need to provide either level_id or wp_fun'
                          info = -1
                          GOTO 9999
                      ENDIF
                  ELSE
                      level => Get_wps(Particles,Particles%level_id)
                      wp    => Get_wps(Particles,Particles%adapt_wpid)
                      DO ip=1,Particles%Npart
                          IF (ABS(level(ip)).GT. &
                              opts%nb_width2*nb_fun(wp(ip),opts%scale_D)) THEN
                              rcp(ip) = 1.3_MK * rcp(ip)
                          ENDIF
                      ENDDO
                      level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
                      wp    => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
                  ENDIF

              D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
              rcp => Set_wps(Particles,Particles%rcp_id)

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

              ! add narrow band potential
              !D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
              IF (Particles%level_id.EQ.0) THEN
                  IF (.NOT. PRESENT(level_fun).OR. .NOT. PRESENT(level_grad_fun)) THEN
                      CALL ppm_write(ppm_rank,caller,&
                          'need both level_grad_fun and level_fun here',info)
                      info = -1
                      GOTO 9999
                  ENDIF
                  xp => Get_xp(Particles)
                  DO ip=1,Particles%Npart
                      nb = opts%nb_width*nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)
                      lev = level_fun(xp(1:ppm_dim,ip))
                      IF (ABS(lev).GT.nb) THEN
                          Psi_global= Psi_global + opts%param_nb*(ABS(lev)-nb)**2
                      ENDIF
                  ENDDO
                  DO ip=1,Particles%Mpart
                      nb = opts%nb_width*nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)
                      lev = level_fun(xp(1:ppm_dim,ip))
                      IF (ABS(lev).GT.nb) THEN
                          IF (lev .GT.0) THEN
                              Gradient_Psi(1:ppm_dim,ip) = Gradient_Psi(1:ppm_dim,ip) &
                                  - 2._MK*opts%param_nb*(lev-nb)&
                                  * level_grad_fun(xp(1:ppm_dim,ip))
                          ELSE
                              Gradient_Psi(1:ppm_dim,ip) = Gradient_Psi(1:ppm_dim,ip) &
                                  - 2._MK*opts%param_nb*(lev+nb)&
                                  * level_grad_fun(xp(1:ppm_dim,ip))
                          ENDIF
                      ENDIF
                  ENDDO
                  xp => Set_xp(Particles,read_only=.TRUE.)

              ELSE

                  level => Get_wps(Particles,Particles%level_id,with_ghosts=.TRUE.)
                  level_grad => Get_wpv(Particles,Particles%level_grad_id,&
                      with_ghosts=.TRUE.)
                  wp    => Get_wps(Particles,Particles%adapt_wpid,with_ghosts=.TRUE.)
                  DO ip=1,Particles%Npart
                      nb = opts%nb_width*nb_fun(wp(ip),opts%scale_D)
                      IF (ABS(level(ip)) .GT. nb) THEN
                          Psi_global= Psi_global + &
                              opts%param_nb* (ABS(level(ip))-nb)**2
                      ENDIF
                  ENDDO
                  DO ip=1,Particles%Mpart
                      nb = opts%nb_width*nb_fun(wp(ip),opts%scale_D)
                      IF (ABS(level(ip)) .GT. nb) THEN
                          IF (level(ip) .GT.0) THEN
                              Gradient_Psi(1:ppm_dim,ip)=Gradient_Psi(1:ppm_dim,ip) &
                                  - 2._MK*opts%param_nb* (level(ip)-nb) * &
                                  level_grad(1:ppm_dim,ip)
                          ELSE
                              Gradient_Psi(1:ppm_dim,ip)=Gradient_Psi(1:ppm_dim,ip) &
                                  - 2._MK*opts%param_nb* (level(ip)+nb) * &
                                  level_grad(1:ppm_dim,ip)
                          ENDIF
                      ENDIF
                  ENDDO
                  level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
                  level_grad => Set_wpv(Particles,&
                      Particles%level_grad_id,read_only=.TRUE.)
                  wp    => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
              ENDIF
              !D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)


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
              Psi_1 = big
              alpha1 = -1._MK
              alpha2 = -1._MK
              step_previous = 0._MK

              !!---------------------------------------------------------------------!
              !! Evaluate potential after different step sizes
              !!---------------------------------------------------------------------!
              linesearch_loop: DO WHILE (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK)

                  !move particles along the gradient with the current step size
                  xp => Get_xp(Particles,with_ghosts=.TRUE.)
                  DO ip=1,Particles%Mpart
                      xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + &
                          (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
                  ENDDO
                  xp => Set_xp(Particles,ghosts_ok=.TRUE.)
                  step_previous = step

                  !re-compute potential
                  D => Get_wps(Particles,Particles%D_id,with_ghosts=.TRUE.)
                  CALL sop_potential_psi(Particles,Psi_global,Psi_max,opts,info)
                  IF (info .NE. 0) THEN
                      CALL ppm_write(ppm_rank,caller,'sop_potential_psi failed.',info)
                      info = -1
                      GOTO 9999
                  ENDIF
                  D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)

                  ! add narrow band potential
                  IF (Particles%level_id.EQ.0) THEN
                      xp => Get_xp(Particles)
                      DO ip=1,Particles%Npart
                          lev = level_fun(xp(1:ppm_dim,ip))
                          nb = opts%nb_width*nb_fun(wp_fun(xp(1:ppm_dim,ip)),opts%scale_D)
                          IF (ABS(lev).GT. nb ) THEN
                              Psi_global= Psi_global + opts%param_nb*(ABS(lev)-nb)**2
                          ENDIF
                      ENDDO
                      xp => Set_xp(Particles,read_only=.TRUE.)
                  ELSE
                      level => Get_wps(Particles,Particles%level_id)
                      wp    => Get_wps(Particles,Particles%adapt_wpid)
                      DO ip=1,Particles%Npart
                          nb = opts%nb_width*nb_fun(wp(ip),opts%scale_D)
                          IF (ABS(level(ip)) .GT. nb) THEN
                              Psi_global= Psi_global + opts%param_nb*&
                                  (ABS(level(ip))-nb)**2
                          ENDIF
                      ENDDO
                      level => Set_wps(Particles,Particles%level_id,read_only=.TRUE.)
                      wp    => Set_wps(Particles,Particles%adapt_wpid,read_only=.TRUE.)
                  ENDIF

                  IF (Psi_global .LT. Psi_global_old) THEN
                      IF (Psi_global .LT. Psi_1) THEN
                          Psi_1 = Psi_global
                          alpha1 = step
                      ENDIF
                      step = 2._MK * step
                      IF (step .GT. step_max) EXIT linesearch_loop
                  ELSE IF (Psi_global .GT. Psi_global_old) THEN
                      Psi_2 = Psi_global
                      alpha2 = step
                      step = 0.5_MK * step
                      IF (step .LT. step_min) EXIT linesearch_loop
                  ELSE
                      step = 10._MK * step
                      IF (step .GT. step_max) EXIT linesearch_loop
                  ENDIF

                  !FIXME
                  IF (ABS(Psi_global_old - Psi_global)/Psi_global_old .LT. 1E-4) &
                      EXIT linesearch_loop

              ENDDO linesearch_loop

              1000 CONTINUE

              IF (alpha1 .LT. 0._MK .OR. alpha2 .LT. 0._MK) THEN
                  step = MAX(step_min,MIN(step,step_max))
              ELSE
                  !Quadratic fit
                  Psi_1 = Psi_1 - Psi_global_old
                  Psi_2 = Psi_2 - Psi_global_old
                  step = 0.5_MK * (alpha1**2 * Psi_2 - alpha2**2 * Psi_1) / &
                      (alpha1*Psi_2 - alpha2*Psi_1)
              ENDIF

              !!---------------------------------------------------------------------!
              !! Choose best step size (from quadratic fit)
              !! Move particles (no need to move the ghosts, since we will have
              !! to get them through a local mapping anyway...)
              !!---------------------------------------------------------------------!
              !Move particles (including ghosts)
              xp => Get_xp(Particles,with_ghosts=.TRUE.)
              DO ip=1,Particles%Mpart
                  xp(1:ppm_dim,ip) = xp(1:ppm_dim,ip) + &
                      (step-step_previous) * Gradient_Psi(1:ppm_dim,ip)
              ENDDO
              xp => Set_xp(Particles)

#if debug_verbosity > 2
#ifdef __MPI
              CALL MPI_Allreduce(step*MAXVAL(ABS(Gradient_Psi(1:ppm_dim, &
                  1:Particles%Npart))),tmpvar2,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
              IF (ppm_rank.EQ.0) THEN
                  WRITE(cbuf,*) 'Moved particles: max displacement= ',tmpvar2
                  CALL ppm_write(ppm_rank,caller,cbuf,info)
              ENDIF
#endif
#endif

              CALL particles_updated_positions(Particles,info)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_updated_positions failed.',info)
                  info = -1
                  GOTO 9999
              ENDIF

              !!---------------------------------------------------------------------!
              !! /end Line search **
              !!---------------------------------------------------------------------!

#if writeout_verbosity > 2
#if debug_verbosity > 2
              !!---------------------------------------------------------------------!
              !! Some writeout
              !!---------------------------------------------------------------------!
              !xp => Get_xp(Particles)
              !D => Get_wps(Particles,Particles%D_id)
              !rcp => Get_wps(Particles,Particles%rcp_id)
              !IF (opts%write_pdb) THEN
                  !CALL sop_io(it_adapt+1000,debugdir,xp,ppm_dim,&
                      !Particles%Npart,Particles%Mpart,info,rcp=rcp,D=D)
                  !IF (info .NE. 0) THEN
                      !CALL ppm_write(ppm_rank,caller,'sop_io failed.',info)
                      !info = -1
                      !GOTO 9999
                  !ENDIF
              !ENDIF
              !IF (opts%write_xyz) THEN
                  !CALL sop_io_xyz(it_adapt+1000,debugdir,xp,ppm_dim,Particles%Npart,&
                      !Particles%Mpart,info,rcp=rcp,D=D)
                  !IF (info .NE. 0) THEN
                      !CALL ppm_write(ppm_rank,caller,'sop_io_xyz failed.',info)
                      !info = -1
                      !GOTO 9999
                  !ENDIF
              !ENDIF
              !xp => Set_xp(Particles,read_only=.TRUE.)
              !D => Set_wps(Particles,Particles%D_id,read_only=.TRUE.)
              !rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
#endif
#endif

#if debug_verbosity > 0
#ifdef __MPI
              CALL MPI_Allreduce(Psi_max,tmpvar2,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
              CALL MPI_Allreduce(Particles%Npart,tmpvari1,1,&
                  MPI_INTEGER,MPI_SUM,ppm_comm,info)
              CALL MPI_Allreduce(Particles%Mpart,tmpvari2,1,&
                  MPI_INTEGER,MPI_SUM,ppm_comm,info)
              IF (ppm_rank.EQ.0) THEN
                  WRITE(cbuf,'(A,I3,2(A,E11.4),A,I4,1X,I4,1X,A,I6,A,I6,A,E9.2)') &
                      'it_adapt= ',it_adapt,&
                      ' Psi_mean= ',Psi_global/REAL(tmpvari1,MK),' Psi_max= ',tmpvar2, &
                      ' Nneigh= ', Particles%nneighmin, Particles%nneighmax, &
                      'Np=',tmpvari1,' Mp=',tmpvari2,' step=',step
                  CALL ppm_write(ppm_rank,caller,cbuf,info)
              ENDIF
#endif
#endif

          ENDDO it_adapt_loop


          !------------------------------------------------------------------
          ! Deallocate D_tilde
          !------------------------------------------------------------------
          IF (PRESENT(wp_fun)) THEN
              CALL particles_allocate_wps(Particles,Particles%Dtilde_id,&
                  info,iopt=ppm_param_dealloc)
              IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'particles_allocate_wps (dealloc) failed',info)
                  info = -1
                  GOTO 9999
              ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Since particles have moved, we need to recompute the neigbour lists
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
          CALL particles_mapping_ghosts(Particles,topo_id,info)
          IF (info .NE. 0) THEN
              CALL ppm_write(ppm_rank,caller,&
                  'particles_apply_bc failed',info)
              info = -1
              GOTO 9999
          ENDIF
          CALL particles_neighlists(Particles,topo_id,info)
          IF (info .NE. 0) THEN
              CALL ppm_write(ppm_rank,caller,&
                  'particles_neighlists failed',info)
              info = -1
              GOTO 9999
          ENDIF

          !------------------------------------------------------------------
          ! In the (unlikely) event that neighmin is less than the critical
          ! number of neighbours needed, go back into the adaptation loop
          !------------------------------------------------------------------
          IF (Particles%nneighmin .LT. opts%nneigh_critical) THEN
              IF (it_adapt .LT. 9999) THEN
                  CALL ppm_write(ppm_rank,caller,&
                      'min nb of neighbr not enough. Going back to adaptation loop',info)
                  GOTO 7099
              ELSE
                  CALL ppm_write(ppm_rank,caller,&
                      'not enough neighbours and max nb of iterations exceeded.',&
                      info)
                  info = -1
                  GOTO 9999
              ENDIF
          ELSE
#if debug_verbosity > 0
              IF (ppm_rank.EQ.0) &
                  CALL ppm_write(ppm_rank,caller,&
                  'enough neighbours, we are good to go.',info)
#endif
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

      END SUBROUTINE DTYPE(sop_gradient_descent_ls)
