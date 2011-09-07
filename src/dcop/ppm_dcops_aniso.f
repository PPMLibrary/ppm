   SUBROUTINE get_grad_aniso(Particles,wp_id,wp_grad_id,order,info,with_ghosts)

    USE ppm_module_map
    USE ppm_module_particles_typedef
    USE ppm_module_write
    USE ppm_module_particles
    IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(INOUT)   :: Particles
    !!! Data structure containing the particles
    INTEGER,                            INTENT(IN)      :: wp_id
    !!! source of wp
    INTEGER,                            INTENT(IN)      :: wp_grad_id
    !!! resulting gradients for all particles 
    INTEGER,    DIMENSION(ppm_dim)   , INTENT(IN)        :: order
    !!! order of approximation in dcops
    LOGICAL, OPTIONAL                                   :: with_ghosts
    !!! do we need the gradients also for the ghosts
    INTEGER,                            INTENT(OUT)     :: info
    !!! Return status, on success 0.

    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: caller = 'get_grad_aniso'
    REAL(MK),     DIMENSION(ppm_dim)           :: coeffs,wp_grad_temp
    INTEGER,      DIMENSION(ppm_dim*ppm_dim)   :: degree
    INTEGER                                    :: eta_id, i, ineigh, iq, ip
    REAL(MK)                                   :: t0
    REAL(MK), DIMENSION(:,:),POINTER           :: inv => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: wp => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: wp_grad => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: eta => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: inv_transpose => NULL()

    info = 0 ! change if error occurs

    CALL substart(caller,t0,info)

    ! check if wpgradid exists, is ok!?
    IF (wp_grad_id.EQ.0 .OR. wp_id .EQ. 0) THEN
       info = ppm_error_error
       CALL ppm_error(ppm_err_alloc,caller,&
            'tried to calculate gradient, but wp_id or wp_grad_id not allocated!',__LINE__,info)
       GOTO 9999
    ENDIF

    ! 1. calculate the reference gradient

    !Compute gradients using PSE kernels
    coeffs=1._MK; degree = 0
    FORALL(i=1:ppm_dim) degree((i-1)*ppm_dim+i)=1 !Gradient
    eta_id = 0

    CALL particles_dcop_define(Particles,eta_id,coeffs,degree,&
       order,ppm_dim,info,name="gradient",vector=.TRUE.)

    IF (info.NE.0) THEN
       info = ppm_error_error
       CALL ppm_error(ppm_err_sub_failed,caller,&
             'particles_dcop_define failed', __LINE__,info)
       GOTO 9999
    ENDIF

    CALL particles_dcop_compute(Particles,eta_id,info)

    IF (info.NE.0) THEN
       info = ppm_error_error
       CALL ppm_error(ppm_err_sub_failed,caller,&
             'particles_dcop_compute failed',__LINE__,info)
       GOTO 9999
    ENDIF

    ! 2. Transform gradient for anisotorpic spaces
    ! grad_aniso = inv^T * grad_iso
   
    wp_grad => Get_wpv(Particles,wp_grad_id)
    eta => Get_dcop(Particles,eta_id)
    wp => Get_wps(Particles,wp_id,with_ghosts=.TRUE.)
    inv => get_wpv(Particles,Particles%G_id,.TRUE.)

    DO ip = 1,Particles%Npart 
         wp_grad(1:ppm_dim,ip) = 0._MK
    ENDDO
    
    ! get the gradient grad_iso
    DO ip=1,Particles%Npart
         DO ineigh=1,Particles%nvlist(ip)
            iq=Particles%vlist(ineigh,ip)
            wp_grad(1:ppm_dim,ip) = wp_grad(1:ppm_dim,ip) + &
              & (wp(iq)-wp(ip)) * eta(1+(ineigh-1)*ppm_dim:ineigh*ppm_dim,ip)
         ENDDO
    ENDDO

    CALL ppm_alloc(inv_transpose,(/ Particles%tensor_length /),ppm_param_alloc_fit,info)

    ! transform into grad_aniso
    DO ip=1,Particles%Npart
         IF (ppm_dim .EQ. 2) THEN
            inv_transpose = (/ inv(1,ip), inv(3,ip), inv(2,ip), inv(4,ip) /)
            wp_grad_temp = wp_grad(1:ppm_dim,ip)
            wp_grad(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_grad(2,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
         ELSE
            inv_transpose = (/ inv(1,ip), inv(4,ip), inv(7,ip), &
                        &      inv(2,ip), inv(5,ip), inv(8,ip), &
                        &      inv(3,ip), inv(6,ip), inv(9,ip) /)
            wp_grad_temp = wp_grad(1:ppm_dim,ip)
            wp_grad(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_grad(2,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
            wp_grad(3,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp)
         ENDIF
    ENDDO

    wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
    eta => Set_dcop(Particles,eta_id)
    wp_grad => Set_wpv(Particles,wp_grad_id)

    IF (PRESENT(with_ghosts)) THEN
        IF (with_ghosts) THEN
            !---------------------------------------------------------------------!
            ! Update ghosts
            !---------------------------------------------------------------------!
            CALL particles_mapping_ghosts(Particles,Particles%active_topoid,info)
            IF (info .NE. 0) THEN
                  CALL ppm_write(ppm_rank,caller,&
                     'particles_mapping_ghosts failed.',info)
                  info = -1
                  GOTO 9999
            ENDIF
         ENDIF
    ENDIF


    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error

   END SUBROUTINE get_grad_aniso


! define a get_hess for anisotropic particles
   SUBROUTINE get_hess_aniso(Particles, hess, info)
   
    USE ppm_module_map
    USE ppm_module_particles_typedef
    USE ppm_module_write
    USE ppm_module_particles

    IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !-------------------------------------------------------------------------
    !  Arguments
    !-------------------------------------------------------------------------
    TYPE(ppm_t_particles), POINTER,     INTENT(IN)    :: Particles
    !!! Data structure containing the particles
    REAL(MK),DIMENSION(:,:), POINTER                  :: hess
    !!! resulting hessian matrix for all particles
    !!! 2D: [h11, h12, h21, h22]
    !!! 3D: [h11, h12, h13, h21, h22, h23, h31, h32, h33]
    INTEGER,                             INTENT(OUT)    :: info
    !!! Return status, on success 0.

    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: caller = 'get_hess_aniso'
    REAL(MK),     DIMENSION(ppm_dim)           :: coeffs,wp_grad_temp
    INTEGER,      DIMENSION(ppm_dim)           :: order
    INTEGER,      DIMENSION(ppm_dim*ppm_dim)   :: degree
    INTEGER                                    :: eta_id, i, ineigh, iq, ip
    REAL(MK)                                   :: t0
    REAL(MK), DIMENSION(:,:),POINTER           :: inv => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: wp => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: wp_grad => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: eta => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: inv_transpose => NULL()

    info = 0 ! change if error occurs
 
    CALL substart(caller,t0,info)
    
    ! CODE HERE

    CALL substop(caller,t0,info)

    9999  CONTINUE ! jump here upon error

END SUBROUTINE get_hess_aniso

