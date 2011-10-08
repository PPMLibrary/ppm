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
    INTEGER,    DIMENSION(ppm_dim)   , INTENT(IN)       :: order
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
    coeffs = 1._MK; degree = 0
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
    
    CALL ppm_alloc(inv_transpose,(/ Particles%tensor_length /),ppm_param_dealloc,info)
    CALL particles_dcop_free(Particles,eta_id,info)

    wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
    eta => Set_dcop(Particles,eta_id)
    wp_grad => Set_wpv(Particles,wp_grad_id)
    inv => set_wpv(Particles,Particles%G_id,read_only=.TRUE.)

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

   SUBROUTINE get_hess_aniso(Particles,wp_id,wp_hess_id,order,info,with_ghosts)
   
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
    INTEGER,                            INTENT(IN)      :: wp_hess_id
    !!! resulting gradients for all particles 
    !!! resulting hessian matrix for all particles
    !!! 2D: [h11, h12, h21, h22]
    !!! 3D: [h11, h12, h13, h21, h22, h23, h31, h32, h33]
    INTEGER,    DIMENSION((ppm_dim-1)*3)   , INTENT(IN) :: order
    !!! order of approximation in dcops
    !!! 1 - xx, 2 - yy, 3 - xy
    !!! 1 - xx, 2 - yy-, 3 - zz, 4 - xy, 5 - xz, 6 - yz
    LOGICAL, OPTIONAL                                   :: with_ghosts
    !!! do we need the gradients also for the ghosts
    INTEGER,                            INTENT(OUT)     :: info
    !!! Return status, on success 0.

    !-------------------------------------------------------------------------
    ! local variables
    !-------------------------------------------------------------------------
    CHARACTER(LEN = ppm_char)                  :: caller = 'get_hess_aniso'
    REAL(MK),     DIMENSION(ppm_dim)           :: wp_grad_temp, wp_grad_temp2, wp_grad_temp3
    REAL(MK),     DIMENSION(((ppm_dim-1)*3))    :: coeffs
    INTEGER, DIMENSION(ppm_dim*((ppm_dim-1)*3)):: degree
    INTEGER                                    :: eta_id, i, ineigh, iq, ip
    REAL(MK)                                   :: t0
    REAL(MK), DIMENSION(:,:),POINTER           :: inv => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: wp => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: wp_hess => NULL()
    REAL(MK), DIMENSION(:,:),POINTER           :: eta => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: hess_matrix_iso => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: hess_matrix_aniso => NULL()
    REAL(MK), DIMENSION(:),  POINTER           :: inv_transpose => NULL()

    info = 0 ! change if error occurs
 
    CALL substart(caller,t0,info)
    
    ! check if wphessid exists, is ok!?
    IF (wp_hess_id.EQ.0 .OR. wp_id .EQ. 0) THEN
       info = ppm_error_error
       CALL ppm_error(ppm_err_alloc,caller,&
            'tried to calculate hessian, but wp_id or wp_hess_id not allocated!',__LINE__,info)
       GOTO 9999
    ENDIF

    ! 1. calculate the reference hessian matrix
    ! we assume continous second derivatives, therefore we can calc only one

    ! Compute hessian using PSE kernels
    coeffs = 1._MK; degree = 0

    ! For both 2D and 3D
    FORALL(i=1:ppm_dim) degree((i-1)*ppm_dim+i)=2 !Second derivative
    
    IF (ppm_dim .EQ. 3) THEN
      !xy
      degree(10:11) = 1
      
      !xz
      degree(13) = 1
      degree(15) = 1
      
      !yz
      degree(17:18) = 1

    ELSE
      !xy
      degree(5:6) = 1

    ENDIF
    eta_id = 0

    CALL particles_dcop_define(Particles,eta_id,coeffs,degree,&
       order,((ppm_dim-1)*3),info,name="hess",vector=.TRUE.)

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

    ! 2. Transform hessian for anisotorpic spaces
    ! hess_aniso = inv^T * hess_iso * inv
   
    wp_hess => Get_wpv(Particles,wp_hess_id)
    eta => Get_dcop(Particles,eta_id)
    wp => Get_wps(Particles,wp_id,with_ghosts=.TRUE.)
    inv => get_wpv(Particles,Particles%G_id,.TRUE.)

    DO ip = 1,Particles%Npart 
         wp_hess(1:ppm_dim**2,ip) = 0._MK
    ENDDO
    
    ! get the hessian grad_iso
    DO ip=1,Particles%Npart
         DO ineigh=1,Particles%nvlist(ip)
            iq=Particles%vlist(ineigh,ip)
            wp_hess(1:((ppm_dim-1)*3),ip) = wp_hess(1:((ppm_dim-1)*3),ip) + &
              & (wp(iq)-wp(ip)) * eta(1+(ineigh-1)*((ppm_dim-1)*3):ineigh*((ppm_dim-1)*3),ip)
         ENDDO
    ENDDO

   CALL ppm_alloc(hess_matrix_aniso,(/ Particles%tensor_length /),ppm_param_alloc_fit,info)
   CALL ppm_alloc(inv_transpose,(/ Particles%tensor_length /),ppm_param_alloc_fit,info)

    ! transform into grad_aniso
    DO ip=1,Particles%Npart
         IF (ppm_dim .EQ. 2) THEN
            inv_transpose = (/ inv(1,ip), inv(3,ip), inv(2,ip), inv(4,ip) /)

            !1st matrix multiplication
            wp_grad_temp =  (/ wp_hess(1,ip) , wp_hess(3,ip) /)
            wp_grad_temp2 = (/ wp_hess(3,ip) , wp_hess(2,ip) /)
            wp_hess(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_hess(2,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
            wp_hess(3,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp2)
            wp_hess(4,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp2)
            
            !2nd matrix multiplication
            wp_grad_temp =  (/ wp_hess(1,ip) , wp_hess(3,ip) /)
            wp_grad_temp2 = (/ wp_hess(2,ip) , wp_hess(4,ip) /)
            wp_hess(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_hess(2,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp2)
            wp_hess(3,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
            wp_hess(4,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp2)

         ELSE
            !3d
            inv_transpose = (/ inv(1,ip), inv(4,ip), inv(7,ip), &
                        &      inv(2,ip), inv(5,ip), inv(8,ip), &
                        &      inv(3,ip), inv(6,ip), inv(9,ip) /)

            !1st matrix multiplication
            wp_grad_temp =  (/ wp_hess(1,ip) , wp_hess(4,ip) , wp_hess(5,ip)/)
            wp_grad_temp2 = (/ wp_hess(4,ip) , wp_hess(2,ip), wp_hess(6,ip)/)
            wp_grad_temp3 = (/ wp_hess(5,ip) , wp_hess(6,ip), wp_hess(3,ip)/)

            wp_hess(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_hess(2,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
            wp_hess(3,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp)

            wp_hess(4,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp2)
            wp_hess(5,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp2)
            wp_hess(6,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp2)

            wp_hess(7,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp2)
            wp_hess(8,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp2)
            wp_hess(9,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp2)
            
            !2nd matrix multiplication
            wp_grad_temp =  (/ wp_hess(1,ip) , wp_hess(4,ip), wp_hess(7,ip)/)
            wp_grad_temp2 = (/ wp_hess(2,ip) , wp_hess(5,ip), wp_hess(8,ip)/)
            wp_grad_temp3 = (/ wp_hess(3,ip) , wp_hess(6,ip), wp_hess(9,ip)/)
            
            wp_hess(1,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp)
            wp_hess(2,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp2)
            wp_hess(3,ip) = SUM(inv_transpose(1:ppm_dim)*wp_grad_temp3)

            wp_hess(4,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp)
            wp_hess(5,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp2)
            wp_hess(6,ip) = SUM(inv_transpose(ppm_dim+1:2*ppm_dim)*wp_grad_temp3)

            wp_hess(7,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp)
            wp_hess(8,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp2)
            wp_hess(9,ip) = SUM(inv_transpose(2*ppm_dim+1:3*ppm_dim)*wp_grad_temp3)

         ENDIF
    ENDDO
    
    CALL ppm_alloc(hess_matrix_aniso,(/ Particles%tensor_length /),ppm_param_dealloc,info)
    CALL ppm_alloc(inv_transpose,(/ Particles%tensor_length /),ppm_param_dealloc,info)
    CALL particles_dcop_free(Particles,eta_id,info)

    wp => Set_wps(Particles,wp_id,read_only=.TRUE.)
    eta => Set_dcop(Particles,eta_id)
    wp_hess => Set_wpv(Particles,wp_hess_id)
    inv => set_wpv(Particles,Particles%G_id,read_only=.TRUE.)

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

END SUBROUTINE get_hess_aniso

