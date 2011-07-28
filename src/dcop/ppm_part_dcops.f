#ifdef __2D
SUBROUTINE ppm_part_dcops_2d(Particles,eta_id,c,info,&
        order_deriv,order_approx,islaplacian,iscartesian,&
        isinterp,Particles_old,nn_sq_id)
#endif

#ifdef __3D
SUBROUTINE ppm_part_dcops_3d(Particles,eta_id,c,info,&
        order_deriv,order_approx,islaplacian,iscartesian,&
        isinterp,Particles_old,nn_sq_id)
#endif

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank
    USE ppm_module_error
    USE ppm_module_particles
    IMPLICIT NONE

    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    TYPE(ppm_t_particles),   POINTER,    INTENT(INOUT)   :: Particles
    !!! particles
    INTEGER,                             INTENT(INOUT)   :: eta_id
    !!! id of the operator kernel
    REAL(MK),                            INTENT(IN   )   :: c
    !!! ratio h/epsilon
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
    INTEGER,                   OPTIONAL, INTENT(IN   )   :: order_approx
    !!! order of approximation
    INTEGER,DIMENSION(ppm_dim),OPTIONAL, INTENT(IN   )   :: order_deriv
    !!! degree of the derivative that we want to approximate
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: islaplacian
    !!! special treatment for the laplacian operator. Default is FALSE
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: iscartesian
    !!! special treatment if particles are on a grid. Default is FALSE
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: isinterp
    !!! special treatment if the operator is interpolating. Default is FALSE
    TYPE(ppm_t_particles),OPTIONAL,POINTER,INTENT(IN   ) :: Particles_old
    !!! Second set of particles. If present, then the kernel uses data
    !!! on the old set (used to compute interpolation kernels e.g.)
    INTEGER,                   OPTIONAL, INTENT(IN   )   :: nn_sq_id
    !!! id where the nearest-neighbour distances within Particles_old 
    !!! are stored (they must have been already computed)

    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    INTEGER                               :: i,j,k,ip,iq,beta(ppm_dim),ineigh
    INTEGER                               :: ncoeff,n_odd
    CHARACTER(LEN = 256)                  :: caller='ppm_part_dcops'
    CHARACTER(LEN = 256)                  :: cbuf
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: expo,byh0powerbeta,byh
    REAL(MK)                              :: sv,dist2
    REAL(MK)                              :: nn_scaled
    REAL(MK), DIMENSION(:),  POINTER      :: b=>NULL(),b_0=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: Z=>NULL()
    REAL(MK),DIMENSION(:)  ,ALLOCATABLE   :: d2_one2all
    REAL(MK),DIMENSION(:,:),ALLOCATABLE   :: dx
    INTEGER, DIMENSION(:,:),ALLOCATABLE   :: alpha,gamma
    INTEGER                               :: order_a
    INTEGER,DIMENSION(ppm_dim)            :: order_d
    INTEGER                               :: degree_poly,sum_order_d
    LOGICAL                               :: check_stability=.FALSE.
    LOGICAL                               :: laplacian,cartesian,interp,adaptive

    REAL(MK), DIMENSION(:,:),  POINTER    :: xp1=>NULL()
    !!! particles (or points) where the operators are computed
    REAL(MK), DIMENSION(:,:),  POINTER    :: xp2=>NULL()
    !!! particles that contain the data to be used (xp1 can be equal to xp2) 
    !!! A case where xp1 .ne. xp2 is for interpolation.
    REAL(MK), DIMENSION(:),  POINTER      :: rcp=>NULL()
    REAL(MK), DIMENSION(:),  POINTER      :: nn_sq=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: eta=>NULL()
    INTEGER, DIMENSION(:),   POINTER      :: nvlist=>NULL()
    INTEGER, DIMENSION(:,:), POINTER      :: vlist=>NULL()
    INTEGER                               :: nneighmin,nneighmax

    !!---------------------------------------------------------------------!
    !! Initialize
    !!---------------------------------------------------------------------!
    info = 0
    order_d=-1

    IF (.NOT.PRESENT(islaplacian)) THEN
        laplacian=.FALSE.
    ELSE
        laplacian=islaplacian
    ENDIF
    IF (.NOT.PRESENT(iscartesian)) THEN
        cartesian=.FALSE.
    ELSE
        cartesian=iscartesian
    ENDIF
    IF (.NOT.PRESENT(isinterp)) THEN
        interp=.FALSE.
    ELSE
        interp=isinterp
    ENDIF
    IF (.NOT.PRESENT(order_deriv)) THEN
        !default is to compute the laplacian operator
        laplacian=.TRUE.
    ELSE
        order_d=order_deriv
    ENDIF
    IF (Particles%adaptive) THEN
        adaptive = .TRUE.
    ELSE
        adaptive = .FALSE.
    ENDIF

    IF (interp .AND. .NOT. PRESENT(Particles_old)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Not enough arguments - need Particles_old for interpolation',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (interp .AND. .NOT. PRESENT(nn_sq_id)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Not enough arguments - need nn_sq_id for interpolation',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    IF (laplacian) THEN
        order_d = 2
        sum_order_d = 2
    ELSE
        sum_order_d = SUM(order_d)
    ENDIF

    IF (MINVAL(order_d) .LT. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Negative degree for derivatives',__LINE__,info)
        GOTO 9999
    ENDIF
    IF (.NOT.PRESENT(order_approx)) THEN
        !default is to compute a second-order approximation
        order_a=2
    ELSE
        order_a=order_approx
    ENDIF

    degree_poly = sum_order_d + order_a - 1

    IF (degree_poly .LT. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Negative degree for polynomial basis',__LINE__,info)
        GOTO 9999
    ENDIF

    IF (cartesian) THEN
        n_odd = 0
        DO i = 1,ppm_dim
            n_odd = n_odd + MOD(order_d(i),2)
        ENDDO
        ncoeff = binomial((sum_order_d-n_odd)/2 + &
            CEILING(order_a/2.0) -1 + ppm_dim,ppm_dim) 
    ELSE 
        ncoeff = binomial(degree_poly+ppm_dim,ppm_dim)
    ENDIF

    IF (ncoeff.LE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Could not compute number of coefficients',__LINE__,info)
        GOTO 9999
    ENDIF

    ALLOCATE(gamma(ppm_dim,ncoeff),alpha(ppm_dim,ncoeff))
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,&
            'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF ! Generate the polynomial basis for particle approximation
    !   alpha is the approximation basis
    !   gamma is the template for DC operators
    ip = 0
    alpha = 0 

    DO i=0,degree_poly

#ifdef __2D
        loopj: DO j=0,i
            !if cartesian, it is not needed to compute coefficients 
            !for which gamma+order_d contains odd elements
            IF (cartesian) THEN 
                IF (MOD(j+order_d(1),2)  .NE.0) CYCLE loopj
                IF (MOD(i-j+order_d(2),2).NE.0) CYCLE loopj
            ENDIF
            ip=ip+1
            alpha(1,ip) = j
            alpha(2,ip) = i-j
        ENDDO loopj
#else
        DO j=0,i
            loopk: DO k=0,i-j
                !if cartesian, it is not needed to compute coefficients 
                !for which gamma+order_d contains odd elements
                IF (cartesian) THEN
                    IF (MOD(j+order_d(1),2)    .NE.0) CYCLE loopk
                    IF (MOD(k+order_d(2),2)    .NE.0) CYCLE loopk
                    IF (MOD(i-j-k+order_d(3),2).NE.0) CYCLE loopk
                ENDIF

                ip=ip+1
                alpha(1,ip) = j
                alpha(2,ip) = k
                alpha(3,ip) = i-j-k
            ENDDO loopk
        ENDDO
#endif
    ENDDO

    IF (ip.NE.ncoeff) THEN
        WRITE(cbuf,*) 'Something wrong when computing coefficients. Theory: ',&
            ncoeff,', we have ',ip
        CALL pwrite(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    ncoeff = ip
    gamma = alpha

    IF (interp) THEN
        nneighmin = Particles%nneighmin_cross
        nneighmax = Particles%nneighmax_cross
    ELSE
        nneighmin = Particles%nneighmin
        nneighmax = Particles%nneighmax
    ENDIF

    IF (nneighmin .LT. ncoeff) THEN
        CALL pwrite(ppm_rank,caller,'Not enough neighbours',info)
        WRITE(cbuf,*) 'For this DC-operator, we need ',&
            ncoeff,' neighbours. We have ',nneighmin
        CALL pwrite(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(d2_one2all(nneighmax))
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(dx(ppm_dim,nneighmax))
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(Z(ncoeff,ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(b(ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(b_0(ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    CALL particles_allocate_wpv(Particles,eta_id,nneighmax,info,&
        zero=.TRUE.,iopt=ppm_param_alloc_grow)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'particles_allocate_wpv failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    !----------------------------------------------------------------------!
    ! Compute diff op weights
    !----------------------------------------------------------------------!
    b_0 = 0._MK

    IF (.NOT.laplacian) THEN
        DO j=1,ncoeff
            IF (MAXVAL(ABS(alpha(:,j)-order_d)).EQ.0) THEN
                b_0(j)=(-1)**(sum_order_d)*factorial_m(order_d,ppm_dim)
            ENDIF
        ENDDO
    ELSE
        !special treatment for the Laplacian
        DO j=1,ncoeff
            IF (SUM(alpha(:,j)).EQ.2 .AND. MAXVAL(alpha(:,j)).EQ.2) THEN
                b_0(j)=2
            ENDIF
        ENDDO
    ENDIF

    !When applicable, and for stability reasons, set the zeroth moment to 5
    IF (.NOT.interp) THEN
        IF (SUM(alpha(1:ppm_dim,1)).EQ. 0 .AND. MOD(sum_order_d,2) .EQ.0) b_0(1) = 5._MK
    ENDIF
    
    IF (.NOT. adaptive) THEN
        byh = 1._MK/Particles%h_avg
    ENDIF

    IF (interp) THEN
        xp1 => Get_xp(Particles)
        xp2 => Get_xp(Particles_old,with_ghosts=.TRUE.)
        nvlist => Particles%nvlist_cross
        vlist => Particles%vlist_cross
    ELSE
        xp1 => Get_xp(Particles,with_ghosts=.TRUE.)
        xp2 => xp1
        nvlist => Particles%nvlist
        vlist => Particles%vlist
    ENDIF
    eta => Get_wpv(Particles,eta_id)
    IF (adaptive) THEN
        rcp => Get_wps(Particles,Particles%rcp_id)
    ENDIF

    IF (interp) THEN
        nn_sq => Get_wps(Particles_old,nn_sq_id)
    ENDIF

    particle_loop: DO ip = 1,Particles%Npart ! loop over all (new) real particles

        Z = 0._MK
        b = b_0
        ! loop over their neighbors
        IF (adaptive) THEN
            byh = 2._MK/rcp(ip)
        ENDIF

        neighbour_loop: DO ineigh = 1,nvlist(ip) 
            iq = vlist(ineigh,ip) ! index in the "old particles" set

            ! distance squared between the new particle and the old ones
            dx(1:ppm_dim,ineigh) = (xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))*byh
            d2_one2all(ineigh) = SUM(dx(1:ppm_dim,ineigh)**2)

            expo = exp(-c**2*d2_one2all(ineigh))

            ! Fill matrix Z
            ! moments in alpha order: (alphas = rows)
            ! coeffs in gamma order (gammas = cols):

            DO j=1,ncoeff
                DO i=1,ncoeff
                    beta = alpha(1:ppm_dim,i) + gamma(1:ppm_dim,j)

                    Z(i,j) = Z(i,j) + &
                        dx(1,ineigh)**beta(1) * &
#ifdef __3D
                        dx(3,ineigh)**beta(3) * &
#endif
                        dx(2,ineigh)**beta(2) * expo
                ENDDO
            ENDDO

        ENDDO neighbour_loop

        IF (interp .and. sum_order_d.eq.0) THEN
            ! Assemble the rhs for the linear system that has to be solved for 
            ! interpolating functions
            DO ineigh = 1,nvlist(ip)
                !rescaled nearest-neighbour distance  
                iq = vlist(ineigh,ip)
                dist2 = SUM((xp1(1:ndim,ip) - xp2(1:ndim,iq))**2)
                eta(ineigh,ip) = &
                    primitive(SQRT(dist2/nn_sq(iq)) / 0.9_MK)

                !reuse the variable to assemble the rhs
                DO j=1,ncoeff
                    b(j) = b(j) - &
                        eta(ineigh,ip)* &
#ifdef __2D
                        dx(1,ineigh)**alpha(1,j) * dx(2,ineigh)**alpha(2,j)
#else
                        dx(1,ineigh)**alpha(1,j) * dx(2,ineigh)**alpha(2,j) * &
                            dx(3,ineigh)**alpha(3,j)
#endif
                ENDDO
            ENDDO
        ENDIF


        CALL solveLSE(Z,b,info)
        ! now b contain the solution to the LSEs A*x=b 

        IF (info .NE. 0) THEN
            CALL pwrite(ppm_rank,caller,'solveLSE failed',info)
            info = -1
            GOTO 9999
        ENDIF

        byh0powerbeta = byh**(sum_order_d)

        !------------------------------------------------------------------!
        ! Compute the operators
        !------------------------------------------------------------------!
        ! loop over old particles
        IF (interp) THEN
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c**2*d2_one2all(ineigh))
                DO j=1,ncoeff

                    eta(ineigh,ip) = eta(ineigh,ip) + b(j)*  &
                        dx(1,ineigh)**gamma(1,j)* &
#ifdef __3D
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)* expo * byh0powerbeta

                        !note: do not factorise out expo, like for eta. 
                        ! eta_interp already contains some data (the primitive &
                        ! function) that should not be multiplied by expo.
                ENDDO
            ENDDO
        ELSE
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c**2*d2_one2all(ineigh))
                DO j=1,ncoeff

                    eta(ineigh,ip) = eta(ineigh,ip) + b(j)*  &
                        dx(1,ineigh)**gamma(1,j)* &
#ifdef __3D
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)
                ENDDO
                eta(ineigh,ip) = eta(ineigh,ip) * expo * byh0powerbeta
            ENDDO
        ENDIF

    ENDDO particle_loop

    xp1 => Set_xp(Particles,read_only=.TRUE.)
    IF (interp) THEN
        xp2 => Set_xp(Particles_old,read_only=.TRUE.)
    ELSE
        xp2 => NULL()
    ENDIF
    eta => Set_wpv(Particles,eta_id)
    IF (adaptive) THEN
        rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
    ENDIF
    IF (interp) THEN
        nn_sq => Set_wps(Particles_old,nn_sq_id,read_only=.TRUE.)
    ENDIF

    !!---------------------------------------------------------------------!
    !! Finalize
    !!---------------------------------------------------------------------!
    DEALLOCATE(Z,STAT=info)
    DEALLOCATE(b,STAT=info)
    DEALLOCATE(b_0,STAT=info)
    DEALLOCATE(d2_one2all,STAT=info)
    DEALLOCATE(dx,STAT=info)
    DEALLOCATE(gamma,STAT=info)
    DEALLOCATE(alpha,STAT=info)

    9999 CONTINUE ! jump here upon error

#ifdef __2D
END SUBROUTINE ppm_part_dcops_2d
#endif

#ifdef __3D
END SUBROUTINE ppm_part_dcops_3d
#endif
