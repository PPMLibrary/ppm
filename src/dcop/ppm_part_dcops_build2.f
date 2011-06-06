#if   __DIM == 2
SUBROUTINE ppm_dcop_compute2d(Particles,eta_id,info,interp,c)
#elif __DIM == 3
SUBROUTINE ppm_dcop_compute3d(Particles,eta_id,info,interp,c)
#endif
    !!! Computes generalized DC operators
    !!! if the optional argument interp is true, the routine uses
    !!! one set of particles (Particles_cross) as input data and
    !!! compute the differential opearator on the new set
    !!! If it has to do be interpolating (if the quantity to be computed
    !!! is the field itself, ie. if the degree of one term of the 
    !!! differential operator is zero), then the nearest-neighbour 
    !!! distances within Particles_cross must have been already computed.
    USE ppm_module_error
    USE ppm_module_particles_typedef
    USE ppm_module_particles
    IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION 
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    TYPE(ppm_t_particles),   POINTER,    INTENT(INOUT)   :: Particles
    !!! particles
    INTEGER,                             INTENT(IN   )   :: eta_id
    !!! id of the operator kernel
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
    LOGICAL,                   OPTIONAL, INTENT(IN   )   :: interp
    !!! true if the operator is to be computed using interpolating methods
    !!! (with data stored in another set of particles Particles%Particles_cross)
    REAL(MK),                  OPTIONAL, INTENT(IN   )   :: c
    !!! true if the operator is to be computed using interpolating methods
    !!! (with data stored in another set of particles Particles%Particles_cross)

    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    INTEGER                               :: i,j,k,ip,iq,beta(ppm_dim),ineigh
    INTEGER                               :: ncoeff,n_odd
    CHARACTER(LEN = 256)                  :: caller='ppm_part_dcops'
    CHARACTER(LEN = 256)                  :: cbuf
    CHARACTER(LEN = 32)                   :: myformat
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: expo,byh
    REAL(MK),DIMENSION(:),ALLOCATABLE     :: byh0powerbeta
    REAL(MK)                              :: sv,dist2
    REAL(MK)                              :: nn_scaled
    REAL(MK), DIMENSION(:,:),POINTER      :: b=>NULL(),b_0=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: Z=>NULL()
    REAL(MK),DIMENSION(:)  ,ALLOCATABLE   :: d2_one2all
    REAL(MK),DIMENSION(:,:),ALLOCATABLE   :: dx
    INTEGER, DIMENSION(:,:),ALLOCATABLE   :: alpha,gamma
    INTEGER                               :: order_a
    INTEGER,DIMENSION(ppm_dim)            :: degree
    INTEGER,DIMENSION(:),ALLOCATABLE      :: sum_degree
    REAL(MK),DIMENSION(:),POINTER         :: coeffs=>NULL()
    INTEGER,DIMENSION(3)                  :: ldc
    INTEGER                               :: degree_poly
    LOGICAL                               :: cartesian,isinterp,adaptive

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
    INTEGER                               :: nterms
    !!! number of terms of the differential operator
    TYPE(ppm_t_particles),POINTER         :: Particles_cross=>NULL()
    REAL(MK)                              :: c_value

    !!---------------------------------------------------------------------!
    !! Initialize
    !!---------------------------------------------------------------------!
    info = 0 ! change if error occurs
    CALL substart(caller,t0,info)

    IF (PRESENT(c)) THEN
        c_value=c
    ELSE
        c_value=1._MK
    ENDIF
    IF (.NOT.PRESENT(interp)) THEN
        isinterp=.FALSE.
    ELSE
        isinterp=interp
    ENDIF
    IF (isinterp) THEN
        IF (.NOT. ASSOCIATED(Particles%Particles_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Need to specify which set of particles &
                &   (particles_cross) should be used for interpolation',&
                & __LINE__,info)
            GOTO 9999
        ENDIF
        IF (.NOT. (Particles%neighlists_cross)) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Please compute xset neighbor lists first',&
                __LINE__,info)
            GOTO 9999
        ENDIF
        Particles_cross => Particles%Particles_cross
    ENDIF

    nterms = Particles%ops%desc(eta_id)%nterms
    ALLOCATE(sum_degree(nterms),byh0powerbeta(nterms))

    cartesian=Particles%cartesian
    IF (cartesian .AND. nterms .GT. 1) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Case where nterms>1 is not yet implemented for Cartesian particles',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    IF (Particles%adaptive) THEN
        adaptive = .TRUE.
    ELSE
        adaptive = .FALSE.
    ENDIF

    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        IF (Particles_cross%nn_sq_id .EQ. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                & 'need to call particles_nearest_neighbors first',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF


    !Determine number of coefficients needed for this operator
    DO i=1,nterms
        degree = Particles%ops%desc(eta_id)%degree(1+(i-1)*ppm_dim:i*ppm_dim)
        sum_degree(i) = SUM(degree)
        order_a=Particles%ops%desc(eta_id)%order(i)

        degree_poly = sum_degree(i) + order_a - 1

        IF (degree_poly .LT. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                'Negative degree for polynomial basis',__LINE__,info)
            GOTO 9999
        ENDIF

        IF (cartesian) THEN
            n_odd = 0
            DO j = 1,ppm_dim
                n_odd = n_odd + MOD(degree(j),2)
            ENDDO
            ncoeff = binomial((sum_degree(i)-n_odd)/2 + &
                CEILING(order_a/2.0) -1 + ppm_dim,ppm_dim) 
        ELSE 
            ncoeff = binomial(degree_poly+ppm_dim,ppm_dim)
        ENDIF

        !write(*,*) 'degree: ',degree, sum_degree(i)
        !write(*,*) 'degree_poly: ',degree_poly, ncoeff
    ENDDO

    IF (ncoeff.LE.0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Could not compute number of coefficients',__LINE__,info)
        GOTO 9999
    ENDIF

    ALLOCATE(gamma(ppm_dim,ncoeff),alpha(ppm_dim,ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF ! Generate the polynomial basis for particle approximation
    !   alpha is the approximation basis
    !   gamma is the template for DC operators
    ip = 0
    alpha = 0 

    DO i=0,degree_poly

#if   __DIM == 2
        loopj: DO j=0,i
            !if cartesian, it is not needed to compute coefficients 
            !for which gamma+order_d contains odd elements
            IF (cartesian) THEN 
                IF (MOD(j+degree(1),2)  .NE.0) CYCLE loopj
                IF (MOD(i-j+degree(2),2).NE.0) CYCLE loopj
            ENDIF
            ip=ip+1
            alpha(1,ip) = j
            alpha(2,ip) = i-j
        ENDDO loopj
#elif __DIM == 3
        DO j=0,i
            loopk: DO k=0,i-j
                !if cartesian, it is not needed to compute coefficients 
                !for which gamma+order_d contains odd elements
                IF (cartesian) THEN
                    IF (MOD(j+degree(1),2)    .NE.0) CYCLE loopk
                    IF (MOD(k+degree(2),2)    .NE.0) CYCLE loopk
                    IF (MOD(i-j-k+degree(3),2).NE.0) CYCLE loopk
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
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    ncoeff = ip
    gamma = alpha

    IF (isinterp) THEN
        nneighmin = Particles%nneighmin_cross
        nneighmax = Particles%nneighmax_cross
    ELSE
        nneighmin = Particles%nneighmin
        nneighmax = Particles%nneighmax
    ENDIF

    IF (nneighmin .LT. ncoeff) THEN
        CALL ppm_write(ppm_rank,caller,'Not enough neighbours',info)
        WRITE(cbuf,*) 'For this DC-operator, we need ',&
            ncoeff,' neighbours. We have ',nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(d2_one2all(nneighmax),dx(ppm_dim,nneighmax),Z(ncoeff,ncoeff),&
        b(ncoeff,nterms),b_0(ncoeff,nterms),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF

    ldc(1) = nneighmax; ldc(2) = Particles%Npart
    CALL ppm_alloc(Particles%ops%ker(eta_id)%vec,ldc,ppm_param_alloc_grow,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF
    !----------------------------------------------------------------------!
    ! Compute diff op weights
    !----------------------------------------------------------------------!
    b_0 = 0._MK

    DO i=1,nterms
        degree = Particles%ops%desc(eta_id)%degree(1+(i-1)*ppm_dim:i*ppm_dim)
        DO j=1,ncoeff
            IF (MAXVAL(ABS(alpha(:,j)-degree)).EQ.0) THEN
                b_0(j,i)=(-1)**(sum_degree(i))*factorial_m(degree,ppm_dim)
            ENDIF
        ENDDO

        !When applicable, and for stability reasons, set the zeroth moment to 5
        IF (.NOT.isinterp) THEN
            IF (SUM(alpha(1:ppm_dim,1)).EQ. 0 .AND. MOD(sum_degree(i),2) .EQ.0)&
                b_0(1,i) = 5._MK
        ENDIF
    ENDDO
    
    IF (.NOT. adaptive) THEN
        byh = 1._MK/Particles%h_avg
    ENDIF

    IF (isinterp) THEN
        xp1 => Get_xp(Particles)
        xp2 => Get_xp(Particles_cross,with_ghosts=.TRUE.)
        nvlist => Particles%nvlist_cross
        vlist => Particles%vlist_cross
    ELSE
        xp1 => Get_xp(Particles,with_ghosts=.TRUE.)
        xp2 => xp1
        nvlist => Particles%nvlist
        vlist => Particles%vlist
    ENDIF
    eta => Get_dcop(Particles,eta_id)
    FORALL(i=1:nneighmax,j=1:Particles%Npart) eta(i,j)=0._MK
    coeffs => Particles%ops%desc(eta_id)%coeffs(1:nterms)

    IF (adaptive) THEN
        rcp => Get_wps(Particles,Particles%rcp_id)
    ENDIF

    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        nn_sq => Get_wps(Particles_cross,Particles_cross%nn_sq_id)
        !!! nearest-neighbour distances within Particles_cross 
        !!! (they must have been already computed)
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

            expo = exp(-c_value**2*d2_one2all(ineigh))

            ! Fill matrix Z
            ! moments in alpha order: (alphas = rows)
            ! coefficients in gamma order (gammas = cols):

            DO j=1,ncoeff
                DO i=1,ncoeff
                    beta = alpha(1:ppm_dim,i) + gamma(1:ppm_dim,j)

                    Z(i,j) = Z(i,j) + &
                        dx(1,ineigh)**beta(1) * &
#if   __DIM == 3
                        dx(3,ineigh)**beta(3) * &
#endif
                        dx(2,ineigh)**beta(2) * expo
                ENDDO
            ENDDO

        ENDDO neighbour_loop

        IF (isinterp) THEN
            DO i=1,nterms
                IF (sum_degree(i).NE.0) CYCLE
                ! Assemble the rhs for the linear system that has to be solved for 
                ! interpolating functions
                DO ineigh = 1,nvlist(ip)
                    !rescaled nearest-neighbour distance  
                    iq = vlist(ineigh,ip)
                    dist2 = SUM((xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))**2)
                    eta(ineigh,ip) = &
                        primitive(SQRT(dist2/nn_sq(iq)) / 0.9_MK)

                    !reuse the variable to assemble the rhs
                    DO j=1,ncoeff
                        b(j,i) = b(j,i) - &
                            eta(ineigh,ip)* &
#if   __DIM == 2
                            dx(1,ineigh)**alpha(1,j) * dx(2,ineigh)**alpha(2,j)
#elif __DIM == 3
                            dx(1,ineigh)**alpha(1,j) * dx(2,ineigh)**alpha(2,j) * &
                                dx(3,ineigh)**alpha(3,j)
#endif
                    ENDDO
                ENDDO
            ENDDO
        ENDIF


        CALL solveLSE_n(Z,b,nterms,info)
        ! now b contain the solutions to the LSEs A*x_i=b_i for i=1:nterms
        IF (info .NE. 0) THEN
            !writes the coordinate of the stencil that lead to the error
            IF (ppm_dim .EQ. 2 ) THEN
                myformat = TRIM(ADJUSTL('(2(E30.22))'))
            ELSE
                myformat = TRIM(ADJUSTL('(3(E30.22))'))
            ENDIF
            WRITE(9000,myformat) xp1(1:ppm_dim,ip)
            DO ineigh = 1,nvlist(ip)
                WRITE(9000,myformat) xp2(1:ppm_dim,vlist(ineigh,ip))
            ENDDO
            CALL ppm_write(ppm_rank,caller,&
                'stencil written in file fort.9000',info)

            WRITE(9003,myformat) xp1(1:ppm_dim,ip)*byh
            DO ineigh = 1,nvlist(ip)
                WRITE(9003,myformat) xp2(1:ppm_dim,vlist(ineigh,ip))*byh
            ENDDO
            CALL ppm_write(ppm_rank,caller,&
                'h-scaled stencil written in file fort.9003',info)


            info = ppm_error_error
            CALL ppm_error(999,caller,'Failed to solve the LSE', __LINE__,info)
            GOTO 9999
        ENDIF

        DO i=1,nterms 
            byh0powerbeta(i) = byh**(sum_degree(i))
        ENDDO

        !------------------------------------------------------------------!
        ! Compute the operators
        !------------------------------------------------------------------!
        ! loop over old particles
        IF (isinterp) THEN
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c_value**2*d2_one2all(ineigh))
                DO j=1,ncoeff

                    eta(ineigh,ip) = eta(ineigh,ip) + &
                        SUM(b(j,1:nterms)*coeffs*byh0powerbeta)* &
                        dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)* expo

                        !note: do not factorise out expo, like for eta. 
                        ! eta_interp already contains some data (the primitive &
                        ! function) that should not be multiplied by expo.
                ENDDO
            ENDDO
        ELSE
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c_value**2*d2_one2all(ineigh))
                DO j=1,ncoeff

                    eta(ineigh,ip) = eta(ineigh,ip) +&
                        SUM(b(j,1:nterms)*coeffs*byh0powerbeta)* &
                        dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)
                ENDDO
                eta(ineigh,ip) = eta(ineigh,ip) * expo
            ENDDO
        ENDIF

    ENDDO particle_loop

    coeffs=> NULL()
    xp1 => Set_xp(Particles,read_only=.TRUE.)
    IF (isinterp) THEN
        xp2 => Set_xp(Particles_cross,read_only=.TRUE.)
    ELSE
        xp2 => NULL()
    ENDIF
    eta => Set_dcop(Particles,eta_id)
    IF (adaptive) THEN
        rcp => Set_wps(Particles,Particles%rcp_id,read_only=.TRUE.)
    ENDIF
    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        nn_sq => Set_wps(Particles_cross,Particles_cross%nn_sq_id,&
            read_only=.TRUE.)
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
    DEALLOCATE(sum_degree,STAT=info)
    DEALLOCATE(byh0powerbeta,STAT=info)

    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

#if   __DIM == 2
END SUBROUTINE ppm_dcop_compute2d
#elif __DIM == 3
END SUBROUTINE ppm_dcop_compute3d
#endif

#undef __DIM
