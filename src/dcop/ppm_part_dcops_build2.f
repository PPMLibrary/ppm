#if   __DIM == 2
SUBROUTINE ppm_dcop_compute2d(Particles,eta_id,info,interp,c,min_sv)
#elif __DIM == 3
SUBROUTINE ppm_dcop_compute3d(Particles,eta_id,info,interp,c,min_sv)
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
    !!! ratio h/epsilon
    REAL(MK),                  OPTIONAL, INTENT(  OUT)   :: min_sv
    !!! if present, compute the singular value decomposition of the 
    !!! vandermonde matrix for each operator and return the smallest one

    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    INTEGER                               :: i,j,k,ip,iq,beta(ppm_dim),ineigh
    INTEGER                               :: ncoeff,n_odd,np_target
    CHARACTER(LEN = 256)                  :: caller='ppm_dcop_compute'
    CHARACTER(LEN = 256)                  :: cbuf
    CHARACTER(LEN = 32)                   :: myformat
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: expo,byh
    REAL(MK),DIMENSION(:),ALLOCATABLE     :: byh0powerbeta
    REAL(MK)                              :: sv,dist2
    REAL(MK)                              :: nn_scaled
    REAL(MK), DIMENSION(:,:),POINTER      :: b=>NULL(),b_0=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: Z=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: Z_copy=>NULL()
    REAL(MK),DIMENSION(:)  ,ALLOCATABLE   :: d2_one2all
    REAL(MK),DIMENSION(:,:),ALLOCATABLE   :: dx
    INTEGER, DIMENSION(:,:),ALLOCATABLE   :: alpha,gamma
    INTEGER                               :: order_a
    INTEGER,DIMENSION(ppm_dim)            :: degree
    INTEGER,DIMENSION(:),ALLOCATABLE      :: sum_degree
    REAL(MK),DIMENSION(:),POINTER         :: coeffs=>NULL()
    INTEGER,DIMENSION(3)                  :: ldc
    INTEGER                               :: degree_poly
    LOGICAL                               :: cartesian,isinterp,adaptive,vector
    LOGICAL                               :: with_ghosts

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
    REAL(MK)                              :: min_sv_p
    
    ! haeckic: anisotropic additional variables
    LOGICAL                               :: anisotropic
    REAL(MK), DIMENSION(:,:),  POINTER    :: inv=>NULL()
    REAL(MK), DIMENSION(ppm_dim)          :: dx_temp

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
    
    ! haeckic
    IF (Particles%anisotropic) THEN
        anisotropic = .TRUE.
    ELSE
        anisotropic = .FALSE.
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
#ifdef __MKL
    !only used to compute the SVD, mainly for debugging.
    IF(PRESENT(min_sv)) THEN
        ALLOCATE(Z_copy(ncoeff,ncoeff))
    ENDIF
#endif

    vector =  Particles%ops%desc(eta_id)%vector
    !if true, then each term of the differential opearator is stored as one
    !component in eta. This is used when computing e.g. the gradient opearator.
    !if false, the same input parameters would yield an operator approximating
    ! the divergence operator.
    with_ghosts =  Particles%ops%desc(eta_id)%with_ghosts
    !if true, then the operator should be computed for ghost particles too. 
    !Note that the resulting values will be wrong for the ghost particles
    !that have some neighbours outside the ghost layers. Some of these particles
    !may also not have enough neighbours for the Vandermonde matrix to be
    !invertible. These particles will be skipped without raising a warning.

    IF (vector) THEN
        ldc(1) = nneighmax*nterms; 
    ELSE
        ldc(1) = nneighmax; 
    ENDIF
    IF (with_ghosts) THEN
        np_target = Particles%Mpart
    ELSE
        np_target = Particles%Npart
    ENDIF
    ldc(2) = np_target
    CALL ppm_alloc(Particles%ops%ker(eta_id)%vec,ldc,ppm_param_alloc_grow,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF
    eta => Get_dcop(Particles,eta_id,with_ghosts=with_ghosts)
    FORALL(i=1:ldc(1),j=1:np_target) eta(i,j)=0._MK

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
    
    ! haeckic
    IF (.NOT. adaptive .AND. .NOT. anisotropic) THEN
        byh = 1._MK/Particles%h_avg
    ENDIF

    IF (isinterp) THEN
        xp1 => Get_xp(Particles,with_ghosts=with_ghosts)
        xp2 => Get_xp(Particles_cross,with_ghosts=.TRUE.)
        nvlist => Particles%nvlist_cross
        vlist => Particles%vlist_cross
    ELSE
        xp1 => Get_xp(Particles,with_ghosts=.TRUE.)
        xp2 => xp1
        nvlist => Particles%nvlist
        vlist => Particles%vlist
    ENDIF
    coeffs => Particles%ops%desc(eta_id)%coeffs(1:nterms)

   ! haeckic
    IF (adaptive) THEN
        rcp => Get_wps(Particles,Particles%rcp_id,with_ghosts=with_ghosts)
    ENDIF
    IF (anisotropic) THEN
        inv => Get_wpv(Particles,Particles%G_id,with_ghosts=with_ghosts)
    ENDIF

    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        nn_sq => Get_wps(Particles_cross,Particles_cross%nn_sq_id,&
            with_ghosts=.TRUE.)
        !!! nearest-neighbour distances within Particles_cross 
        !!! (they must have been already computed)
    ENDIF

    particle_loop: DO ip = 1,np_target ! loop over all target particles

        IF (ip .GT. Particles%Npart .AND. nvlist(ip).LT.ncoeff) THEN
            !not enough neigbours for this ghost particle - skip it
            CYCLE particle_loop
        ENDIF

        Z = 0._MK
        b = b_0
        ! loop over their neighbors
        IF (adaptive) THEN
            ! haeckic changed this to 1.0
            byh = 1._MK/rcp(ip)
        ENDIF

        neighbour_loop: DO ineigh = 1,nvlist(ip) 
            iq = vlist(ineigh,ip) ! index in the "old particles" set

            ! haeckic: CHECK THE 2.0_mk 
            ! distance squared between the new particle and the old ones
            IF (anisotropic) THEN
               ! there is a 2 over rcp(ip)??
               dx_temp(1:ppm_dim) = (xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))
               
               ! Transformation
               dx(1,ineigh) = SUM(inv(1:ppm_dim,ip)*dx_temp(1:ppm_dim))
               dx(2,ineigh) = SUM(inv(ppm_dim+1:2*ppm_dim,ip)*dx_temp(1:ppm_dim))
               IF (ppm_dim .EQ. 3) THEN
                  dx(3,ineigh) = SUM(inv(2*ppm_dim+1:3*ppm_dim,ip)*dx_temp(1:ppm_dim))
               ENDIF

            ELSE
               dx(1:ppm_dim,ineigh) = (xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))*byh
            ENDIF

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

                    ! haeckic
                    IF (anisotropic) THEN
                       dx_temp(1:ppm_dim) = (xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))
               
                        ! Transformation
                        dist2 = SUM(inv(1:ppm_dim,ip)*dx_temp(1:ppm_dim))**2
                        dist2 = dist2 + SUM(inv(ppm_dim+1:2*ppm_dim,ip)*dx_temp(1:ppm_dim))**2
                        IF (ppm_dim .EQ. 3) THEN
                           dist2 = dist2 + SUM(inv(2*ppm_dim+1:3*ppm_dim,ip)*dx_temp(1:ppm_dim))**2
                        ENDIF
                        
                    ELSE
                    
                        dist2 = SUM((xp1(1:ppm_dim,ip) - xp2(1:ppm_dim,iq))**2)
                    
                    ENDIF
                    
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

#ifdef __MKL
        IF (PRESENT(min_sv)) THEN
            Z_copy = Z
            CALL ppm_matrix_svd(Z_copy,ncoeff,ncoeff,info,min_sv_p)
            IF (info .NE. 0) THEN
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
                DO i=1,ncoeff
                    WRITE(9002,*) (Z(i,j), j=1,ncoeff)
                ENDDO
                CALL ppm_write(ppm_rank,caller,&
                    'Z matrix written in file fort.9002',info)
                WRITE(9003,myformat) xp1(1:ppm_dim,ip)*byh
                DO ineigh = 1,nvlist(ip)
                    WRITE(9003,myformat) xp2(1:ppm_dim,vlist(ineigh,ip))*byh
                ENDDO
                CALL ppm_write(ppm_rank,caller,&
                    'h-scaled stencil written in file fort.9003',info)
                info = ppm_error_error
                CALL ppm_error(ppm_err_test_fail,caller,&
                    &  'ppm_matrix_svd failed',__LINE__,info)
                GOTO 9999
            ENDIF
            min_sv = MIN(min_sv,min_sv_p)
        ENDIF
#endif

        CALL solveLSE_n(Z,b,nterms,info)
        ! now b contain the solutions to the LSEs A*x_i=b_i for i=1:nterms
        IF (info .NE. 0) THEN
            IF (ip .GT. Particles%Npart) THEN
                !ignore error in matrix inversion for ghost particles
                ! simply skip it
                CYCLE particle_loop
            ENDIF
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
            CALL ppm_error(ppm_err_sub_failed,caller,&
                'Failed to solve the LSE', __LINE__,info)
            GOTO 9999
        ENDIF

         ! haeckic
        IF (anisotropic) THEN
            byh0powerbeta = 1.0_mk
        ELSE
            DO i=1,nterms 
                  byh0powerbeta(i) = byh**(sum_degree(i))
            ENDDO        
        ENDIF


        !------------------------------------------------------------------!
        ! Compute the operators
        !------------------------------------------------------------------!
        ! loop over old particles
        IF (isinterp) THEN
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c_value**2*d2_one2all(ineigh))
                DO j=1,ncoeff
                    IF (vector) THEN
                        eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) = &
                            eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) + &
                            b(j,1:nterms)*coeffs*byh0powerbeta* &
                            dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                            dx(3,ineigh)**gamma(3,j)* &
#endif
                            dx(2,ineigh)**gamma(2,j)* expo
                    ELSE
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
                    ENDIF
                ENDDO
            ENDDO
        ELSE
            DO ineigh = 1,nvlist(ip) 
                expo = exp(-c_value**2*d2_one2all(ineigh))
                DO j=1,ncoeff
                    IF (vector) THEN
                        eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) = &
                            eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) + &
                            b(j,1:nterms)*coeffs*byh0powerbeta* &
                        dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)
                    ELSE

                    eta(ineigh,ip) = eta(ineigh,ip) +&
                        SUM(b(j,1:nterms)*coeffs*byh0powerbeta)* &
                        dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                        dx(3,ineigh)**gamma(3,j)* &
#endif
                        dx(2,ineigh)**gamma(2,j)
                    ENDIF
                ENDDO
                IF (vector) THEN
                    eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) = &
                        eta(1+(ineigh-1)*nterms:ineigh*nterms,ip) * expo
                ELSE
                    eta(ineigh,ip) = eta(ineigh,ip) * expo
                ENDIF
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
    ! haeckic
    IF (anisotropic) THEN
        inv => Set_wpv(Particles,Particles%G_id,read_only=.TRUE.)
    ENDIF
    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        nn_sq => Set_wps(Particles_cross,Particles_cross%nn_sq_id,&
            read_only=.TRUE.)
    ENDIF

    !!---------------------------------------------------------------------!
    !! Finalize
    !!---------------------------------------------------------------------!
    DEALLOCATE(Z,b,b_0,d2_one2all,dx,gamma,alpha,sum_degree,byh0powerbeta)
#ifdef __MKL
    IF(PRESENT(min_sv)) THEN
        DEALLOCATE(Z_copy)
    ENDIF
#endif

    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

#if   __DIM == 2
END SUBROUTINE ppm_dcop_compute2d
#elif __DIM == 3
END SUBROUTINE ppm_dcop_compute3d
#endif

#undef __DIM
