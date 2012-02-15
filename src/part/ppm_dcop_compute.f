#if   __DIM == 2
SUBROUTINE DTYPE(ppm_dcop_compute2d)(Pc,op_id,info,c,min_sv)
#elif __DIM == 3
SUBROUTINE DTYPE(ppm_dcop_compute3d)(Pc,op_id,info,c,min_sv)
#endif
    !!! Computes generalized DC operators
    !!! if the optional argument interp is true, the routine uses
    !!! one set of particles (Pc2) as input data and
    !!! compute the differential opearator on the new set
    !!! If it has to do be interpolating (if the quantity to be computed
    !!! is the field itself, ie. if the degree of one term of the 
    !!! differential operator is zero), then the nearest-neighbour 
    !!! distances within Pc2 must have been already computed.
    USE ppm_module_error
    IMPLICIT NONE

    DEFINE_MK()
    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    CLASS(DTYPE(ppm_t_particles))                        :: Pc
    !!! particles
    INTEGER,                             INTENT(IN   )   :: op_id
    !!! id of the operator kernel
    INTEGER,                             INTENT(  OUT)   :: info
    !!! non-zero on output if some error occurred
    !---------------------------------------------------------
    ! Optional arguments
    !---------------------------------------------------------
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
    INTEGER                               :: nterms
    !!! number of terms of the differential operator
    REAL(MK)                              :: c_value
    REAL(MK)                              :: min_sv_p

    TYPE(DTYPE(ppm_t_sop)),POINTER  :: Pc2 => NULL()
    TYPE(DTYPE(ppm_t_operator)),POINTER   :: op  => NULL()
    TYPE(DTYPE(ppm_t_neighlist)),POINTER  :: Nlist => NULL()

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

    op => Pc%ops%vec(op_id)%t

    isinterp = op%flags(ppm_ops_interp)
    vector = op%flags(ppm_ops_vector)
    !if true, then each term of the differential opearator is stored as one
    !component in eta. This is used when computing e.g. the gradient opearator.
    !if false, the same input parameters would yield an operator approximating
    ! the divergence operator.
    with_ghosts = op%flags(ppm_ops_inc_ghosts)
    !if true, then the operator should be computed for ghost particles too. 
    !Note that the resulting values will be wrong for the ghost particles
    !that have some neighbours outside the ghost layers. Some of these particles
    !may also not have enough neighbours for the Vandermonde matrix to be
    !invertible. These particles will be skipped without raising a warning.


    Nlist => Pc%neighs%vec(op%neigh_id)%t
    IF (.NOT. Nlist%uptodate) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Need to specify which set of particles &
            &   (Pc2) should be used for interpolation',&
            & __LINE__,info)
        GOTO 9999
    ENDIF

    nterms = op%desc%nterms
    ALLOCATE(sum_degree(nterms),byh0powerbeta(nterms))

    cartesian=Pc%flags(ppm_part_cartesian)
    IF (cartesian .AND. nterms .GT. 1) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Case where nterms>1 is not yet implemented for Cartesian particles',&
            __LINE__,info)
        GOTO 9999
    ENDIF

    !IF (Pc%adaptive) THEN
        !adaptive = .TRUE.
    !ELSE
        !adaptive = .FALSE.
    !ENDIF

    SELECT TYPE(Pc)
    TYPE IS (DTYPE(ppm_t_sop))
        CALL Pc%get(rcp,Pc%rcp_id,with_ghosts=with_ghosts)
        adaptive = .TRUE.
        IF (isinterp) THEN
            Pc2 => Pc%set_aPc%vec(op%P_id)%t
        ENDIF

        IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
            CALL Pc2%get(nn_sq,Pc2%nn_sq_id,with_ghosts=.TRUE.)
            !!! nearest-neighbour distances within Pc2 
            !!! (they must have been already computed)
        ENDIF
    CLASS DEFAULT
        byh = 1._MK/Pc%h_avg
        adaptive = .FALSE.
    END SELECT


    IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
        IF (Pc2%nn_sq_id .EQ. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,caller,&
                & 'need to call particles_nearest_neighbors first',__LINE__,info)
            GOTO 9999
        ENDIF
    ENDIF


    !Determine number of coefficients needed for this operator
    DO i=1,nterms
        degree = op%desc%degree(1+(i-1)*ppm_dim:i*ppm_dim)
        sum_degree(i) = SUM(degree)
        order_a=op%desc%order(i)

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

    IF (Nlist%nneighmin .LT. ncoeff) THEN
        CALL ppm_write(ppm_rank,caller,'Not enough neighbours',info)
        WRITE(cbuf,*) 'For this DC-operator, we need ',&
            ncoeff,' neighbours. We have ',Nlist%nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF


    ALLOCATE(d2_one2all(Nlist%nneighmax),dx(ppm_dim,Nlist%nneighmax),&
        Z(ncoeff,ncoeff),b(ncoeff,nterms),b_0(ncoeff,nterms),STAT=info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF

    !only used to compute the SVD, mainly for debugging.
    IF(PRESENT(min_sv)) THEN
        min_sv = HUGE(1._MK)
        ALLOCATE(Z_copy(ncoeff,ncoeff))
    ENDIF

    IF (vector) THEN
        ldc(1) = Nlist%nneighmax*nterms; 
    ELSE
        ldc(1) = Nlist%nneighmax; 
    ENDIF
    IF (with_ghosts) THEN
        np_target = Pc%Mpart
    ELSE
        np_target = Pc%Npart
    ENDIF
    ldc(2) = np_target
    CALL ppm_alloc(op%ker,ldc,ppm_param_alloc_grow,info)
    IF (info .NE. 0) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_alloc,caller,'Allocation failed',__LINE__,info)
        GOTO 9999
    ENDIF
    eta => Pc%get_dcop(op_id,with_ghosts=with_ghosts)
    FORALL(i=1:ldc(1),j=1:np_target) eta(i,j)=0._MK

    !----------------------------------------------------------------------!
    ! Compute diff op weights
    !----------------------------------------------------------------------!
    b_0 = 0._MK

    DO i=1,nterms
        degree = op%desc%degree(1+(i-1)*ppm_dim:i*ppm_dim)
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
    
    IF (isinterp) THEN
        CALL Pc%get_xp(xp1,with_ghosts=with_ghosts)
        CALL Pc2%get_xp(xp2,with_ghosts=.TRUE.)
    ELSE
        CALL Pc%get_xp(xp1,with_ghosts=.TRUE.)
        xp2 => xp1
    ENDIF
    nvlist => Nlist%nvlist
    vlist => Nlist%vlist
    coeffs => op%desc%coeffs(1:nterms)

    particle_loop: DO ip = 1,np_target ! loop over all target particles

        IF (ip .GT. Pc%Npart .AND. nvlist(ip).LT.ncoeff) THEN
            !not enough neigbours for this ghost particle - skip it
            CYCLE particle_loop
        ENDIF

        Z = 0._MK
        b = b_0
        ! loop over their neighbors
        IF (adaptive) THEN
            !byh = 2._MK/rcp(ip)
            byh = 0.5_MK/rcp(ip)
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
                        DTYPE(primitive)(SQRT(dist2/nn_sq(iq)) / 0.9_MK)

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

        CALL DTYPE(solveLSE_n)(Z,b,nterms,info)
        ! now b contain the solutions to the LSEs A*x_i=b_i for i=1:nterms
        IF (info .NE. 0) THEN
            IF (ip .GT. Pc%Npart) THEN
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
    CALL Pc%set_xp(xp1,read_only=.TRUE.)
    IF (isinterp) THEN
        CALL Pc2%set_xp(xp2,read_only=.TRUE.)
    ELSE
        xp2 => NULL()
    ENDIF
    eta => Pc%set_dcop(op_id)

    SELECT TYPE (Pc)
    TYPE IS (DTYPE(ppm_t_sop))
        CALL Pc%set(rcp,Pc%rcp_id,read_only=.TRUE.)
        IF (isinterp .AND. MINVAL(sum_degree).EQ.0) THEN
            CALL Pc2%set(nn_sq,Pc2%nn_sq_id,read_only=.TRUE.)
        ENDIF
    END SELECT

    !!---------------------------------------------------------------------!
    !! Finalize
    !!---------------------------------------------------------------------!
    DEALLOCATE(Z,b,b_0,d2_one2all,dx,gamma,alpha,sum_degree,byh0powerbeta)

    IF(PRESENT(min_sv)) THEN
        DEALLOCATE(Z_copy)
    ENDIF

    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error

#if   __DIM == 2
END SUBROUTINE DTYPE(ppm_dcop_compute2d)
#elif __DIM == 3
END SUBROUTINE DTYPE(ppm_dcop_compute3d)
#endif

#undef __DIM