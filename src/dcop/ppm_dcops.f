#if   __DIM  == 2
#if   __MODE == __UNIFORM
#if   __KIND == __SINGLE_PRECISION
SUBROUTINE ppm_dcops_unif_2d_s(xp,h,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE ppm_dcops_unif_2d_d(xp,h,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#endif
#elif __MODE == __ADAPTIVE
#if __KIND   == __SINGLE_PRECISION
SUBROUTINE ppm_dcops_adapt_2d_s(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE ppm_dcops_adapt_2d_d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#endif
#endif

#elif __DIM  == 3
#if   __MODE == __UNIFORM
#if   __KIND == __SINGLE_PRECISION
SUBROUTINE ppm_dcops_unif_3d_s(xp,h,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE ppm_dcops_unif_3d_d(xp,h,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#endif
#elif __MODE == __ADAPTIVE
#if __KIND   == __SINGLE_PRECISION
SUBROUTINE ppm_dcops_adapt_3d_s(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#elif __KIND == __DOUBLE_PRECISION
SUBROUTINE ppm_dcops_adapt_3d_d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#endif
#endif
#endif

    USE ppm_module_data, ONLY: ppm_dim,ppm_rank
    USE ppm_module_error
    IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    !---------------------------------------------------------
    ! arguments
    !---------------------------------------------------------
    REAL(MK), DIMENSION(:,:),POINTER,    INTENT(IN   )   :: xp
    !!! position of particles
#if   __MODE == __UNIFORM
    REAL(MK),                            INTENT(IN   )   :: h
    !!! characteristic inter-particle spacing
#elif __MODE == __ADAPTIVE
    REAL(MK), DIMENSION(:), POINTER,     INTENT(IN   )   :: rcp
    !!! local characteristic inter-particle spacing
#endif
    INTEGER,                             INTENT(IN   )   :: Npart
    !!! number of real particles
    INTEGER,                             INTENT(IN   )   :: Mpart
    !!! number of (real+ghosts) particles
    INTEGER,  DIMENSION(:),  POINTER,    INTENT(IN   )   :: nvlist
    !!! number of neighbours for each particle
    INTEGER,  DIMENSION(:,:),POINTER,    INTENT(IN   )   :: vlist
    !!! verlet lists
    REAL(MK), DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT)   :: eta
    !!! operator kernel
    REAL(MK),                            INTENT(IN   )   :: c
    !!! ratio h/epsilon
    INTEGER,                             INTENT(IN   )   :: nneighmin
    !!! minimum number of neighbours
    INTEGER,                             INTENT(IN   )   :: nneighmax
    !!! maximum number of neighbours
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

    !---------------------------------------------------------
    ! local variables
    !---------------------------------------------------------
    INTEGER                               :: i,j,k,ip,iq,beta(ppm_dim),ineigh
    INTEGER                               :: ncoeff
    CHARACTER(LEN = 256)                  :: caller='ppm_dcops'
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
    INTEGER                               :: order_a,n_odd
    INTEGER,DIMENSION(ppm_dim)            :: order_d
    INTEGER                               :: degree_poly,sum_order_d
    LOGICAL                               :: check_stability=.FALSE.
    LOGICAL                               :: laplacian,cartesian

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
    IF (.NOT.PRESENT(order_deriv)) THEN
        !default is to compute the laplacian operator
        laplacian=.TRUE.
    ELSE
        order_d=order_deriv
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
    alpha=0 !initialised to negative values - used for error checking only

    DO i=0,degree_poly

#if __DIM == 2
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
#elif __DIM == 3
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
    ncoeff = ip
    gamma = alpha

    IF (nneighmin .LT. ncoeff) THEN
        CALL ppm_write(ppm_rank,caller,'Not enough neighbours',info)
        WRITE(cbuf,*) 'For this DC-operator, we need ',&
            ncoeff,' neighbours. We have ',nneighmin
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(d2_one2all(nneighmax))
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(dx(ppm_dim,nneighmax))
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(Z(ncoeff,ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(b(ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    ALLOCATE(b_0(ncoeff),STAT=info)
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF

    IF (ALLOCATED(eta)) DEALLOCATE(eta)
    ALLOCATE(eta(nneighmax,Npart))
    IF (info .NE. 0) THEN
        CALL ppm_write(ppm_rank,caller,'Allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    DO ip=1,Npart
        DO ineigh=1,nneighmax
            eta(ineigh,ip) = 0._MK
        ENDDO
    ENDDO
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
    IF (SUM(alpha(1:ppm_dim,1)).EQ. 0 .AND. MOD(sum_order_d,2) .EQ.0) b_0(1) = 5._MK


#if __MODE == __UNIFORM
    byh = 1._MK/h
#endif
    particle_loop: DO ip = 1,Npart ! loop over all (new) real particles

        Z = 0._MK
        b = b_0
        ! loop over their neighbors
#if __MODE == __ADAPTIVE
    byh = 2._MK/rcp(ip)
#endif

        neighbour_loop: DO ineigh = 1,nvlist(ip) 
            iq = vlist(ineigh,ip) ! index in the "old particles" set

            ! distance squared between the new particle and the old ones
            dx(1:ppm_dim,ineigh) = (xp(1:ppm_dim,ip) - xp(1:ppm_dim,iq))*byh
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
#if __DIM == 3
                        dx(3,ineigh)**beta(3) * &
#endif
                        dx(2,ineigh)**beta(2) * expo
                ENDDO
            ENDDO

        ENDDO neighbour_loop

        CALL solveLSE(Z,b,info)
        ! now b contain the solution to the LSEs A*x=b 

        IF (info .NE. 0) THEN
                    CALL ppm_write(ppm_rank,caller,'solveLSE failed',info)
                    info = -1
                    GOTO 9999
        ENDIF

        byh0powerbeta = byh**(sum_order_d)

        !------------------------------------------------------------------!
        ! Compute the laplacian operator eta
        !------------------------------------------------------------------!
        ! loop over old particles
        neighbour_loop2: DO ineigh = 1,nvlist(ip) 
            expo = exp(-c**2*d2_one2all(ineigh))
            DO j=1,ncoeff

                eta(ineigh,ip) = eta(ineigh,ip) + b(j)*  &
                    dx(1,ineigh)**gamma(1,j)* &
#if __DIM == 3
                    dx(3,ineigh)**gamma(3,j)* &
#endif
                    dx(2,ineigh)**gamma(2,j)
            ENDDO

            eta(ineigh,ip) = eta(ineigh,ip)*&
                exp(-c**2 * d2_one2all(ineigh))*byh0powerbeta
        ENDDO neighbour_loop2

    ENDDO particle_loop

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

#if   __DIM  == 2
#if   __MODE == __UNIFORM
#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE ppm_dcops_unif_2d_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_dcops_unif_2d_d
#endif
#elif __MODE == __ADAPTIVE
#if __KIND   == __SINGLE_PRECISION
END SUBROUTINE ppm_dcops_adapt_2d_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_dcops_adapt_2d_d
#endif
#endif

#elif __DIM  == 3
#if   __MODE == __UNIFORM
#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE ppm_dcops_unif_3d_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_dcops_unif_3d_d
#endif
#elif __MODE == __ADAPTIVE
#if __KIND   == __SINGLE_PRECISION
END SUBROUTINE ppm_dcops_adapt_3d_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_dcops_adapt_3d_d
#endif
#endif
#endif

#undef __KIND
#undef __MODE
#undef __DIM
