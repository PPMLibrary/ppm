
!-------------------------------------------------------------
! Primitive function as defined in Chen et al., 
!             Int. J. Numer. Meth. Engng 2003; 56:935–960.
! (here, the quartic spline)
!-------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
FUNCTION primitive_s(x)

    IMPLICIT NONE
    INTEGER, PARAMETER :: MK = ppm_kind_single
    !arguments
    REAL(MK), INTENT(IN) :: x
    REAL(MK)             :: primitive_s
    IF(x.GT.1._MK) THEN 
        primitive_s = 0._MK 
    ELSE
        primitive_s = 1._MK + x**2 * (-6._MK + x * (8._MK -3._MK * x))
    ENDIF

#elif  __KIND == __DOUBLE_PRECISION
FUNCTION primitive_d(x)

    IMPLICIT NONE
    INTEGER, PARAMETER :: MK = ppm_kind_double
    !arguments
    REAL(MK), INTENT(IN) :: x
    REAL(MK)             :: primitive_d
    IF(x.GT.1._MK) THEN 
        primitive_d = 0._MK 
    ELSE
        primitive_d = 1._MK + x**2 * (-6._MK + x * (8._MK -3._MK * x))
    ENDIF

#endif
END FUNCTION


#if    __KIND == __SINGLE_PRECISION
SUBROUTINE solveLSEs(A,x_or_b,info)
#elif  __KIND == __DOUBLE_PRECISION
SUBROUTINE solveLSEd(A,x_or_b,info)
#endif

    !=======================================================================!
    ! solves the LSE A*x=b
    ! if necessary, singularities are removed &
    ! (WARNING: no check is performed on the 
    ! validity of the solution after that)
    ! x_or_b is expected to contain b on input
    ! x_or_b contains x on output
    ! A is not altered

#ifdef __MKL
    USE mkl95_lapack
    USE mkl95_blas
#endif
    IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION 
    INTEGER, PARAMETER :: MK = ppm_kind_single
#elif  __KIND == __DOUBLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    INTEGER , DIMENSION (:)  , POINTER   :: indx,valid,indxnew,roworder
    REAL(MK), DIMENSION (:)  , POINTER   :: bnew,bnew_2,exact_b,exact_b_2,real_b
    REAL(MK), DIMENSION (:,:), POINTER   :: bnew_n,exact_b_n
    REAL(MK), DIMENSION (:,:), POINTER   :: Anew,Acopy
    REAL(MK)                             :: tolerance_lse = 1e-1

    ! arguments
    REAL(MK), DIMENSION (:,:), POINTER,INTENT(IN   ):: A 
    REAL(MK), DIMENSION (:)  , POINTER  :: x_or_b
    INTEGER, INTENT(OUT)                  :: info
    ! local variables
    INTEGER                               :: n,nnew
    INTEGER                               :: i,j,inew,jnew,itemp
    LOGICAL                               :: check = .TRUE.
    CHARACTER(LEN = 256),PARAMETER    :: caller = 'solveLSE'
    CHARACTER(LEN = 256)              :: cbuf
    REAL(MK)                            :: closetozero

    !=======================================================================!
    ! init
    info = 0
    n = SIZE(A,1)
    !closetozero = SQRT(EPSILON(closetozero))
    closetozero = 100._MK*(EPSILON(closetozero))

    ALLOCATE(indx(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(valid(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(Acopy(n,n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF
    Acopy = A

    IF (check) THEN
        ALLOCATE(exact_b(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        exact_b = x_or_b

        ALLOCATE(real_b(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !=======================================================================!
    ! compute the LU factorization of A
#ifdef __MKL
    CALL getrf(Acopy,indx,info)
#else
    CALL dgetrf(n,n,Acopy,n,indx,info)
#endif
    ! if info = i > 0, then U_ii = 0 and the matrix is singular
    ! this case will be dealt with below when checking for singularities
    IF (info .LT. 0) THEN
        WRITE(cbuf,'(A,I2)') ' getrf failed with info = ', info
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    !=======================================================================!
    ! check for singularities

    ! mark redundant eqs/coefficients
    valid = 1
    DO j=1,n
        IF (ABS(Acopy(j,j)) .LT. closetozero) THEN
            valid(j) = 0
        ENDIF
    ENDDO

    nnew = SUM(valid) ! size of reduced system

    ! if necessary, remove redundant equations/coefficients
    IF (nnew .NE. n) THEN

        ALLOCATE(roworder(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        DO i= 1,n
            roworder(i) = i
        ENDDO
        DO i= 1,n
            itemp = roworder(i)
            roworder(i) = roworder(indx(i))
            roworder(indx(i)) = itemp
        ENDDO

        !====================================================================!
        ! create new LSE

        ! allocation
        ALLOCATE(bnew(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(indxnew(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(Anew(nnew,nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        ! delete singular rows/columns in A
        ! delete corresponding entries in b
        inew = 0
        DO i = 1,n

            IF (valid(i) .EQ. 1) THEN

                inew = inew + 1
                bnew(inew) = x_or_b(roworder(i))

                jnew = 0
                DO  j = 1,n
                    IF (valid(j) .EQ. 1) THEN

                        jnew = jnew + 1
                        Anew(inew,jnew) = A(roworder(i),roworder(j))

                    ENDIF
                ENDDO

            ENDIF

        ENDDO

        !====================================================================!
        ! solve new LSE
#ifdef __MKL
        CALL getrf(Anew,indxnew,info)
#else
        CALL dgetrf(nnew,nnew,Anew,nnew,indxnew,info)
#endif
        IF (info .NE. 0) THEN
            !get the value of info
            WRITE(cbuf,'(A,I2)') ' getrf new failed with info = ', info

            IF (info.GT.0) THEN
                OPEN(81,FILE='DumpA.dat',FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
                DO i=1,n
                    WRITE(81,'(10E20.12)') A(i,:)
                ENDDO
                CLOSE(81)
            ENDIF

            !print the error message
            CALL ppm_write(ppm_rank,caller,cbuf,info)

            info = -1
            GOTO 9999
        ENDIF

#ifdef __MKL
        CALL getrs(Anew,indxnew,bnew,'N',info)
#else
        CALL dgetrs('N',nnew,1,Anew,nnew,indxnew,bnew,nnew,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(cbuf,'(A,I2)') ' getrs new failed with info = ', info
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF

        !====================================================================!
        ! sort solution (coefficients) into old places, remaining places are
        ! set to zero 
        jnew = 0
        DO j = 1,n
            x_or_b(j) = 0._MK
        ENDDO
        DO j = 1,n

            IF (valid(j) .EQ. 1) THEN
                jnew = jnew + 1
                x_or_b(roworder(j)) = bnew(jnew)
            ENDIF

        ENDDO

        DEALLOCATE(roworder)
        DEALLOCATE(indxnew)
        DEALLOCATE(bnew)
        DEALLOCATE(Anew)

    ELSE ! no singularities

        !====================================================================!
        ! solve originial LSE
#ifdef __MKL
        CALL getrs(Acopy,indx,x_or_b,'N',info)
#else
        CALL dgetrs('N',n,1,Acopy,n,indx,x_or_b,n,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': getrs failed'
            info = -1
            GOTO 9999
        ENDIF

    ENDIF

    !=======================================================================!
    ! check if equations were solved right (if you trust in LAPACK and this 
    ! code PUT THIS IN IF BRANCH FOR THE REDUCED SYSTEM)
    IF (check) THEN
#ifdef __MKL
        CALL gemv(A,x_or_b,real_b,1._MK,0._MK)
#else
        CALL dgemv('N',SIZE(A,1),SIZE(A,2),1._MK,A,&
            SIZE(A,1),x_or_b,1,0._MK,real_b,1)
#endif
        IF (SUM(ABS(real_b - exact_b)) .GT. tolerance_LSE) THEN
            WRITE(*,*)'WARNING from ', TRIM(caller),': no solution!'
            CALL ppm_write(ppm_rank,caller,&
                'Error in moment conditions too large.',info)
            WRITE(cbuf,'(A,E12.3)')'The l1-norm of the error is ',&
                SUM(ABS(real_b - exact_b))
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            DO i=1,n
                DO j=1,n
                    WRITE(9001,'(1(E30.22))') A(i,j)
                ENDDO
            ENDDO
            CALL ppm_write(ppm_rank,caller,&
                'moment matrix written in file fort.9001',info)
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !=======================================================================!
    ! dealloc
    DEALLOCATE(valid)
    DEALLOCATE(indx)
    DEALLOCATE(Acopy)
    IF (check) THEN
        DEALLOCATE(exact_b)
        DEALLOCATE(real_b)
    ENDIF

    9999 CONTINUE ! jump here upon error

#if    __KIND == __SINGLE_PRECISION
END SUBROUTINE solveLSEs
#elif  __KIND == __DOUBLE_PRECISION
END SUBROUTINE solveLSEd
#endif

#if    __KIND == __SINGLE_PRECISION
SUBROUTINE solveLSE_2s(A,x_or_b,x_or_b_2,info)
#elif  __KIND == __DOUBLE_PRECISION
SUBROUTINE solveLSE_2d(A,x_or_b,x_or_b_2,info)
#endif

    !=======================================================================!
    ! solves the LSE A*x1=b1 and A*x2=b2
    ! if necessary, singularities are removed 
    ! (WARNING: no check is performed on the 
    ! validity of the solution after that)
    ! x_or_b is expected to contain b on input
    ! x_or_b contains x on output
    ! A is not altered

#ifdef __MKL
    USE mkl95_lapack
    USE mkl95_blas
#endif
    IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#else
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    ! arguments
    REAL(MK), DIMENSION (:,:), POINTER,INTENT(IN   )  :: A 
    REAL(MK), DIMENSION (:)  , POINTER  :: x_or_b, x_or_b_2
    INTEGER, INTENT(OUT)                  :: info
    ! local variables
    INTEGER , DIMENSION (:)  , POINTER   :: indx,valid,indxnew,roworder
    REAL(MK), DIMENSION (:)  , POINTER   :: bnew,bnew_2,exact_b,exact_b_2,real_b
    REAL(MK), DIMENSION (:,:), POINTER   :: bnew_n,exact_b_n
    REAL(MK), DIMENSION (:,:), POINTER   :: Anew,Acopy
    REAL(MK)                             :: tolerance_lse = 1e-1
    INTEGER                               :: n,nnew
    INTEGER                               :: i,j,inew,jnew,itemp
    LOGICAL                               :: check = .TRUE.
    CHARACTER(LEN = 256),PARAMETER    :: caller = 'solveLSE_2'
    CHARACTER(LEN = 256)              :: cbuf
    REAL(MK)                            :: closetozero

    !=======================================================================!
    ! init
    info = 0
    n = SIZE(A,1)
    !closetozero = SQRT(EPSILON(closetozero))
    closetozero = 100._MK*(EPSILON(closetozero))

    ALLOCATE(indx(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(valid(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(Acopy(n,n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF
    Acopy = A

    IF (check) THEN
        ALLOCATE(exact_b(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        exact_b = x_or_b

        ALLOCATE(exact_b_2(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        exact_b_2 = x_or_b_2

        ALLOCATE(real_b(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !=======================================================================!
    ! compute the LU factorization of A
#ifdef __MKL
    CALL getrf(Acopy,indx,info)
#else
    CALL dgetrf(n,n,Acopy,n,indx,info)
#endif
    ! if info = i > 0, then U_ii = 0 and the matrix is singular
    ! this case will be dealt with below when checking for singularities
    IF (info .LT. 0) THEN
        WRITE(cbuf,'(A,I2)') ' getrf failed with info = ', info
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    !=======================================================================!
    ! check for singularities

    ! mark redundant eqs/coefficients
    valid = 1
    DO j=1,n
        IF (ABS(Acopy(j,j)) .LT. closetozero) THEN
            valid(j) = 0
        ENDIF
    ENDDO

    nnew = SUM(valid) ! size of reduced system

    ! if necessary, remove redundant equations/coefficients
    IF (nnew .NE. n) THEN

        ALLOCATE(roworder(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        DO i= 1,n
            roworder(i) = i
        ENDDO
        DO i= 1,n
            itemp = roworder(i)
            roworder(i) = roworder(indx(i))
            roworder(indx(i)) = itemp
        ENDDO

        !====================================================================!
        ! create new LSE

        ! allocation
        ALLOCATE(bnew(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(bnew_2(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(indxnew(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(Anew(nnew,nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        ! delete singular rows/columns in A
        ! delete corresponding entries in b
        inew = 0
        DO i = 1,n

            IF (valid(i) .EQ. 1) THEN

                inew = inew + 1
                bnew(inew)   = x_or_b(roworder(i))
                bnew_2(inew) = x_or_b_2(roworder(i))

                jnew = 0
                DO  j = 1,n
                    IF (valid(j) .EQ. 1) THEN

                        jnew = jnew + 1
                        Anew(inew,jnew) = A(roworder(i),roworder(j))

                    ENDIF
                ENDDO

            ENDIF

        ENDDO

        !====================================================================!
        ! solve new LSE
#ifdef __MKL
        CALL getrf(Anew,indxnew,info)
#else
        CALL dgetrf(nnew,nnew,Anew,nnew,indxnew,info)
#endif
        IF (info .NE. 0) THEN
            !get the value of info
            WRITE(cbuf,'(A,I2)') ' getrf new failed with info = ', info

            IF (info.GT.0) THEN
                OPEN(81,FILE='DumpA.dat',FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
                DO i=1,n
                    WRITE(81,'(10E20.12)') A(i,:)
                ENDDO
                CLOSE(81)
            ENDIF

            !print the error message
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF

#ifdef __MKL
        CALL getrs(Anew,indxnew,bnew,'N',info)
#else
        CALL dgetrs('N',nnew,1,Anew,nnew,indxnew,bnew,nnew,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(cbuf,'(A,I2)') ' getrs new failed with info = ', info
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF
#ifdef __MKL
        CALL getrs(Anew,indxnew,bnew_2,'N',info)
#else
        CALL dgetrs('N',nnew,1,Anew,nnew,indxnew,bnew_2,nnew,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(cbuf,'(A,I2)') ' getrs new failed with info = ', info
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF

        !====================================================================!
        ! sort solution (coefficients) into old places, remaining places are
        ! set to zero 
        jnew = 0
        DO j = 1,n
            x_or_b(j) = 0._MK
            x_or_b_2(j) = 0._MK
        ENDDO
        DO j = 1,n

            IF (valid(j) .EQ. 1) THEN
                jnew = jnew + 1
                x_or_b(roworder(j)) = bnew(jnew)
                x_or_b_2(roworder(j)) = bnew_2(jnew)
            ENDIF

        ENDDO

        DEALLOCATE(roworder)
        DEALLOCATE(indxnew)
        DEALLOCATE(bnew)
        DEALLOCATE(bnew_2)
        DEALLOCATE(Anew)

    ELSE ! no singularities

        !====================================================================!
        ! solve originial LSE
#ifdef __MKL
        CALL getrs(Acopy,indx,x_or_b,'N',info)
#else
        CALL dgetrs('N',n,1,Acopy,n,indx,x_or_b,n,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': getrs failed'
            info = -1
            GOTO 9999
        ENDIF
#ifdef __MKL
        CALL getrs(Acopy,indx,x_or_b_2,'N',info)
#else
        CALL dgetrs('N',n,1,Acopy,n,indx,x_or_b_2,n,info)
#endif
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': getrs failed'
            info = -1
            GOTO 9999
        ENDIF

    ENDIF

    !=======================================================================!
    ! check if equations were solved right (if you trust in LAPACK and this 
    ! code PUT THIS IN IF BRANCH FOR THE REDUCED SYSTEM)
    IF (check) THEN
#ifdef __MKL
        CALL gemv(A,x_or_b,real_b,1._MK,0._MK)
#else
        CALL dgemv('N',SIZE(A,1),SIZE(A,2),1._MK,A,&
            SIZE(A,1),x_or_b,1,0._MK,real_b,1)
#endif
        IF (SUM(ABS(real_b - exact_b)) .GT. tolerance_LSE) THEN
            WRITE(*,*)'WARNING from ', TRIM(caller),': no solution!'
            CALL ppm_write(ppm_rank,caller,&
                'Error in moment conditions 1 too large.',info)
            WRITE(cbuf,'(A,E12.3)')'The l1-norm of the error is ',&
                SUM(ABS(real_b - exact_b))
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            DO i=1,n
                DO j=1,n
                    WRITE(9001,'(1(E30.22))') A(i,j)
                ENDDO
            ENDDO
            CALL ppm_write(ppm_rank,caller,&
                'moment matrix written in file fort.9001',info)
            info = -1
            GOTO 9999
        ENDIF
#ifdef __MKL
        CALL gemv(A,x_or_b_2,real_b,1._MK,0._MK)
#else
        CALL dgemv('N',SIZE(A,1),SIZE(A,2),1._MK,A,&
            SIZE(A,1),x_or_b_2,1,0._MK,real_b,1)
#endif
        IF (SUM(ABS(real_b - exact_b_2)) .GT. tolerance_LSE) THEN
            WRITE(*,*)'WARNING from ', TRIM(caller),': no solution!'
            CALL ppm_write(ppm_rank,caller,&
                'Error in moment conditions 2 too large.',info)
            WRITE(cbuf,'(A,E12.3)')'The l1-norm of the error is ',&
                SUM(ABS(real_b - exact_b_2))
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !=======================================================================!
    ! dealloc
    DEALLOCATE(valid)
    DEALLOCATE(indx)
    DEALLOCATE(Acopy)
    IF (check) THEN
        DEALLOCATE(exact_b)
        DEALLOCATE(exact_b_2)
        DEALLOCATE(real_b)
    ENDIF

    9999 CONTINUE ! jump here upon error

#if    __KIND == __SINGLE_PRECISION
END SUBROUTINE solveLSE_2s
#elif  __KIND == __DOUBLE_PRECISION
END SUBROUTINE solveLSE_2d
#endif

#if    __KIND == __SINGLE_PRECISION
SUBROUTINE solveLSE_ns(A,x_or_b,n_eq,info)
#elif  __KIND == __DOUBLE_PRECISION
SUBROUTINE solveLSE_nd(A,x_or_b,n_eq,info)
#endif

    !=======================================================================!
    ! solves the LSE A*x_i=b_i  for i=1..n_eq
    ! if necessary, singularities are removed 
    ! (WARNING: no check is performed on the 
    ! validity of the solution after that)
    ! x_or_b is expected to contain b on input
    ! x_or_b contains x on output
    ! A is not altered

#ifdef __MKL
    USE mkl95_lapack
    USE mkl95_blas
#endif
    IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#else
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
    ! arguments
    REAL(MK), DIMENSION (:,:), POINTER,INTENT(IN   )  :: A 
    REAL(MK), DIMENSION (:,:), POINTER,INTENT(INOUT)  :: x_or_b
    INTEGER, INTENT(IN)                                 :: n_eq
    INTEGER, INTENT(OUT)                                :: info

    ! local variables
    INTEGER , DIMENSION (:)  , POINTER   :: indx,valid,indxnew,roworder
    REAL(MK), DIMENSION (:)  , POINTER   :: bnew,bnew_2,exact_b,exact_b_2,real_b
    REAL(MK), DIMENSION (:,:), POINTER   :: bnew_n,exact_b_n
    REAL(MK), DIMENSION (:,:), POINTER   :: Anew,Acopy
    REAL(MK)                             :: tolerance_lse = 1e-1
    INTEGER                               :: n,nnew
    INTEGER                               :: i,j,k,inew,jnew,itemp
    LOGICAL                               :: check = .TRUE.
    CHARACTER(LEN = 256),PARAMETER        :: caller = 'solveLSE_n'
    CHARACTER(LEN = 256)                  :: cbuf
    REAL(MK)                            :: closetozero

    !=======================================================================!
    ! init
    info = 0
    n = SIZE(A,1)
    !closetozero = SQRT(EPSILON(closetozero))
    closetozero = 100._MK*(EPSILON(closetozero))

    ALLOCATE(indx(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(valid(n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF

    ALLOCATE(Acopy(n,n),STAT=info)
    IF (info .NE. 0) THEN
        WRITE(*,*)caller,': allocation failed.'
        info = -1
        GOTO 9999
    ENDIF
    Acopy = A

    IF (check) THEN
        ALLOCATE(exact_b_n(n,n_eq),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        exact_b_n = x_or_b

        ALLOCATE(real_b(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
    ENDIF

    !=======================================================================!
    ! compute the LU factorization of A
#ifdef __MKL
    CALL getrf(Acopy,indx,info)
#else
    CALL dgetrf(n,n,Acopy,n,indx,info)
#endif
    ! if info = i > 0, then U_ii = 0 and the matrix is singular
    ! this case will be dealt with below when checking for singularities
    IF (info .LT. 0) THEN
        WRITE(cbuf,'(A,I2)') ' getrf failed with info = ', info
        CALL ppm_write(ppm_rank,caller,cbuf,info)
        info = -1
        GOTO 9999
    ENDIF

    !=======================================================================!
    ! check for singularities

    ! mark redundant eqs/coefficients
    valid = 1
    DO j=1,n
        IF (ABS(Acopy(j,j)) .LT. closetozero) THEN
            valid(j) = 0
        ENDIF
    ENDDO

    nnew = SUM(valid) ! size of reduced system

    ! if necessary, remove redundant equations/coefficients
    IF (nnew .NE. n) THEN

        ALLOCATE(roworder(n),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        DO i= 1,n
            roworder(i) = i
        ENDDO
        DO i= 1,n
            itemp = roworder(i)
            roworder(i) = roworder(indx(i))
            roworder(indx(i)) = itemp
        ENDDO

        !====================================================================!
        ! create new LSE

        ! allocation
        ALLOCATE(bnew_n(nnew,n_eq),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(indxnew(nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF
        ALLOCATE(Anew(nnew,nnew),STAT=info)
        IF (info .NE. 0) THEN
            WRITE(*,*)caller,': allocation failed.'
            info = -1
            GOTO 9999
        ENDIF

        ! delete singular rows/columns in A
        ! delete corresponding entries in b
        inew = 0
        DO i = 1,n

            IF (valid(i) .EQ. 1) THEN

                inew = inew + 1
                bnew_n(inew,1:n_eq)   = x_or_b(roworder(i),1:n_eq)

                jnew = 0
                DO  j = 1,n
                    IF (valid(j) .EQ. 1) THEN

                        jnew = jnew + 1
                        Anew(inew,jnew) = A(roworder(i),roworder(j))

                    ENDIF
                ENDDO

            ENDIF

        ENDDO

        !====================================================================!
        ! solve new LSE
#ifdef __MKL
        CALL getrf(Anew,indxnew,info)
#else
        CALL dgetrf(nnew,nnew,Anew,nnew,indxnew,info)
#endif
        IF (info .NE. 0) THEN
            !get the value of info
            WRITE(cbuf,'(A,I2)') ' getrf new failed with info = ', info

            IF (info.GT.0) THEN
                OPEN(81,FILE='DumpA.dat',FORM='FORMATTED',STATUS='REPLACE',IOSTAT=info)
                DO i=1,n
                    WRITE(81,'(10E20.12)') A(i,:)
                ENDDO
                CLOSE(81)
            ENDIF

            !print the error message
            CALL ppm_write(ppm_rank,caller,cbuf,info)
            info = -1
            GOTO 9999
        ENDIF

        DO k=1,n_eq
#ifdef __MKL
            CALL getrs(Anew,indxnew,bnew_n(:,k),'N',info)
#else
            CALL dgetrs('N',nnew,1,Anew,nnew,indxnew,bnew_n(:,k),nnew,info)
#endif
            IF (info .NE. 0) THEN
                WRITE(cbuf,'(A,I2)') ' getrs new failed with info = ', info
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                info = -1
                GOTO 9999
            ENDIF
        ENDDO

        !====================================================================!
        ! sort solution (coefficients) into old places, remaining places are
        ! set to zero 
        jnew = 0
        DO j = 1,n
            x_or_b(j,1:n_eq) = 0._MK
        ENDDO
        DO j = 1,n

            IF (valid(j) .EQ. 1) THEN
                jnew = jnew + 1
                DO k=1,n_eq
                    x_or_b(roworder(j),k) = bnew_n(jnew,k)
                ENDDO
            ENDIF

        ENDDO

        DEALLOCATE(roworder)
        DEALLOCATE(indxnew)
        DEALLOCATE(bnew_n)
        DEALLOCATE(Anew)

    ELSE ! no singularities

        !====================================================================!
        ! solve originial LSE
        DO k=1,n_eq
#ifdef __MKL
            CALL getrs(Acopy,indx,x_or_b(:,k),'N',info)
#else
            CALL dgetrs('N',n,1,Acopy,n,indx,x_or_b(:,k),n,info)
#endif
            IF (info .NE. 0) THEN
                WRITE(cbuf,'(A,I2)') ' getrs new failed with info = ', info
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                info = -1
                GOTO 9999
            ENDIF
        ENDDO

    ENDIF

    !=======================================================================!
    ! check if equations were solved right (if you trust in LAPACK and this 
    ! code PUT THIS IN IF BRANCH FOR THE REDUCED SYSTEM)
    IF (check) THEN
        DO k=1,n_eq
#ifdef __MKL
            CALL gemv(A,x_or_b(:,k),real_b,1._MK,0._MK)
#else
            CALL dgemv('N',SIZE(A,1),SIZE(A,2),1._MK,A,&
                SIZE(A,1),x_or_b(:,k),1,0._MK,real_b,1)
#endif
            
            IF (SUM(ABS(real_b - exact_b_n(:,k))) .GT. tolerance_LSE) THEN
                WRITE(*,*)'WARNING from ', TRIM(caller),': no solution!'
                CALL ppm_write(ppm_rank,caller,&
                    'Error in moment conditions 1 too large.',info)
                WRITE(cbuf,'(A,E12.3)')'The l1-norm of the error is ',&
                    SUM(ABS(real_b - exact_b_n(:,k)))
                CALL ppm_write(ppm_rank,caller,cbuf,info)
                DO i=1,n
                    DO j=1,n
                        WRITE(9001,'(1(E30.22))') A(i,j)
                    ENDDO
                ENDDO
                CALL ppm_write(ppm_rank,caller,&
                    'moment matrix written in file fort.9001',info)
                info = -1
                GOTO 9999
            ENDIF
        ENDDO
    ENDIF

    !=======================================================================!
    ! dealloc
    DEALLOCATE(valid)
    DEALLOCATE(indx)
    DEALLOCATE(Acopy)
    IF (check) THEN
        DEALLOCATE(exact_b_n)
        DEALLOCATE(real_b)
    ENDIF

    9999 CONTINUE ! jump here upon error


#if    __KIND == __SINGLE_PRECISION
END SUBROUTINE solveLSE_ns
#elif  __KIND == __DOUBLE_PRECISION
END SUBROUTINE solveLSE_nd
#endif

