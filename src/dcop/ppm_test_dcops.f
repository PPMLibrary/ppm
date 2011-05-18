!!!
!!! Routine that tests whether the DC operators are accurate.
!!!
#ifdef __2D
SUBROUTINE ppm_test_dcops_2d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        nneighmin,nneighmax,c,order_deriv,order_approx,f0_fun,df0_fun,&
        info,islaplacian,iscartesian)
#endif
#ifdef __3D
SUBROUTINE ppm_test_dcops_3d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        nneighmin,nneighmax,c,order_deriv,order_approx,f0_fun,df0_fun,&
        info,islaplacian,iscartesian)
#endif

    USE ppm_module_data, ONLY: ppm_mpi_kind,ppm_comm
    IMPLICIT NONE
#ifdef __MPI
#include "mpif.h"
#endif

#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
    INTEGER, PARAMETER :: MK = ppm_kind_single
#else
    INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

    ! arguments
    REAL(MK), DIMENSION(:,:),POINTER,    INTENT(IN   )   :: xp
    REAL(MK), DIMENSION(:),  POINTER,    INTENT(IN   )   :: rcp
    INTEGER,                             INTENT(IN   )   :: Npart
    INTEGER,                             INTENT(IN   )   :: Mpart
    INTEGER,                             INTENT(IN   )   :: nneighmin
    INTEGER,                             INTENT(IN   )   :: nneighmax
    INTEGER,  DIMENSION(:),  POINTER,    INTENT(IN   )   :: nvlist
    INTEGER,  DIMENSION(:,:),POINTER,    INTENT(IN   )   :: vlist
    REAL(MK),                            INTENT(IN   )   :: c
    INTEGER,DIMENSION(ppm_dim),          INTENT(IN   )   :: order_deriv
    !degree of the derivative that we want to approximate
    INTEGER,OPTIONAL,                    INTENT(IN   )   :: order_approx
    !order of the approximation
    INTEGER,                             INTENT(  OUT)   :: info

    ! Optional arguments
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: islaplacian
    !special treatment for the laplacian operator
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: iscartesian
    !special case where all particles are on a cartesian grid

! argument-functions need an interface
    INTERFACE
        !Test function
        FUNCTION f0_fun(x)
            USE ppm_module_data, ONLY:ppm_dim
            USE ppm_module_typedef
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
            INTEGER, PARAMETER :: MK = ppm_kind_single
#else
            INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK)                                         :: f0_fun
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: x
        END FUNCTION f0_fun
        !Exact derivative
        FUNCTION df0_fun(x,orderderiv)
            USE ppm_module_data, ONLY:ppm_dim
            USE ppm_module_typedef
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
            INTEGER, PARAMETER :: MK = ppm_kind_single
#else
            INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            REAL(MK)                                         :: df0_fun
            REAL(MK),DIMENSION(ppm_dim),INTENT(IN)           :: x
            INTEGER, DIMENSION(ppm_dim),INTENT(IN),OPTIONAL  :: orderderiv
        END FUNCTION df0_fun
    END INTERFACE


    ! local variables
    INTEGER                               :: i,j,ip,iq,ineigh
    INTEGER                               :: di,direction
    CHARACTER (LEN = 256)                 :: cbuf,caller='ppm_test_dcops'
    REAL(MK),DIMENSION(49)                :: kh,sum1
    INTEGER                               :: beta,sig
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: lambda_max
    REAL(MK)                              :: norm_inf,norm_L2
    REAL(MK)                              :: err_inf,err_L2
    REAL(MK), DIMENSION(:,:),ALLOCATABLE  :: eta
    REAL(MK), DIMENSION(:)  ,ALLOCATABLE  :: err
    REAL(MK), DIMENSION(:)  ,ALLOCATABLE  :: phi
    REAL(MK), DIMENSION(:)  ,ALLOCATABLE  :: phi_deriv
    LOGICAL                               :: laplacian,cartesian


    !-------------------------------------------------------------------------
    !  Initialise
    !-------------------------------------------------------------------------
    info = 0
    CALL substart(caller,t0,info)
    laplacian = .FALSE.
    IF (PRESENT(islaplacian)) THEN
        IF (islaplacian) THEN
            laplacian = .TRUE.
        ENDIF
    ENDIF
    cartesian = .FALSE.
    IF (PRESENT(iscartesian)) THEN
        IF (iscartesian) THEN
            cartesian = .TRUE.
        ENDIF
    ENDIF

    IF(laplacian) THEN
        WRITE(cbuf,'(A,A,A,I2,A,L)') 'Testing - Degree of derivative: ','laplacian',&
            ' Approximation order: ',order_approx,' Cartesian = ',cartesian
    ELSE
#ifdef __2D
        WRITE(cbuf,'(A,2(I2,1X),A,I2,A,L)') 'Testing - Degree of derivative: ',order_deriv,&
            ' Approximation order: ',order_approx,' Cartesian = ',cartesian
#else
        WRITE(cbuf,'(A,3(I2,1X),A,I2,A,L)') 'Testing - Degree of derivative: ',order_deriv,&
            ' Approximation order: ',order_approx,' Cartesian = ',cartesian
#endif
    ENDIF
    CALL pwrite(ppm_rank,caller,cbuf,info)

    !-------------------------------------------------------------------------
    !  Compute DC operators
    !-------------------------------------------------------------------------
#ifdef __2D
    CALL ppm_dcops_2d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#else
    CALL ppm_dcops_3d(xp,rcp,Npart,Mpart,nvlist,vlist,&
        eta,c,nneighmin,nneighmax,info,&
        order_deriv,order_approx,islaplacian,iscartesian)
#endif
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'ppm_dcops failed',info)
        info = -1
        GOTO 9999
    ENDIF

    !-------------------------------------------------------------------------
    !  Allocate arrays for derivatives and errors
    !-------------------------------------------------------------------------
    ALLOCATE(phi(Mpart),phi_deriv(Npart),err(Npart),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    DO ip=1,Mpart
            phi(ip) = f0_fun(xp(1:ppm_dim,ip))
    ENDDO

    !-------------------------------------------------------------------------
    !  Compute derivatives using DC operators
    !-------------------------------------------------------------------------
    beta = SUM(order_deriv)
    IF (laplacian) beta = 2

    sig = 2*MOD(beta,2)-1 ! positive/negative if beta is odd/even
    phi_deriv = 0._MK

    DO ip=1,Npart
        DO ineigh=1,nvlist(ip)
            iq=vlist(ineigh,ip)
            phi_deriv(ip) = phi_deriv(ip)+        &
                (phi(iq)+sig*phi(ip))*eta(ineigh,ip)
        ENDDO
    ENDDO

    !-------------------------------------------------------------------------
    !  Compute errors
    !-------------------------------------------------------------------------
    DO ip=1,Npart
        err(ip) = df0_fun(xp(1:ppm_dim,ip),order_deriv)
    ENDDO

    norm_inf=MAXVAL(err(1:Npart))
    norm_L2=SUM(err(1:Npart)**2)
    DO ip=1,Npart
        err(ip) = phi_deriv(ip) - df0_fun(xp(1:ppm_dim,ip),order_deriv)
    ENDDO
    err_inf=MAXVAL(ABS(err(1:Npart)))
    err_L2=SUM(err(1:Npart)**2)
    
#ifdef __MPI
    CALL MPI_Allreduce(err_inf,err_inf,1,ppm_mpi_kind,MPI_SUM,ppm_comm,info)
    CALL MPI_Allreduce(err_L2,err_L2,1,ppm_mpi_kind,MPI_SUM,ppm_comm,info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'MPI_Allreduce failed',info)
        info=-1
        GOTO 9999
    ENDIF
#endif
    err_inf = err_inf / norm_inf
    err_L2 = SQRT(err_L2 / norm_L2)

    WRITE(cbuf,'(A,E12.6,A,E12.6)') 'L_infinity error: ',err_inf,' L2 error: ',err_L2
    CALL pwrite(ppm_rank,caller,cbuf,info)


    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error


#ifdef __2D
END SUBROUTINE ppm_test_dcops_2d
#endif
#ifdef __3D
END SUBROUTINE ppm_test_dcops_3d
#endif
