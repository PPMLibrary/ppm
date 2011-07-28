!!!
!!! Routine that tests whether the DC operators are accurate.
!!!
#ifdef __2D
SUBROUTINE ppm_test_part_dcops_2d(Particles,eta_id,c,order_deriv,order_approx,&
        f0_fun,df0_fun,info,islaplacian,iscartesian,isinterp,&
        Particles_old,nn_sq_id)
#endif
#ifdef __3D
SUBROUTINE ppm_test_part_dcops_3d(Particles,eta_id,c,order_deriv,order_approx,&
        f0_fun,df0_fun,info,islaplacian,iscartesian,isinterp,&
        Particles_old,nn_sq_id)
#endif

    USE ppm_module_data, ONLY: ppm_mpi_kind,ppm_comm
    USE ppm_module_particles
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
    TYPE(ppm_t_particles), POINTER,      INTENT(INOUT)   :: Particles
    !!! Particles
    INTEGER,                             INTENT(INOUT)   :: eta_id
    !!! id where eta is stored
    REAL(MK),                            INTENT(IN   )   :: c
    !!! ratio h/epsilon
    INTEGER,DIMENSION(ppm_dim),          INTENT(IN   )   :: order_deriv
    !!! degree of the derivative that we want to approximate
    INTEGER,                             INTENT(IN   )   :: order_approx
    !!! order of the approximation
    INTEGER,                             INTENT(  OUT)   :: info

    ! Optional arguments
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: islaplacian
    !!! special treatment for the laplacian operator
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: iscartesian
    !!! special case where all particles are on a cartesian grid
    LOGICAL,OPTIONAL,                    INTENT(IN   )   :: isinterp
    !!! special case for interpolating operators 
    TYPE(ppm_t_particles),OPTIONAL,POINTER,INTENT(INOUT) :: Particles_old
    !!! second set of Particles
    INTEGER,OPTIONAL,                    INTENT(IN   )   :: nn_sq_id
    !!! id where the nearest-neighbour distances within Particles_old
    !!! are stored

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
    CHARACTER (LEN = 256)                 :: cbuf,caller='ppm_test_part_dcops'
    REAL(MK),DIMENSION(49)                :: kh,sum1
    INTEGER                               :: beta,sig
    REAL(KIND(1.D0))                      :: t0
    REAL(MK)                              :: lambda_max
    REAL(MK)                              :: norm_inf,norm_L2
    REAL(MK)                              :: err_inf,err_L2
    REAL(MK), DIMENSION(:)  ,ALLOCATABLE  :: err
    REAL(MK), DIMENSION(:)  ,POINTER      :: phi1=>NULL()
    REAL(MK), DIMENSION(:)  ,POINTER      :: phi2=>NULL()
    REAL(MK), DIMENSION(:)  ,ALLOCATABLE  :: phi_deriv
    LOGICAL                               :: laplacian,cartesian
    REAL(MK), DIMENSION(:,:),POINTER      :: eta=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: xp1=>NULL()
    REAL(MK), DIMENSION(:,:),POINTER      :: xp2=>NULL()
    REAL(MK), DIMENSION(:),  POINTER      :: rcp=>NULL()
    INTEGER,  DIMENSION(:,:),POINTER      :: vlist=>NULL()
    INTEGER,  DIMENSION(:  ),POINTER      :: nvlist=>NULL()
    INTEGER                               :: Npart,Mpart


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
    IF (isinterp .AND. .NOT. PRESENT(Particles_old)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Not enough arguments - need Particles_old for interpolation',&
            __LINE__,info)
        GOTO 9999
    ENDIF
    IF (isinterp .AND. .NOT. PRESENT(nn_sq_id)) THEN
        info = ppm_error_error
        CALL ppm_error(ppm_err_argument,caller,&
            'Not enough arguments - need nn_sq_id for interpolation',&
            __LINE__,info)
        GOTO 9999
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
    CALL ppm_part_dcops_2d(Particles,eta_id,c,info,&
        order_deriv,order_approx,islaplacian,iscartesian,&
        isinterp,Particles_old,nn_sq_id)
#else
    CALL ppm_part_dcops_3d(Particles,eta_id,c,info,&
        order_deriv,order_approx,islaplacian,iscartesian,&
        isinterp,Particles_old,nn_sq_id)
#endif
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'ppm_dcops failed',info)
        info = -1
        GOTO 9999
    ENDIF


    Npart = Particles%Npart
    Mpart = Particles%Mpart
    IF (isinterp) THEN
        nvlist => Particles%nvlist_cross
        vlist  => Particles%vlist_cross
        xp1 => Get_xp(Particles)
        xp2 => Get_xp(Particles_old,with_ghosts=.TRUE.)
    ELSE
        nvlist => Particles%nvlist
        vlist  => Particles%vlist
        xp1 => Get_xp(Particles,with_ghosts=.TRUE.)
        xp2 => xp1
    ENDIF
    eta => Get_wpv(Particles,eta_id)

    !-------------------------------------------------------------------------
    !  Allocate arrays for derivatives and errors
    !-------------------------------------------------------------------------
    ALLOCATE(phi1(Mpart),phi_deriv(Npart),err(Npart),STAT=info)
    IF (info .NE. 0) THEN
        CALL pwrite(ppm_rank,caller,'allocation failed.',info)
        info = -1
        GOTO 9999
    ENDIF
    IF (isinterp) THEN
        ALLOCATE(phi2(Particles_old%Mpart),STAT=info)
    ELSE
        phi2 => phi1
    ENDIF

    DO ip=1,Mpart
            phi1(ip) = f0_fun(xp1(1:ppm_dim,ip))
    ENDDO
    IF (isinterp) THEN
        DO ip=1,Particles_old%Mpart
            phi2(ip) = f0_fun(xp2(1:ppm_dim,ip))
        ENDDO
    ENDIF

    !-------------------------------------------------------------------------
    !  Compute derivatives using DC operators
    !-------------------------------------------------------------------------
    beta = SUM(order_deriv)
    IF (laplacian) beta = 2

    sig = 2*MOD(beta,2)-1 ! positive/negative if beta is odd/even
    phi_deriv = 0._MK

    IF (isinterp .and. SUM(order_deriv).eq.0) sig = 0._MK

    DO ip=1,Npart
        DO ineigh=1,nvlist(ip)
            iq=vlist(ineigh,ip)
            phi_deriv(ip) = phi_deriv(ip)+        &
                (phi2(iq)+sig*phi1(ip))*eta(ineigh,ip)
        ENDDO
    ENDDO

    !-------------------------------------------------------------------------
    !  Compute errors
    !-------------------------------------------------------------------------
    DO ip=1,Npart
        err(ip) = df0_fun(xp1(1:ppm_dim,ip),order_deriv)
    ENDDO

    !DO ip=1,Npart
        !write(199,'(3(E12.4,2X))') xp1(1:ppm_dim,ip),phi1(ip)
        !write(200,'(3(E12.4,2X))') xp1(1:ppm_dim,ip),err(ip)
    !ENDDO

    norm_inf=MAXVAL(err(1:Npart))
    norm_L2=SUM(err(1:Npart)**2)
    DO ip=1,Npart
        err(ip) = phi_deriv(ip) - df0_fun(xp1(1:ppm_dim,ip),order_deriv)
    ENDDO
    err_inf=MAXVAL(ABS(err(1:Npart)))
    err_L2=SUM(err(1:Npart)**2)

    !DO ip=1,Npart
        !write(201,'(3(E12.4,2X))') xp1(1:ppm_dim,ip),phi_deriv(ip)
        !write(202,'(3(E12.4,2X))') xp1(1:ppm_dim,ip),err(ip)
    !ENDDO

    !stop

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
    IF (err_inf .GT. 1._MK) THEN
        CALL pwrite(ppm_rank,caller,' ** FAILED ** ',info)
    ENDIF


    nvlist => NULL()
    vlist => NULL()
    IF (isinterp) THEN
        xp1 => Set_xp(Particles)
        xp2 => Set_xp(Particles_old,read_only=.TRUE.)
    ELSE
        xp1 => Set_xp(Particles,read_only=.TRUE.)
        xp2 => NULL()
    ENDIF
    eta => Set_wpv(Particles,eta_id,read_only=.TRUE.)

    DEALLOCATE(phi1,phi_deriv,err)
    IF (isinterp) DEALLOCATE(phi2)
    phi1 => NULL()
    phi2 => NULL()

    CALL substop(caller,t0,info)

    9999 CONTINUE ! jump here upon error


#ifdef __2D
END SUBROUTINE ppm_test_part_dcops_2d
#endif
#ifdef __3D
END SUBROUTINE ppm_test_part_dcops_3d
#endif
