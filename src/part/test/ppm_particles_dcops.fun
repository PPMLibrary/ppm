test_suite ppm_particles_dcops

USE ppm_module_interfaces
USE ppm_module_particles_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_operator_typedef
USE ppm_module_sop
USE ppm_module_data
USE ppm_module_io_vtk

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: tol=EPSILON(1.0_MK)*100
REAL(MK),PARAMETER              :: pi = aCOS(-1.0_MK)
INTEGER,PARAMETER               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc,topoid
INTEGER                         :: np_global = 3000
REAL(MK),PARAMETER              :: cutoff = 0.15_MK
REAL(MK),DIMENSION(:,:),POINTER :: xp=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,len_phys
INTEGER                         :: i,j,k,ip,nterms
INTEGER                         :: wp1_id=0, dwp1_id=0, wp2_id=0, op_id=0
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost
TYPE(ppm_t_particles_d),TARGET  :: Part1
TYPE(ppm_t_sop_d),      TARGET  :: Part1_a
TYPE(ppm_t_field)               :: SField1,SField2,SField3,SField4
TYPE(ppm_t_field)               :: VFieldD,VField2,VField3,VField5
INTEGER                         :: seedsize
INTEGER, DIMENSION(:),POINTER   :: nvlist=>NULL()
INTEGER, DIMENSION(:,:),POINTER :: vlist=>NULL()
INTEGER, DIMENSION(:),ALLOCATABLE:: seed,degree,order
REAL(MK),DIMENSION(:,:),ALLOCATABLE:: randn
REAL(MK),DIMENSION(:),ALLOCATABLE:: coeffs
INTEGER, DIMENSION(:),  POINTER :: gi => NULL()
REAL(MK),DIMENSION(:),  POINTER :: wp_1r => NULL()
REAL(MK),DIMENSION(:,:),POINTER :: wp_2r => NULL()
REAL(MK)                        :: tol_error,err
TYPE(ppm_t_operator)            :: Op
TYPE(ppm_t_options_op)          :: opts_op
CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
CLASS(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
CLASS(ppm_t_discr_data),POINTER :: prop => NULL()

    init

        USE ppm_module_init
        USE ppm_module_mktopo
        USE ppm_module_util_qsort
        start_subroutine("init")

        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),STAT=info)

        min_phys(1:ndim) = 0.0_MK
        max_phys(1:ndim) = 1.0_MK
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1.0_MK)))+10
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

        CALL RANDOM_SEED(size=seedsize)
        ALLOCATE(seed(seedsize))
        DO i=1,seedsize
            seed(i)=10+i*i !*(rank+1)
        ENDDO
        CALL RANDOM_SEED(put=seed)
        ALLOCATE(randn(ndim,np_global))
        CALL RANDOM_NUMBER(randn)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        !--------------------------
        !Define Fields
        !--------------------------
        CALL SField1%create(1,info,name="F_sca1") !scalar field
        CALL SField2%create(1,info,name="F_sca2") !scalar field
        CALL SField3%create(1,info,name="F_sca3") !scalar field

        CALL VFieldD%create(ndim,info,name="F_vecDim") !vector field
        CALL VField2%create(2,info,name="F_vec2") !vector field
        CALL VField3%create(3,info,name="F_vec3") !vector field

        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)

        !initialize particles on a grid
        CALL Part1%initialize(np_global,info,topoid=topoid, &
        distrib=ppm_param_part_init_cartesian)

        CALL Part1%create_neighlist(Part1,info)

        CALL Part1%comp_global_index(info)

        CALL Part1%map(info,global=.TRUE.,topoid=topoid)

        CALL Part1%map_ghosts(info)

        CALL SField1%discretize_on(Part1,info)
        CALL SField2%discretize_on(Part1,info)
        CALL SField3%discretize_on(Part1,info)

        CALL VFieldD%discretize_on(Part1,info)
        CALL VField2%discretize_on(Part1,info)
        CALL VField3%discretize_on(Part1,info)

        !Perturb the  particles positions
        CALL Part1%set_cutoff(3.0_MK * Part1%h_avg,info)

        ALLOCATE(wp_2r(ndim,Part1%Npart))

        CALL Part1%get(Part1%gi,gi,info,read_only=.TRUE.)

        FORALL (ip=1:Part1%Npart) wp_2r(1:ndim,ip) = randn(1:ndim,gi(ip))

        wp_2r = 0.5_MK * (wp_2r - 0.5_MK) * Part1%h_avg

        CALL Part1%move(wp_2r,info)

        DEALLOCATE(wp_2r)

        CALL Part1%apply_bc(info)
        CALL Part1%map(info)

        foreach p in particles(Part1) with positions(x) sca_fields(S1=SField1,S2=SField2,S3=SField3) vec_fields(VD=VFieldD,V2=VField2,V3=VField3)
            S1_p = f0_test(x_p(1:ndim),ndim)
            S2_p = 0.0_MK !f0_test(x_p(1:ndim),ndim)
            S3_p = -8.0_MK * pi*pi * f0_test(x_p(1:ndim),ndim)
            VD_p(1:ndim)    = -10.0_MK
            V2_p(1:2)       = -10.0_MK
            V3_p(1:3)       = -10.0_MK
        end foreach

        end_subroutine()
    end init


    finalize
        USE ppm_module_finalize

        CALL DCop%destroy(info)
        CALL Op%destroy(info)
        CALL Part1%destroy(info)
        CALL SField1%destroy(info)
        CALL SField2%destroy(info)
        CALL SField3%destroy(info)
        CALL SField4%destroy(info)
        CALL VFieldD%destroy(info)
        CALL VField2%destroy(info)
        CALL VField3%destroy(info)
!        CALL VField5%destroy(info)
        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,len_phys)
        DEALLOCATE(randn)

    end finalize


    setup

    end setup


    teardown

        IF (ALLOCATED(degree)) DEALLOCATE(degree,coeffs,order)

    end teardown


    test Laplacian
        start_subroutine("test_laplacian")
        tol_error = 2e-2

        CALL Part1%set_cutoff(3.0_MK * Part1%h_avg,info)
        Assert_Equal(info,0)

        CALL Part1%map_ghosts(info)
        Assert_Equal(info,0)

        CALL Part1%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        ALLOCATE(degree(nterms*ndim),coeffs(nterms),order(nterms))
        IF (ndim .EQ. 2) THEN
               degree =  (/2,0,   0,2/)
        else
               degree =  (/2,0,0, 0,2,0, 0,0,2/)
        ENDIF
        coeffs = 1.0_MK

        CALL Op%create(ndim,coeffs,degree,info,name="Laplacian")
        Assert_Equal(info,0)

        CALL opts_op%create(ppm_param_op_dcpse,info,order=2,c=1.4D0)
        or_fail("failed to initialize option object for operator")


        CALL Op%discretize_on(Part1,DCop,opts_op,info)
        Assert_Equal(info,0)
        Assert_True(ASSOCIATED(DCop))
        !CALL Op%discretize_on(Part1,PSEop,info,method="PSE")
        !Assert_Equal(info,0)
        !Assert_True(ASSOCIATED(PSEop))

        CALL DCop%compute(SField1,SField2,info)
        Assert_Equal(info,0)
        !testing output on a field that is not yet created nor discretized
        CALL DCop%compute(SField1,SField4,info)
        Assert_Equal(info,0)

        !CALL ppm_vtk_particles("output",Part1,info)

        Assert_True(inf_error(Part1,SField1,SField2,DCop).LT.tol_error)
        end_subroutine()
    end test

    test Gradient
        start_subroutine("test_gradient")

        tol_error = 1e-2

        CALL Part1%set_cutoff(3.0_MK * Part1%h_avg,info)
        Assert_Equal(info,0)

        CALL Part1%map_ghosts(info)
        Assert_Equal(info,0)

        CALL Part1%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        ALLOCATE(degree(nterms*ndim),coeffs(nterms),order(nterms))

        IF (ndim .EQ. 2) THEN
               degree =  (/1,0,   0,1/)
        else
               degree =  (/1,0,0, 0,1,0, 0,0,1/)
        ENDIF
        coeffs = 1.0_MK

        CALL Op%create(ndim,coeffs,degree,info,name="Gradient")
        Assert_Equal(info,0)

        CALL opts_op%create(ppm_param_op_dcpse,info,order=2,&
        c=1.4D0,vector=.TRUE.)
        or_fail("failed to initialize option object for operator")

        CALL Op%discretize_on(Part1,DCop,opts_op,info)
        Assert_Equal(info,0)
        Assert_True(ASSOCIATED(DCop))

        CALL DCop%compute(SField1,VFieldD,info)
        Assert_Equal(info,0)

        Assert_True(inf_error(Part1,SField1,VFieldD,DCop).LT.tol_error)

        end_subroutine()
    end test
!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
PURE FUNCTION f0_test(pos,ndim)

    REAL(MK)                              :: f0_test
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    IF (ndim .EQ. 2) THEN
        f0_test =  SIN(2._MK*pi*pos(1)) * COS(2._MK*pi*pos(2))
    else
        f0_test =  SIN(2._MK*pi*pos(1)) * COS(2._MK*pi*pos(2)) * &
            & SIN(2._MK*pi*pos(3))
    ENDIF

    return

END FUNCTION f0_test


!-------------------------------------------------------------
! derivatives of the test function
!-------------------------------------------------------------
PURE FUNCTION df0_test(pos,order_deriv,ndim)

    REAL(MK)                                  :: df0_test
    INTEGER                     , INTENT(IN)  :: ndim
    REAL(MK), DIMENSION(ppm_dim), INTENT(IN)  :: pos
    INTEGER,  DIMENSION(ppm_dim), INTENT(IN)  :: order_deriv

    select case (order_deriv(1))
    case (0)
        df0_test =    sin (2._MK*pi * pos(1))
    case (1)
        df0_test =    2._MK*pi*COS(2._MK*pi * pos(1))
    case (2)
        df0_test =  -(2._MK*pi)**2*SIN(2._MK*pi * pos(1))
    case (3)
        df0_test =  -(2._MK*pi)**3*COS(2._MK*pi * pos(1))
    case (4)
        df0_test =   (2._MK*pi)**4*SIN(2._MK*pi * pos(1))
    case default
        df0_test =  0.0_MK
    endselect

    select case (order_deriv(2))
    case (0)
        df0_test =   df0_test * cos (2._MK*pi * pos(2))
    case (1)
        df0_test =   df0_test * (-2._MK*pi)*SIN(2._MK*pi * pos(2))
    case (2)
        df0_test =  df0_test * (-(2._MK*pi)**2)*COS(2._MK*pi * pos(2))
    case (3)
        df0_test =  df0_test * ( (2._MK*pi)**3)*SIN(2._MK*pi * pos(2))
    case (4)
        df0_test =  df0_test * ( (2._MK*pi)**4)*COS(2._MK*pi * pos(2))
    case default
        df0_test =  0.0_MK
    endselect

    IF (ndim .EQ. 3 ) THEN
        select case (order_deriv(3))
        case (0)
            df0_test =   df0_test * sin (2._MK*pi * pos(3))
        case (1)
            df0_test =   df0_test * (2._MK*pi)*COS(2._MK*pi * pos(3))
        case (2)
            df0_test =  df0_test * (-(2._MK*pi)**2)*SIN(2._MK*pi * pos(3))
        case (3)
            df0_test =  df0_test * (-(2._MK*pi)**3)*COS(2._MK*pi * pos(3))
        case (4)
            df0_test =  df0_test * ( (2._MK*pi)**4)*SIN(2._MK*pi * pos(3))
        case default
            df0_test =  0.0_MK
        endselect
    ENDIF

END FUNCTION df0_test


!-------------------------------------------------------------
! Compute the infinity norm of the error
!-------------------------------------------------------------
function inf_error(Part1,Field1,Field2,Op)
    USE ppm_module_data
    TYPE(ppm_t_particles_d),TARGET   :: Part1
    TYPE(ppm_t_field)                :: Field1,Field2
    CLASS(ppm_t_operator_discr)       :: Op
    INTEGER                          :: ip,nterms
    REAL(MK), DIMENSION(:),  POINTER :: wp_1 => NULL(), dwp_1 => NULL()
    REAL(MK), DIMENSION(:,:),POINTER :: xp=>NULL(), wp_2=>NULL(), dwp_2=>NULL()
    REAL(MK), DIMENSION(:), ALLOCATABLE :: err,exact
    REAL(MK)                         :: linf,inf_error,coeff
    LOGICAL                          :: input_vec,output_vec
    INTEGER,DIMENSION(:),ALLOCATABLE :: degree,order
    INTEGER,DIMENSION(ndim)          :: dg
    character(len=100)               :: fname

    IF (op%flags(ppm_ops_vector)) THEN
        CALL Part1%get(Field2,dwp_2,info,read_only=.TRUE.)
        IF (info.NE.0) WRITE(*,*) "Failed te get Field2"
        output_vec = .TRUE.
    else
        CALL Part1%get(Field2,dwp_1,info,read_only=.TRUE.)
        IF (info.NE.0) WRITE(*,*) "Failed te get Field2"
        output_vec = .FALSE.
    ENDIF

    IF (Field1%lda.EQ.1) THEN
        CALL Part1%get(Field1,wp_1,info,read_only=.TRUE.)
        IF (info.NE.0) WRITE(*,*) "Failed te get Field1"
        input_vec = .FALSE.
    else
        CALL Part1%get(Field1,wp_2,info,read_only=.TRUE.)
        IF (info.NE.0) WRITE(*,*) "Failed te get Field1"
        input_vec = .TRUE.
    ENDIF

    nterms = op%op_ptr%nterms
    ALLOCATE(err(nterms),exact(nterms),degree(nterms*ndim),order(nterms))

    CALL Part1%get_xp(xp,info)
        IF (info.NE.0) WRITE(*,*) "Could not get xp"

    err = 0.0_MK
    linf = 0.0_MK
    DO ip=1,Part1%Npart
        exact = 0.0_MK
        DO i=1,nterms
            coeff = op%op_ptr%coeffs(i)
            dg = op%op_ptr%degree(1+(i-1)*ndim:i*ndim)
            IF (output_vec) THEN
                exact(i) = coeff*df0_test(xp(1:ndim,ip),dg,ndim)
            else
                exact(1) = exact(1) + coeff*df0_test(xp(1:ndim,ip),dg,ndim)
            ENDIF
        ENDDO
        IF (output_vec) THEN
            DO i=1,nterms
                err(i) = MAX(err(i),ABS(dwp_2(i,ip) - exact(i)))
                !REMOVME
                !dwp_2(i,ip) = err(i)
            ENDDO
        else
            err(1) = MAX(err(1),ABS(dwp_1(ip) - exact(1)))
        ENDIF
        linf = MAX(linf,MAXVAL(ABS(exact)))
    ENDDO

#ifdef __MPI
    CALL MPI_Allreduce(linf,linf,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
    CALL MPI_Allreduce(MAXVAL(err),inf_error,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
#else
    inf_error = MAXVAL(err)
#endif
    inf_error = inf_error/linf

    IF (ppm_rank.EQ.0) &
        WRITE(*,*) '[',ppm_rank,']','Error is ',inf_error

    DEALLOCATE(err,exact,degree,order)
END FUNCTION inf_error

end test_suite
