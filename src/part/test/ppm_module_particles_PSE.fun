test_suite ppm_module_particles_PSE

USE ppm_module_particles_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_operator_typedef
USE ppm_module_interfaces
USE ppm_module_data

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: tol=EPSILON(1.0_MK)*100
REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
REAL(MK),PARAMETER              :: skin = 0.0_MK
INTEGER,PARAMETER               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc,topoid
INTEGER                         :: np_global = 3000
REAL(MK),PARAMETER              :: cutoff = 0.15_MK
REAL(MK),DIMENSION(:,:),POINTER :: xp=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: min_phys=>NULL(),max_phys=>NULL()
REAL(MK),DIMENSION(:  ),POINTER :: len_phys=>NULL()
INTEGER                         :: i,j,k,ip,wp_id
INTEGER                         :: nstep
INTEGER,DIMENSION(3)            :: ldc
INTEGER, DIMENSION(6)           :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost=>NULL()
INTEGER                         :: isymm = 0
REAL(MK)                        :: t0,t1,t2,t3
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),ALLOCATABLE :: seed
INTEGER, DIMENSION(:),POINTER   :: nvlist=>NULL()
INTEGER, DIMENSION(:,:),POINTER :: vlist=>NULL()
REAL(MK)                        :: err

INTEGER, DIMENSION(:), POINTER                 :: wp_1i => NULL()
INTEGER, DIMENSION(:,:), POINTER               :: wp_2i => NULL()
INTEGER(ppm_kind_int64),DIMENSION(:),  POINTER :: wp_1li => NULL()
INTEGER(ppm_kind_int64),DIMENSION(:,:),POINTER :: wp_2li => NULL()
REAL(MK), DIMENSION(:),   POINTER              :: wp_1r => NULL()
REAL(MK), DIMENSION(:,:), POINTER              :: wp_2r => NULL()
COMPLEX(MK), DIMENSION(:),   POINTER           :: wp_1c => NULL()
COMPLEX(MK), DIMENSION(:,:), POINTER           :: wp_2c => NULL()
LOGICAL, DIMENSION(:),   POINTER               :: wp_1l => NULL()

INTEGER, DIMENSION(:),ALLOCATABLE              :: degree,order
REAL(ppm_kind_double),DIMENSION(:),ALLOCATABLE :: coeffs
INTEGER                                        :: nterms

    init

        USE ppm_module_init
        USE ppm_module_mktopo

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
            seed(i)=10+i*i*(rank+1)
        ENDDO
        CALL RANDOM_SEED(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    end init


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,len_phys)

    end finalize


    setup


    end setup


    teardown

    end teardown

    test PSE_client
        TYPE(ppm_t_particles_d),TARGET  :: Part1
        TYPE(ppm_t_field)               :: Field1
        TYPE(ppm_t_field)               :: Field2
        TYPE(ppm_t_operator)            :: Laplacian
        CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
        CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
        TYPE(ppm_t_options_op)          :: opts_op

        start_subroutine("test_PSE")
        !--------------------------
        !Define Fields
        !--------------------------
        CALL Field1%create(5,info,name="Concentration") !vector field
        Assert_Equal(info,0)

        CALL Part1%initialize(np_global,info,topoid=topoid)
        Assert_Equal(info,0)

        CALL Part1%create_neighlist(Part1,info)
        Assert_Equal(info,0)

        CALL Part1%set_cutoff(3.0_MK * Part1%h_avg,info)
        Assert_Equal(info,0)

        ALLOCATE(wp_2r(ndim,Part1%Npart))
        CALL RANDOM_NUMBER(wp_2r)
        wp_2r = (wp_2r - 0.5_MK) * Part1%h_avg * 0.15_MK
        CALL Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        DEALLOCATE(wp_2r)
        CALL Part1%apply_bc(info)
        Assert_Equal(info,0)
        CALL Part1%map(info,global=.TRUE.,topoid=topoid)
        Assert_Equal(info,0)

        CALL Field1%discretize_on(Part1,info)
        Assert_Equal(info,0)


        !Initialize concentrations
        foreach p in particles(Part1) with positions(x) vec_fields(w=Field1)
            w_p(1) = 1.0_MK
            w_p(2) = -10.0_MK
            w_p(3) = COS(x_p(1))
            w_p(4) = SIN(x_p(2))
            w_p(5) = SQRT(SUM(x_p(1:ndim)**2))
        end foreach

        !CALL Part1%print_info(info)

        CALL Part1%map_ghosts(info)
        Assert_Equal(info,0)

        CALL Part1%comp_neighlist(info)
        Assert_Equal(info,0)


        nterms=ndim
        ALLOCATE(degree(nterms*ndim),coeffs(nterms),order(nterms))
        IF (ndim .EQ. 2) then
               degree =  (/2,0,   0,2/)
        else
               degree =  (/2,0,0, 0,2,0, 0,0,2/)
        ENDIF
        coeffs = 1.0_MK

        CALL Laplacian%create(ndim,coeffs,degree,info,name="Laplacian")
        Assert_Equal(info,0)

        CALL opts_op%create(ppm_param_op_dcpse,info,order=2,c=0.5D0)
            or_fail("failed to initialize option object for operator")

        CALL Laplacian%discretize_on(Part1,DCop,opts_op,info)
        Assert_Equal(info,0)
        Assert_True(ASSOCIATED(DCop))
        !CALL Laplacian%discretize_on(Part1,PSEop,info,method="PSE")
        !Assert_Equal(info,0)
        !Assert_True(ASSOCIATED(PSEop))

        CALL DCop%compute(Field1,Field2,info)
        Assert_Equal(info,0)


        CALL DCop%destroy(info)
        Assert_Equal(info,0)
        CALL Laplacian%destroy(info)
        Assert_Equal(info,0)
        CALL Part1%destroy(info)
        Assert_Equal(info,0)
        CALL Field1%destroy(info)
        Assert_Equal(info,0)
        DEALLOCATE(degree,coeffs,order)
        CALL opts_op%destroy(info)
            or_fail("failed to destroy opts_op")
        end_subroutine()
    end test
!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
PURE FUNCTION f0_test(pos,ndim)

    REAL(MK)                              :: f0_test
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    f0_test =  SIN(2._MK*pi*pos(1)) * COS(2._MK*pi*pos(2)) * &
        & SIN(2._MK*pi*pos(ndim))

    return

END FUNCTION f0_test



end test_suite
