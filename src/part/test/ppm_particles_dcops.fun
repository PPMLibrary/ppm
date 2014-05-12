test_suite ppm_particles_dcops

use ppm_module_interfaces
use ppm_module_particles_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_operator_typedef
use ppm_module_sop
use ppm_module_data
use ppm_module_io_vtk

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = acos(-1._mk)
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc,topoid
integer                         :: np_global = 3000
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys,len_phys
integer                         :: i,j,k,ip,nterms
integer                         :: wp1_id=0, dwp1_id=0, wp2_id=0, op_id=0
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
type(ppm_t_particles_d),TARGET  :: Part1
type(ppm_t_sop_d),      TARGET  :: Part1_a
type(ppm_t_field)               :: SField1,SField2,SField3,SField4
type(ppm_t_field)               :: VFieldD,VField2,VField3,VField5
integer                         :: seedsize
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
integer, dimension(:),allocatable:: seed,degree,order
real(mk),dimension(:,:),allocatable:: randn
real(mk),dimension(:),allocatable:: coeffs
integer, dimension(:),  pointer :: gi => NULL(), energy => null(),energyt => null()
real(mk),dimension(:),  pointer :: wp_1r => NULL()
real(mk),dimension(:,:),pointer :: wp_2r => NULL()
real(mk)                        :: tol_error,err
type(ppm_t_operator)            :: Op
TYPE(ppm_t_options_op)          :: opts_op
CLASS(ppm_t_operator_discr),POINTER   :: DCop => NULL()
CLASS(ppm_t_operator_discr),POINTER   :: PSEop => NULL()
class(ppm_t_neighlist_d_),POINTER :: Nlist => NULL()
class(ppm_t_discr_data),POINTER :: prop => NULL()

    init

        use ppm_module_init
        use ppm_module_mktopo
        use ppm_module_util_qsort
        start_subroutine("init")

        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = int(log10(epsilon(1._mk)))+10
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i !*(rank+1)
        enddo
        call random_seed(put=seed)
        allocate(randn(ndim,np_global))
        call random_number(randn)

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
        call SField1%create(1,info,name="F_sca1") !scalar field
        call SField2%create(1,info,name="F_sca2") !scalar field
        call SField3%create(1,info,name="F_sca3") !scalar field

        call VFieldD%create(ndim,info,name="F_vecDim") !vector field
        call VField2%create(2,info,name="F_vec2") !vector field
        call VField3%create(3,info,name="F_vec3") !vector field

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)

        !initialize particles on a grid
        call Part1%initialize(np_global,info,topoid=topoid, &
        distrib=ppm_param_part_init_cartesian)

        CALL Part1%create_neighlist(Part1,info)

        call Part1%comp_global_index(info)

        call Part1%map(info,global=.true.,topoid=topoid)

        call Part1%map_ghosts(info)

        if (ppm_rank.eq.0) then
        allocate(energy(5))
        energy(1)=1000
        energy(2)=-150
        energy(3)=300
        energy(4)=120
        energy(5)=50
        stdout("energy1",'energy')
        call ppm_util_qsort(energy,energyt,info)
        stdout("energy2",'energy')
        stdout("energyt",'energyt')
        endif

        CALL MPI_FINALIZE(INFO)
        Stop


        call SField1%discretize_on(Part1,info)
        call SField2%discretize_on(Part1,info)
        call SField3%discretize_on(Part1,info)

        call VFieldD%discretize_on(Part1,info)
        call VField2%discretize_on(Part1,info)
        call VField3%discretize_on(Part1,info)

        !Perturb the  particles positions
        call Part1%set_cutoff(3._mk * Part1%h_avg,info)

        allocate(wp_2r(ndim,Part1%Npart))

        call Part1%get(Part1%gi,gi,info,read_only=.true.)

        FORALL (ip=1:Part1%Npart) wp_2r(1:ndim,ip) = randn(1:ndim,gi(ip))

        wp_2r = 0.5_mk * (wp_2r - 0.5_mk) * Part1%h_avg

        call Part1%move(wp_2r,info)

        deallocate(wp_2r)

        call Part1%apply_bc(info)
        call Part1%map(info)

        foreach p in particles(Part1) with positions(x) sca_fields(S1=SField1,S2=SField2,S3=SField3) vec_fields(VD=VFieldD,V2=VField2,V3=VField3)
            S1_p = f0_test(x_p(1:ndim),ndim)
            S2_p = 0._mk !f0_test(x_p(1:ndim),ndim)
            S3_p = -8._mk * pi*pi * f0_test(x_p(1:ndim),ndim)
            VD_p(1:ndim)    = -10._mk
            V2_p(1:2)       = -10._mk
            V3_p(1:3)       = -10._mk
        end foreach

        end_subroutine()
    end init


    finalize
        use ppm_module_finalize

        call DCop%destroy(info)
        call Op%destroy(info)
        call Part1%destroy(info)
        call SField1%destroy(info)
        call SField2%destroy(info)
        call SField3%destroy(info)
        call SField4%destroy(info)
        call VFieldD%destroy(info)
        call VField2%destroy(info)
        call VField3%destroy(info)
!        call VField5%destroy(info)
        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys)
        deallocate(randn)

    end finalize


    setup

    end setup


    teardown

        if (allocated(degree)) deallocate(degree,coeffs,order)

    end teardown


    test Laplacian
        start_subroutine("test_laplacian")
        tol_error = 2e-2

        call Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/2,0,   0,2/)
        else
               degree =  (/2,0,0, 0,2,0, 0,0,2/)
        endif
        coeffs = 1.0_mk

        call Op%create(ndim,coeffs,degree,info,name="Laplacian")
        Assert_Equal(info,0)

        call opts_op%create(ppm_param_op_dcpse,info,order=2,c=1.4D0)
        or_fail("failed to initialize option object for operator")


        call Op%discretize_on(Part1,DCop,opts_op,info)
        Assert_Equal(info,0)
        Assert_True(associated(DCop))
        !call Op%discretize_on(Part1,PSEop,info,method="PSE")
        !Assert_Equal(info,0)
        !Assert_True(associated(PSEop))

        call DCop%compute(SField1,SField2,info)
        Assert_Equal(info,0)
        !testing output on a field that is not yet created nor discretized
        call DCop%compute(SField1,SField4,info)
        Assert_Equal(info,0)

        !call ppm_vtk_particles("output",Part1,info)

        Assert_True(inf_error(Part1,SField1,SField2,DCop).LT.tol_error)
        end_subroutine()
    end test

    test Gradient
        start_subroutine("test_gradient")

        tol_error = 1e-2

        call Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        call Part1%map_ghosts(info)
        Assert_Equal(info,0)

        call Part1%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))

        if (ndim .eq. 2) then
               degree =  (/1,0,   0,1/)
        else
               degree =  (/1,0,0, 0,1,0, 0,0,1/)
        endif
        coeffs = 1.0_mk

        call Op%create(ndim,coeffs,degree,info,name="Gradient")
        Assert_Equal(info,0)

        call opts_op%create(ppm_param_op_dcpse,info,order=2,&
        c=1.4D0,vector=.true.)
        or_fail("failed to initialize option object for operator")

        call Op%discretize_on(Part1,DCop,opts_op,info)
        Assert_Equal(info,0)
        Assert_True(associated(DCop))

        call DCop%compute(SField1,VFieldD,info)
        Assert_Equal(info,0)

        Assert_True(inf_error(Part1,SField1,VFieldD,DCop).LT.tol_error)

        end_subroutine()
    end test
!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f0_test(pos,ndim)

    real(mk)                              :: f0_test
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    if (ndim .eq. 2) then
        f0_test =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2))
    else
        f0_test =  sin(2._mk*pi*pos(1)) * cos(2._mk*pi*pos(2)) * &
            & sin(2._mk*pi*pos(3))
    endif

    return

end function f0_test


!-------------------------------------------------------------
! derivatives of the test function
!-------------------------------------------------------------
pure function df0_test(pos,order_deriv,ndim)

    real(mk)                                  :: df0_test
    integer                     , intent(in)  :: ndim
    real(mk), dimension(ppm_dim), intent(in)  :: pos
    integer,  dimension(ppm_dim), intent(in)  :: order_deriv

    select case (order_deriv(1))
    case (0)
        df0_test =    sin (2._mk*pi * pos(1))
    case (1)
        df0_test =    2._mk*pi*cos(2._mk*pi * pos(1))
    case (2)
        df0_test =  -(2._mk*pi)**2*sin(2._mk*pi * pos(1))
    case (3)
        df0_test =  -(2._mk*pi)**3*cos(2._mk*pi * pos(1))
    case (4)
        df0_test =   (2._mk*pi)**4*sin(2._mk*pi * pos(1))
    case default
        df0_test =  0._mk
    endselect

    select case (order_deriv(2))
    case (0)
        df0_test =   df0_test * cos (2._mk*pi * pos(2))
    case (1)
        df0_test =   df0_test * (-2._mk*pi)*sin(2._mk*pi * pos(2))
    case (2)
        df0_test =  df0_test * (-(2._mk*pi)**2)*cos(2._mk*pi * pos(2))
    case (3)
        df0_test =  df0_test * ( (2._mk*pi)**3)*sin(2._mk*pi * pos(2))
    case (4)
        df0_test =  df0_test * ( (2._mk*pi)**4)*cos(2._mk*pi * pos(2))
    case default
        df0_test =  0._mk
    endselect

    if (ndim .eq. 3 ) then
        select case (order_deriv(3))
        case (0)
            df0_test =   df0_test * sin (2._mk*pi * pos(3))
        case (1)
            df0_test =   df0_test * (2._mk*pi)*cos(2._mk*pi * pos(3))
        case (2)
            df0_test =  df0_test * (-(2._mk*pi)**2)*sin(2._mk*pi * pos(3))
        case (3)
            df0_test =  df0_test * (-(2._mk*pi)**3)*cos(2._mk*pi * pos(3))
        case (4)
            df0_test =  df0_test * ( (2._mk*pi)**4)*sin(2._mk*pi * pos(3))
        case default
            df0_test =  0._mk
        endselect
    endif

end function df0_test


!-------------------------------------------------------------
! Compute the infinity norm of the error
!-------------------------------------------------------------
function inf_error(Part1,Field1,Field2,Op)
    use ppm_module_data
    type(ppm_t_particles_d),TARGET   :: Part1
    type(ppm_t_field)                :: Field1,Field2
    class(ppm_t_operator_discr)       :: Op
    integer                          :: ip,nterms
    real(mk), dimension(:),  pointer :: wp_1 => NULL(), dwp_1 => NULL()
    real(mk), dimension(:,:),pointer :: xp=>NULL(), wp_2=>NULL(), dwp_2=>NULL()
    real(mk), dimension(:), allocatable :: err,exact
    real(mk)                         :: linf,inf_error,coeff
    logical                          :: input_vec,output_vec
    integer,dimension(:),allocatable :: degree,order
    integer,dimension(ndim)          :: dg
    character(len=100)               :: fname

    if (op%flags(ppm_ops_vector)) then
        call Part1%get(Field2,dwp_2,info,read_only=.true.)
        if (info.ne.0) write(*,*) "Failed te get Field2"
        output_vec = .true.
    else
        call Part1%get(Field2,dwp_1,info,read_only=.true.)
        if (info.ne.0) write(*,*) "Failed te get Field2"
        output_vec = .false.
    endif

    if (Field1%lda.EQ.1) THEN
        call Part1%get(Field1,wp_1,info,read_only=.true.)
        if (info.ne.0) write(*,*) "Failed te get Field1"
        input_vec = .false.
    else
        call Part1%get(Field1,wp_2,info,read_only=.true.)
        if (info.ne.0) write(*,*) "Failed te get Field1"
        input_vec = .true.
    endif

    nterms = op%op_ptr%nterms
    allocate(err(nterms),exact(nterms),degree(nterms*ndim),order(nterms))

    call Part1%get_xp(xp,info)
        if (info.ne.0) write(*,*) "Could not get xp"

    err = 0._mk
    linf = 0._mk
    do ip=1,Part1%Npart
        exact = 0._mk
        do i=1,nterms
            coeff = op%op_ptr%coeffs(i)
            dg = op%op_ptr%degree(1+(i-1)*ndim:i*ndim)
            if (output_vec) then
                exact(i) = coeff*df0_test(xp(1:ndim,ip),dg,ndim)
            else
                exact(1) = exact(1) + coeff*df0_test(xp(1:ndim,ip),dg,ndim)
            endif
        enddo
        if (output_vec) then
            do i=1,nterms
                err(i) = MAX(err(i),abs(dwp_2(i,ip) - exact(i)))
                !REMOVME
                !dwp_2(i,ip) = err(i)
            enddo
        else
            err(1) = MAX(err(1),abs(dwp_1(ip) - exact(1)))
        endif
        linf = MAX(linf,MAXVAL(abs(exact)))
    enddo

#ifdef __MPI
    call MPI_Allreduce(linf,linf,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
    call MPI_Allreduce(maxval(err),inf_error,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
#else
    inf_error = maxval(err)
#endif
    inf_error = inf_error/linf

    if (ppm_rank.eq.0) &
        write(*,*) '[',ppm_rank,']','Error is ',inf_error

    deallocate(err,exact,degree,order)
end function inf_error

end test_suite
