test_suite ppm_particles_dcops

use ppm_module_particles_typedef
use ppm_module_typedef
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
type(ppm_t_particles_d)         :: Pc
type(ppm_t_sop_d)               :: Pc_a
integer                         :: seedsize
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()
integer, dimension(:),allocatable:: seed,degree,order
real(mk),dimension(:,:),allocatable:: randn
real(mk),dimension(:),allocatable:: coeffs
integer, dimension(:),  pointer :: gi => NULL()
real(mk),dimension(:),  pointer :: wp_1r => NULL()
real(mk),dimension(:,:),pointer :: wp_2r => NULL()
real(mk)                        :: tol_error,err

    init

        use ppm_module_init
        use ppm_module_mktopo
        
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

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
        !initialize particles on a grid
        call Pc%initialize(np_global,info,topoid=topoid, &
            distrib=ppm_param_part_init_cartesian)

        call Pc%comp_global_index(info)

        call Pc%map(info,global=.true.,topoid=topoid)
        !Set up a small random displacement field 
        call Pc%create_prop(wp2_id,ppm_type_real,info,3,name='disp')
        call Pc%get(wp_2r,wp2_id)

        call Pc%get(gi,Pc%gi_id)
        FORALL (ip=1:Pc%Npart) wp_2r(1:ndim,ip) = randn(1:ndim,gi(ip))
        wp_2r = 0.5_mk * (wp_2r - 0.5_mk) * Pc%h_avg
        call Pc%set(gi,Pc%gi_id,read_only=.true.)

        !Perturb the  particles positions
        call Pc%move(wp_2r,info)
        call Pc%set(wp_2r,wp2_id)
        call Pc%apply_bc(info)
        call Pc%map(info)

        call Pc%create_prop(wp1_id,ppm_type_real,info,1,name='testf_sca')
        call Pc%create_prop(wp2_id,ppm_type_real,info,ndim,name='testf_vec')
        call Pc%get(wp_1r,wp1_id)
        call Pc%get(wp_2r,wp2_id)
        call Pc%get_xp(xp)
        DO ip=1,Pc%Npart
            wp_1r(ip) = f0_test(xp(1:ndim,ip),ndim)
            wp_2r(1:ndim,ip) = f0_test(xp(1:ndim,ip),ndim)
        ENDDO
        call Pc%set_xp(xp,read_only=.true.)
        call Pc%set(wp_1r,wp1_id)
        call Pc%set(wp_2r,wp2_id)
    end init


    finalize
        use ppm_module_finalize

        call Pc%destroy(info)
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

        tol_error = 2e-2

        call Pc%set_cutoff(3._mk * Pc%h_avg,info) 
        Assert_Equal(info,0)

        call Pc%map_ghosts(info)
        Assert_Equal(info,0)

        call Pc%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        if (ndim .eq. 2) then
               degree =  (/2,0,   0,2/)
        else 
               degree =  (/2,0,0, 0,2,0, 0,0,2/)
        endif
        coeffs = 1.0_mk
        order =  2

        call Pc%create_op(op_id,nterms,coeffs,degree,order,info,name='Laplacian')
        Assert_Equal(info,0)
        call Pc%comp_op(op_id,info)
        Assert_Equal(info,0)

        call Pc%map_ghosts(info)

        call Pc%apply_op(wp1_id,dwp1_id,op_id,info)
        Assert_Equal(info,0)
        Assert_True(inf_error(Pc,wp1_id,dwp1_id,op_id).LT.tol_error)

    end test

    test Gradient

        tol_error = 1e-2

        call Pc%set_cutoff(3._mk * Pc%h_avg,info) 
        Assert_Equal(info,0)

        call Pc%map_ghosts(info)
        Assert_Equal(info,0)

        call Pc%comp_neighlist(info)
        Assert_Equal(info,0)

        nterms=ndim
        allocate(degree(nterms*ndim),coeffs(nterms),order(nterms))
        call Pc%map_ghosts(info)

        if (ndim .eq. 2) then
               degree =  (/1,0,   0,1/)
        else 
               degree =  (/1,0,0, 0,1,0, 0,0,1/)
        endif
        coeffs = 1.0_mk
        order =  2
        call Pc%create_op(op_id,nterms,coeffs,degree,order,info,&
            name='gradient',vector=.true.)
        Assert_Equal(info,0)

        call Pc%comp_op(op_id,info)
        Assert_Equal(info,0)
        call Pc%apply_op(wp2_id,dwp1_id,op_id,info)
        Assert_Equal(info,0)

        Assert_True(inf_error(Pc,wp2_id,dwp1_id,op_id).LT.tol_error)

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
function inf_error(Pc,wp_id,dwp_id,op_id)
    use ppm_module_data
    type(ppm_t_particles_d)          :: Pc
    integer                          :: wp_id,dwp_id,op_id
    integer                          :: ip,nterms
    real(mk), dimension(:),  pointer :: wp_1, dwp_1
    real(mk), dimension(:,:),pointer :: xp, wp_2, dwp_2
    real(mk), dimension(:), allocatable :: err,exact
    real(mk)                         :: linf,inf_error,coeff
    logical                          :: input_vec,output_vec
    integer,dimension(:),allocatable :: degree,order
    integer,dimension(ndim)          :: dg
    character(len=100)               :: fname

    associate (op => Pc%ops%vec(op_id)%t)
        if (op%flags(ppm_ops_vector)) then
            call Pc%get(dwp_2,dwp_id)
            output_vec = .true.
        else
            call Pc%get(dwp_1,dwp_id)
            output_vec = .false.
        endif
        if (Pc%props%vec(wp_id)%t%lda.EQ.1) then
            call Pc%get(wp_1,wp_id)
            output_vec = .false.
        else
            call Pc%get(wp_2,wp_id)
            input_vec = .true.
        endif

        nterms = op%desc%nterms
        allocate(err(nterms),exact(nterms),degree(nterms*ndim),order(nterms))

        call Pc%get_xp(xp)
                
        err = 0._mk
        linf = 0._mk
        do ip=1,Pc%Npart
            exact = 0._mk
            do i=1,nterms
                coeff = op%desc%coeffs(i)
                dg = op%desc%degree(1+(i-1)*ndim:i*ndim)
                if (output_vec) then
                    exact(i) = coeff*df0_test(xp(1:ndim,ip),dg,ndim)
                else
                    exact(1) = exact(1) + coeff*df0_test(xp(1:ndim,ip),dg,ndim)
                endif
            enddo
            if (output_vec) then
                do i=1,op%desc%nterms
                    err(i) = MAX(err(i),abs(dwp_2(i,ip) - exact(i)))
                    !REMOVME
                    dwp_2(i,ip) = err(i)
                enddo
            else
                err(1) = MAX(err(1),abs(dwp_1(ip) - exact(1)))
            endif
            linf = MAX(linf,MAXVAL(abs(exact)))
        enddo

    end associate

!    write(fname,'(A,I0)') 'test_grad',ppm_nproc
!    CALL ppm_vtk_particle_cloud(TRIM(ADJUSTL(fname)),Pc,info)

#ifdef __MPI
    call MPI_Allreduce(linf,linf,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
    call MPI_Allreduce(maxval(err),inf_error,1,ppm_mpi_kind,MPI_MAX,ppm_comm,info)
#elif
    inf_error = maxval(err)
#endif
    inf_error = inf_error/linf

!    if (ppm_rank.eq.0) &
!        write(*,*) '[',ppm_rank,']','Error is ',inf_error

    deallocate(err,exact,degree,order)
end function inf_error

end test_suite
