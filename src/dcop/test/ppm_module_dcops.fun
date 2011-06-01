test_suite ppm_module_dcops
use ppm_module_particles
use ppm_module_particles_typedef
#include "../../ppm_define.h"

#ifdef __MPI
    INCLUDE "mpif.h"
#endif


integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid,nneigh_theo
integer                         :: np_global = 10000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: rcp,wp
integer                         :: i,j,k,isum1,isum2,ip,wp_id,eta_id,eta2_id
real(mk)                        :: rsum1,rsum2
integer                         :: nstep
real(mk),dimension(:),pointer   :: delta
integer,dimension(3)            :: ldc
integer,dimension(:),allocatable:: degree,degree2,order,order2
real(mk),dimension(:),allocatable:: coeffs,coeffs2
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
character(len=ppm_char)         :: dirname
integer                         :: isymm = 0
logical                         :: lsymm = .false.,ok
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles),pointer   :: Particles=>NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()

    init

        use ppm_module_typedef
        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         delta(ndim),stat=info)
        
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
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,bcdef,cutoff,cost,info)
    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys,delta)

    end finalize


    setup


    end setup
        

    teardown
        
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)

    end teardown

    test allocate_operator
        ! test data structure

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_mapping_global(Particles,topoid,info)

        allocate(degree(3*ndim),coeffs(3),order(3),degree2(7*ndim),coeffs2(7),order2(7))
        if (ndim .eq. 2) then
               degree =  (/1,0,    1,1,     0,1  /)
               degree2=  (/1,8,   0,1,   3,3,   1,1,   0,7,   0,2,   3,3/)
        else 
               degree =  (/1,0,0,  1,1,1,  0,0,1 /)
               degree2=  (/1,2,9, 1,1,8, 3,3,1, 1,1,1, 3,0,7, 2,0,2, 3,3,3/)
        endif
        coeffs = (/2.0_mk, 1.0_mk, -3.2_mk/)
        order =  (/2,       1,       1  /)
        coeffs2 = (/0.1_mk, -1.0_mk, -3.8_mk, -3.3_mk, 0.001_mk, 10._mk, 4._mk/)
        order2=  (/2,       2,       2,        2,      1,        2,      1/)

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test")
        Assert_Equal(info,0)
        Assert_False(eta_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta_id)%degree-degree)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta_id)%coeffs-coeffs)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,1)

        eta2_id = 0
        call particles_dcop_define(Particles,eta2_id,coeffs2,degree2,order,7,info,name="test2")
        Assert_False(info.eq.0)
        call particles_dcop_define(Particles,eta2_id,coeffs2,degree2,order2,7,info,name="test2")
        Assert_Equal(info,0)
        Assert_False(eta2_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta2_id)%degree-degree2)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta2_id)%coeffs-coeffs2)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,2)

        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%ops%nb_ops,1)
        Assert_Equal(Particles%ops%max_opsid,2)
        
        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test1")
        Assert_Equal(info,0)
        Assert_False(eta_id.le.0)
        Assert_Equal(sum(abs(Particles%ops%desc(eta_id)%degree-degree)),0)
        Assert_Equal_Within(sum(abs(Particles%ops%desc(eta_id)%coeffs-coeffs)),0,1e-5)
        Assert_Equal(Particles%ops%nb_ops,2)
        Assert_Equal(Particles%ops%max_opsid,2)

        call particles_dcop_free(Particles,eta_id,info)
        Assert_Equal(info,0)
        call particles_dcop_free(Particles,eta2_id,info)
        Assert_Equal(info,0)

        deallocate(degree,coeffs,order,degree2,coeffs2,order2)
    end test

    test compute_operator
        ! test if we can compute the dc operators

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_random,topoid)
        Assert_Equal(info,0)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)

        allocate(degree(3*ndim),coeffs(3),order(3),degree2(7*ndim),coeffs2(7),order2(7))
        if (ndim .eq. 2) then
               degree =  (/0,0,    1,1,     0,1  /)
               degree2=  (/0,0,   0,1,   2,1,   2,2,   0,2,   1,2,   1,1/)
        else 
               degree =  (/0,0,0,  1,1,1,  2,0,1 /)
               degree2=  (/0,0,0, 2,1,0, 0,1,1, 2,1,1, 1,0,1, 2,2,2, 1,1,1/)
        endif
        coeffs = (/2.0_mk, 1.0_mk, -3.2_mk/)
        order =  (/2,       1,       1  /)
        coeffs2 = (/0.1_mk, -1.0_mk, -3.8_mk, -3.3_mk, 0.001_mk, 10._mk, 4._mk/)
        order2=  (/2,       2,       2,        2,      1,        2,      1/)

        eta_id = 0
        call particles_dcop_define(Particles,eta_id,coeffs,degree,order,3,info,name="test")
        Assert_Equal(info,0)
        eta2_id = 0
        call particles_dcop_define(Particles,eta2_id,coeffs2,degree2,order2,7,info,name="test2")
        Assert_Equal(info,0)

        call particles_dcop_compute(Particles,eta_id,info)
        Assert_False(info.eq.0)

        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta_id,info)
        Assert_Equal(info,0)
        call particles_dcop_compute(Particles,eta2_id,info)
        Assert_Equal(info,0)

        deallocate(degree,coeffs,order,degree2,coeffs2,order2)
    end test

pure function f0_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk)                                 :: f0_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                                 :: radius,eps

    centre = 0.5_mk
    centre(2) = 0.75_mk
    radius=0.15_mk
    eps = 0.01_mk

    f0_fun = tanh((sqrt(sum((pos(1:ppm_dim)-centre)**2)) - radius)/eps)

end function f0_fun

pure function f0_grad_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk), dimension(ppm_dim)             :: f0_grad_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                                 :: radius,eps,f0,d

    centre = 0.5_mk
    centre(2) = 0.75_mk
    radius=0.15_mk
    eps = 0.01_mk

    d = sqrt(sum((pos(1:ppm_dim)-centre)**2))
    f0 = tanh((d - radius)/eps)
    f0_grad_fun = (1._mk - f0**2) * (pos(1:ppm_dim)-centre)/(eps*d)

end function f0_grad_fun

pure function level0_fun(pos)

    use ppm_module_data, ONLY: ppm_dim
    real(mk)                              :: level0_fun
    real(mk), dimension(ppm_dim), intent(in) :: pos
    real(mk), dimension(ppm_dim)             :: centre
    real(mk)                              :: radius

    centre = 0.5_mk
    centre(2) = 0.75_mk
    radius=0.15_mk

    level0_fun = sqrt(sum((pos(1:ppm_dim)-centre)**2)) - radius

end function level0_fun

end test_suite
