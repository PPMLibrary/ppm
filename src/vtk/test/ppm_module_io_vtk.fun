test_suite ppm_module_io_vtk



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
integer,parameter               :: ndim=3
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank,nproc
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
real(mk),dimension(:  ),pointer :: ghostlayer => NULL()
integer                         :: i,j,k,sum1,sum2
integer, dimension(ndim*2)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3
real(mk)                         :: eps

    init

        use ppm_module_data
        use ppm_module_topo_typedef
        use ppm_module_init

        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostlayer(2*ndim),stat=info)

        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:ndim) = ppm_param_bcdef_periodic

        eps = epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        call ppm_init(ndim,mk,tolexp,0,debug,info,99)

    end init


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)


    end setup


    teardown

        deallocate(seed)

    end teardown

    test vtkparticles
        use ppm_module_topo_typedef
        use ppm_module_interfaces
        use ppm_module_particles_typedef
        use ppm_module_field_typedef
        use ppm_module_mktopo
        use ppm_module_map
        use ppm_module_topo_check
        use ppm_module_util_dbg
        use ppm_module_test

        integer                         :: topoid
        integer                         :: npart
        integer                         :: mpart
        integer                         :: newnpart
        integer                         :: oldip,ip = -1
        real(mk),dimension(:,:),pointer :: xp => NULL()
        real(mk),dimension(:),pointer   :: wp => NULL()
        integer ,dimension(:),pointer   :: ra => NULL()
        real(mk)                        :: cutoff = 0.001_mk
        real(mk)                        :: h
        type(ppm_t_particles_d)         :: particles
        type(ppm_t_field)               :: Field1
        class(ppm_t_discr_data),pointer :: Prop1 => NULL()
        integer                         :: i,j
        CHARACTER(LEN=32)               :: fname

        start_subroutine("vtk")

        npart=100
        bcdef(1:ndim) = ppm_param_bcdef_freespace
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0
        call ppm_mktopo(topoid,decomp,assig,min_phys,&
            max_phys,bcdef,cutoff,cost,info)
            assert_equal(info,0)

        call Field1%create(ndim,info,name="F_vec") !vector field
            assert_equal(info,0)

        !initialize particles on a grid
        call particles%initialize(npart,info,topoid=topoid, &
            distrib=ppm_param_part_init_cartesian)
            assert_equal(info,0)

        !define a property on the particle set by discretizing a Field on it
        call Field1%discretize_on(particles,info)
            assert_equal(info,0)

        !define a property directly
        call particles%create_prop(info,discr_data=Prop1,dtype=ppm_type_real,&
            name='test_r',zero=.true.)
            assert_equal(info,0)

        foreach p in particles(particles) with positions(x) sca_fields(F2=Prop1) vec_fields(F1=Field1)
            F1_p(1:ndim) = 10._mk * x_p(1:ndim)
            F2_p         = 42._mk
        end foreach

        fname = 'test'
        call ppm_vtk_particles(fname,particles,info)
            assert_equal(info,0)

        end_subroutine()
    end test



end test_suite
