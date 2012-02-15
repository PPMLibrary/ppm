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

        use ppm_module_typedef
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
        use ppm_module_typedef
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
        real(mk)                        :: gl = 0.001_mk
        real(mk)                        :: h
        type(ppm_t_particles_d),pointer :: particles
        integer                         :: i,j
        CHARACTER(LEN=32)                       :: fname
        npart=10
        !CALL part_init(xp,npart,min_phys,max_phys,info,&
        !&    ppm_param_part_init_cartesian,0.5_mk)
        allocate(xp(ndim,npart),randnb(npart*ndim),stat=info)
        xp = 0.0_mk
        call random_number(randnb)
        do i=1,npart
          do j=1,ndim
            xp(j,i) = min_phys(j)+ len_phys(j)*randnb((ndim+1)*i-(ndim-j))
          enddo
        enddo
        h = 2.0_mk*(len_phys(1)/(sqrt(real(npart,mk))))
        bcdef(1:ndim) = ppm_param_bcdef_freespace
        
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        call ppm_map_part_global(topoid,xp,npart,info)
        call ppm_map_part_send(npart,newnpart,info)
        call ppm_map_part_pop(xp,ndim,npart,newnpart,info)
        npart=newnpart
        allocate(wp(npart),ra(npart))
        CALL RANDOM_NUMBER(wp)
        ra(1:npart) = rank
        call ppm_map_part_ghost_get(topoid,xp,ndim,npart,0,gl,info)
        call ppm_map_part_push(wp,npart,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(wp,npart,mpart,info)
        call ppm_map_part_pop(xp,ndim,npart,mpart,info)
        
        allocate(particles)
        particles%xp => xp
        particles%np = npart
        particles%nprop = 1
        particles%niprop = 1
        allocate(particles%prop(particles%nprop))
        allocate(particles%iprop(particles%niprop))
        particles%prop(1)%wp => wp
        particles%iprop(1)%wp => ra
        particles%prop(1)%name = 'bla'
        particles%iprop(1)%name = 'sub'
        fname = 'test'
        call ppm_vtk_particles(topoid,particles,fname,info)
        assert_equal(info,0)
    end test

    

end test_suite
