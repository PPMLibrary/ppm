test_suite ppm_module_particles


#include "../../ppm_define.h"

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=2
integer,parameter               :: pdim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: np_global = 1000 !100000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.1_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: ghostlayer
real(mk),dimension(:  ),pointer :: rcp
integer                         :: i,j,k,sum1,sum2,ip
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
integer                         :: isymm = 0
logical                         :: lsymm = .false.,ok
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles),pointer   :: Particles
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
integer, dimension(:),pointer   :: nvlist=>NULL()
integer, dimension(:,:),pointer :: vlist=>NULL()

    init

        use ppm_module_typedef
        use ppm_module_init
        use ppm_module_mktopo
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostlayer(2*ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostlayer(1:2*ndim) = cutoff
        bcdef(1:6) = ppm_param_bcdef_periodic
        
        nullify(xp)

#ifdef __MPI
        comm = mpi_comm_world
        call mpi_comm_rank(comm,rank,info)
        call mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1._MK)))+10
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

        deallocate(min_phys,max_phys,len_phys)

    end finalize


    setup


    end setup
        

    teardown
        
        deallocate(xp,stat=info)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)

    end teardown

    test initialize_cart
        ! test initialization of particles on a grid

        use ppm_module_typedef
        use ppm_module_topo_check

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        Assert_False(Particles%adaptive)
        call MPI_Allreduce(Particles%Npart,npart_g,1,MPI_INTEGER,MPI_SUM,comm,info)
        Assert_Equal(npart_g,np_global)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)
        Assert_Equal(info,0)

    end test

    test initialize_rand
        ! test initialization of particles randomly

        use ppm_module_typedef
        use ppm_module_topo_check

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_random,topoid)
        Assert_Equal(info,0)
        Assert_False(Particles%adaptive)
        call MPI_Allreduce(Particles%Npart,npart_g,1,MPI_INTEGER,MPI_SUM,comm,info)
        Assert_Equal(npart_g,np_global)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)
        Assert_Equal(info,0)

    end test

    test global
        ! test global mapping of particles

        use ppm_module_typedef
        use ppm_module_topo_check

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_mapping_global(Particles,topoid,info)
        call ppm_topo_check(Particles%active_topoid,Particles%xp,Particles%Npart,ok,info)
        Assert_True(ok)

    end test

    test mappings
        ! test partial and ghost mappings for particles

        use ppm_module_typedef
        use ppm_module_topo_check

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_random,topoid)
        Assert_Equal(info,0)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        call ppm_topo_check(Particles%active_topoid,Particles%xp,Particles%Npart,ok,info)
        Assert_True(ok)
        allocate(disp(ndim,Particles%Npart),stat=info)
        Assert_Equal(info,0)
        call random_number(disp)
        Assert_True(Particles%cutoff>0._MK)
        Particles%xp = Particles%xp + (disp-0.5_mk)*Particles%cutoff
        call particles_updated_positions(Particles,info)
        Assert_Equal(info,0)

!try to do a partial mapping while particles are outside the domain
        call particles_mapping_partial(Particles,topoid,info)
        Assert_FALSE(info.eq.0)

        call particles_apply_bc(Particles,topoid,info)
        Assert_Equal(info,0)

!try to do a ghost mapping before the partial mapping
        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_False(info.eq.0)

        call particles_mapping_partial(Particles,topoid,info)
        Assert_Equal(info,0)
        call ppm_topo_check(Particles%active_topoid,Particles%xp,Particles%Npart,ok,info)
        Assert_True(ok)
        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_Equal(info,0)
        deallocate(disp)

    end test

    test neighlists
        ! test construction of neighbour lists for particles 

        use ppm_module_typedef
        use ppm_module_topo_check

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_mapping_global(Particles,topoid,info)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%nneighmin,12)
        call particles_updated_cutoff(Particles,info,cutoff_new=1.5_mk*Particles%h_avg)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(Particles%nneighmin,8)
        call particles_updated_cutoff(Particles,info,cutoff_new=1.1_mk*Particles%h_avg)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(Particles%nneighmin,4)
        call particles_allocate_wps(Particles,Particles%rcp_id,info,with_ghosts=.true.)
        Assert_Equal(info,0)
        rcp => get_wps(Particles,Particles%rcp_id,with_ghosts=.true.)
        xp => get_xp(Particles,with_ghosts=.true.)
        DO ip=1,Particles%Mpart
            rcp(ip) = 2._MK*Particles%h_avg*(abs(cos(xp(1,ip)+sin(xp(2,ip))))+0.5_MK)
        ENDDO
        rcp => NULL()
        xp  => NULL()
        call particles_updated_cutoff(Particles,info)
        Assert_Equal(info,0)
!Try to compute neighbour lists before updating the ghosts
        call particles_neighlists(Particles,topoid,info)
        Assert_False(info.eq.0)

        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles,topoid,info)
        Assert_True(Particles%adaptive)
        Assert_Equal(Particles%nneighmin,4)
        Assert_Equal(Particles%nneighmax,24)
        
    end test

end test_suite
