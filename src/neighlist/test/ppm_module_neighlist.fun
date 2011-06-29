test_suite ppm_module_neighlist
use ppm_module_particles
use ppm_module_particles_typedef

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
integer                         :: np_global = 13000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.05_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys=>NULL(),max_phys=>NULL()
real(mk),dimension(:  ),pointer :: len_phys=>NULL()
integer                         :: i,j,k,isum1,isum2,ip,iq,ineigh
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost=>NULL()
real(mk)                        :: t0,t1,t2,t3
type(ppm_t_particles),pointer   :: Particles=>NULL()
type(ppm_t_particles),pointer   :: Particles2=>NULL()
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed

    init

        use ppm_module_typedef
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
        
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)

    end teardown

    test inhomogeneous_neighlists
        ! very preliminary for the time being. Just checking that it doesnt crash

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        call particles_initialize(Particles2,np_global,info,ppm_param_part_init_cartesian,topoid)
        allocate(disp(ndim,Particles%Npart))
        call random_number(disp)
        disp=0.15_mk*Particles%h_avg*disp
        call particles_move(Particles,disp,info)
        call random_number(disp)
        disp=0.15_mk*Particles%h_avg*disp
        call particles_move(Particles2,disp,info)
        call particles_apply_bc(Particles,topoid,info)
        call particles_apply_bc(Particles2,topoid,info)
        call particles_update_cutoff(Particles,Particles%h_avg*2.1_mk,info)

        call particles_update_cutoff(Particles2,Particles%h_avg*2.1_mk,info)
        call particles_mapping_global(Particles,topoid,info)
        call particles_mapping_global(Particles2,topoid,info)
        call particles_mapping_ghosts(Particles,topoid,info)
        call particles_mapping_ghosts(Particles2,topoid,info)

        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles2,topoid,info)
        Assert_Equal(info,0)


        Assert_False(associated(Particles2%Particles_cross))
        Assert_False(associated(Particles%Particles_cross))
        call particles_neighlists_xset(Particles2,Particles,topoid,info)
        Assert_Equal(info,0)
        Assert_True(associated(Particles2%Particles_cross))
        Particles%D_id = 129
        Assert_Equal(Particles2%Particles_cross%D_id,129)

    end test

end test_suite
