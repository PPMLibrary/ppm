test_suite ppm_module_particles


#include "../../ppm_define.h"

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: tol=epsilon(1._mk)*100
real(mk),parameter              :: pi = 3.1415926535897931_mk
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim=3
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
integer                         :: topoid,nneigh_theo
integer                         :: np_global = 100000
integer                         :: npart_g
real(mk),parameter              :: cutoff = 0.15_mk
real(mk),dimension(:,:),pointer :: xp=>NULL(),disp=>NULL()
real(mk),dimension(:  ),pointer :: min_phys,max_phys
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: rcp,wp
integer                         :: i,j,k,isum1,isum2,ip,wp_id
real(mk)                        :: rsum1,rsum2
integer                         :: nstep
real(mk),dimension(:),pointer   :: delta
integer,dimension(3)            :: ldc
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
character(len=ppm_char)         :: dirname
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
            &         delta(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
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

    test initialize_cart
        ! test initialization of particles on a grid

        use ppm_module_typedef

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        Assert_Equal(Particles%active_topoid,topoid)
        Assert_False(Particles%adaptive)
        Assert_False(Particles%ontopology)
#ifdef __MPI
        call MPI_Allreduce(Particles%Npart,npart_g,1,MPI_INTEGER,MPI_SUM,comm,info)
#else
        npart_g = Particles%Npart
#endif
        Assert_Equal(npart_g,np_global)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)
        Assert_Equal(info,0)

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,&
            topoid,cutoff=cutoff)
        Assert_Equal(info,0)
        Assert_Equal(Particles%cutoff,cutoff)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)
        Assert_Equal(info,0)

    end test

    test initialize_rand
        ! test initialization of particles randomly

        use ppm_module_typedef

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_random,topoid)
        Assert_Equal(info,0)
        Assert_False(Particles%adaptive)
#ifdef __MPI
        call MPI_Allreduce(Particles%Npart,npart_g,1,MPI_INTEGER,MPI_SUM,comm,info)
#else
        npart_g = Particles%Npart
#endif
        Assert_Equal(npart_g,np_global)
        call ppm_alloc_particles(Particles,np_global,ppm_param_dealloc,info)
        Assert_Equal(info,0)

    end test

    test global1
        ! test global mapping of particles with cartesian particles

        use ppm_module_typedef
        use ppm_module_topo_check

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        Assert_False(Particles%ontopology)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%active_topoid,topoid)
        call ppm_topo_check(Particles%active_topoid,Particles%xp,Particles%Npart,ok,info)
        Assert_True(ok)

    end test

    test global2
        ! test global mapping of particles with cartesian particles plus displacement

        use ppm_module_typedef
        use ppm_module_topo_check

        ok=.false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        allocate(disp(ndim,Particles%Npart),stat=info)
        Assert_Equal(info,0)
        call random_number(disp)
        Assert_True(Particles%cutoff>0._MK)
        
        xp => get_xp(Particles)
        Assert_True(associated(xp))
        FORALL(ip=1:Particles%Npart) &
            xp(1:ndim,ip) = xp(1:ndim,ip) + (disp(1:ndim,ip)-0.5_mk)*Particles%cutoff
        xp => set_xp(Particles)
        call particles_apply_bc(Particles,topoid,info)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%active_topoid,topoid)
        call ppm_topo_check(Particles%active_topoid,Particles%xp,Particles%Npart,ok,info)
        Assert_True(ok)

    end test

    test global3
        ! test global mapping of particles with random particles

        use ppm_module_typedef
        use ppm_module_topo_check

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        Assert_False(Particles%ontopology)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        Assert_Equal(Particles%active_topoid,topoid)
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
        
        xp => get_xp(Particles)
        Assert_True(associated(xp))
        FORALL(ip=1:Particles%Npart) &
            xp(1:ndim,ip) = xp(1:ndim,ip) + (disp(1:ndim,ip)-0.5_mk)*Particles%cutoff
        xp => set_xp(Particles,read_only=.true.)
        Assert_False(associated(xp))
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

        wp_id = 0
        call particles_allocate_wps(Particles,wp_id,info,zero=.true.)
        Assert_Equal(info,0)
        wp => get_wps(Particles,wp_id)
        xp => get_xp(Particles)
        Assert_True(associated(wp))
        Assert_True(associated(xp))
        Assert_Equal(maxval(wp(1:Particles%Npart)),0._mk)
        forall (ip=1:Particles%Npart) wp(ip) = cos(product(xp(1:ndim,ip)))
        rsum1 = SUM(wp(1:Particles%Npart))
        isum1 = Particles%Npart 
        wp => set_wps(Particles,wp_id)
        xp => set_xp(Particles,read_only=.true.)
        Assert_False(associated(wp))
        Assert_False(associated(xp))

        nstep = ceiling(minval(len_phys)/(0.8_mk*cutoff))
        delta = len_phys/(1._mk*nstep)
        do i=1,nstep
            ldc(1) = ndim
            ldc(2) = Particles%Npart
            call ppm_alloc(disp,ldc,ppm_param_alloc_grow,info)
            Assert_Equal(info,0)
            forall (ip=1:Particles%Npart) disp(:,ip) = delta
            call particles_move(Particles,disp,info)
            Assert_Equal(info,0)
            call particles_apply_bc(Particles,topoid,info)
            call particles_mapping_partial(Particles,topoid,info)
            Assert_Equal(info,0)
        enddo
        wp => get_wps(Particles,wp_id)
        xp => get_xp(Particles)
        Assert_True(associated(wp))
        Assert_True(associated(xp))
        forall (ip=1:Particles%Npart) wp(ip) = cos(product(xp(1:ndim,ip)))
        rsum2 = SUM(wp(1:Particles%Npart))
        isum2 = Particles%Npart 
        wp => set_wps(Particles,wp_id)
        xp => set_xp(Particles,read_only=.true.)
        Assert_Equal_Within(rsum1/rsum2,1._mk,tol)
        Assert_Equal(isum1,isum2)


        deallocate(disp)

    end test

    test neighlists
        ! test construction of neighbour lists for particles 

        use ppm_module_typedef
        use ppm_module_topo_check

        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        Assert_Equal(info,0)
        call particles_mapping_global(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_mapping_ghosts(Particles,topoid,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)

        if (ndim.eq.2) then
            nneigh_theo = 12
        else
            nneigh_theo = 32
        endif
        Assert_Equal(Particles%nneighmin,nneigh_theo)

        call particles_update_cutoff(Particles,1.5_mk*Particles%h_avg,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)
        if (ndim.eq.2) then
            nneigh_theo = 8
        else
            nneigh_theo = 18
        endif
        Assert_Equal(Particles%nneighmin,nneigh_theo)

        call particles_update_cutoff(Particles,1.1_mk*Particles%h_avg,info)
        Assert_Equal(info,0)
        call particles_neighlists(Particles,topoid,info)
        Assert_Equal(info,0)
        if (ndim.eq.2) then
            nneigh_theo = 4
        else
            nneigh_theo = 6
        endif
        Assert_Equal(Particles%nneighmin,nneigh_theo)

        call particles_allocate_wps(Particles,Particles%rcp_id,info,with_ghosts=.true.)
        Assert_Equal(info,0)
        rcp => get_wps(Particles,Particles%rcp_id,with_ghosts=.true.)
        xp => get_xp(Particles,with_ghosts=.true.)
        DO ip=1,Particles%Mpart
            rcp(ip) = 2.01_MK*Particles%h_avg*(abs(cos(xp(1,ip)+sin(xp(2,ip))))+0.5_MK)
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
        if (ndim.eq.2) then
            nneigh_theo = 4
        else
            nneigh_theo = 6
        endif

!        write(121+rank,*) Particles%xp(:,MINLOC(Particles%nvlist(1:Particles%Npart),1))
!        ip = MINLOC(Particles%nvlist(1:Particles%Npart),1)
!        write(*,*) 'ip = ',ip, 'Min = ',MINVAL(Particles%nvlist(1:Particles%Npart))
!        do i = 1, MINVAL(Particles%nvlist(1:Particles%Npart))
!            write(121+rank,*) Particles%xp(:,Particles%vlist(i,ip))
!        enddo

        Assert_Equal(Particles%nneighmin,nneigh_theo)
        if (ndim.eq.2) then
            nneigh_theo = 28
        else
            nneigh_theo = 112
        endif
        Assert_Equal(Particles%nneighmax,nneigh_theo)

    end test

    test io
        ! test i/o routines

        use ppm_module_typedef

        ok = .false.
        call particles_initialize(Particles,np_global,info,ppm_param_part_init_cartesian,topoid)
        write(dirname,*) './'
        call particles_io_xyz(Particles,0,dirname,info)
        Assert_Equal(info,0)
    end test


end test_suite
