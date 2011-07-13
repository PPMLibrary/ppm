test_suite ppm_module_map_part



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
real(mk)                        :: tol,min_rcp,max_rcp
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: np = 100000
integer                         :: mp
integer                         :: newnp
real(mk),dimension(:,:),pointer :: xp
real(mk),dimension(:  ),pointer :: rcp
real(mk),dimension(:,:),pointer :: wp
real(mk),dimension(:  ),pointer :: min_phys,max_phys,h,p_h
real(mk),dimension(:  ),pointer :: len_phys
real(mk),dimension(:  ),pointer :: ghostlayer
integer, dimension(:  ),pointer :: ghostsize
integer                         :: i,j,k,sum1,sum2
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
integer, dimension(:  ),pointer :: nm
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok
real(mk)                         :: t0,t1,t2,t3

    init

        use ppm_module_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),ghostlayer(2*ndim),&
            &         nm(ndim),h(ndim),p_h(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        ghostlayer(1:2*ndim) = max_rcp
        bcdef(1:6) = ppm_param_bcdef_periodic
        
        nullify(xp,rcp,wp)

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

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize


    setup

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        allocate(randnb((1+ndim)*np),stat=info)
        do i=1,seedsize
            seed(i)=10+i*i*(rank+1)
        enddo
        call random_seed(put=seed)
        call random_number(randnb)
        
        allocate(xp(ndim,np),rcp(np),wp(pdim,np),stat=info)

    end setup
        

    teardown
        
        deallocate(xp,rcp,wp,stat=info)
        deallocate(seed,randnb)

    end teardown

    test global
        ! test global mapping

        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check

        !----------------
        ! create particles
        !----------------

        xp = 0.0_mk
        rcp = 0.0_mk

        !p_h = len_phys / real(npgrid,mk)
        !do j=1,npgrid
        !    do i=1,npgrid
        !        p_i = i + (j-1)*npgrid
        !        xp(1,p_i) = min_phys(1)+real(i-1,mk)*p_h(1)
        !        xp(2,p_i) = min_phys(2)+real(j-1,mk)*p_h(2)
        !        rcp(p_i) = min_rcp + (max_rcp-min_rcp)*randnb(p_i)
        !        do k=1,pdim
        !            wp(k,i) = rcp(i)*REAL(k,MK)
        !        enddo
        !    enddo
        !enddo
        do i=1,np
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        call ppm_map_part_global(topoid,xp,np,info)
        call ppm_map_part_push(rcp,np,info)
        call ppm_map_part_push(wp,pdim,np,info)
        call ppm_map_part_send(np,newnp,info)
        call ppm_map_part_pop(wp,pdim,np,newnp,info)
        call ppm_map_part_pop(rcp,np,newnp,info)
        call ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp

        call ppm_topo_check(topoid,xp,np,ok,info)

        assert_true(ok)

    end test


    test partial
        ! test partial mapping

        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check

        !----------------
        ! create particles
        !----------------

        xp = 0.0_mk
        rcp = 0.0_mk

        do i=1,np
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        call ppm_map_part_global(topoid,xp,np,info)
        call ppm_map_part_push(rcp,np,info)
        call ppm_map_part_push(wp,pdim,np,info)
        call ppm_map_part_send(np,newnp,info)
        call ppm_map_part_pop(wp,pdim,np,newnp,info)
        call ppm_map_part_pop(rcp,np,newnp,info)
        call ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp

        ! move all particles

        do i=1,np
            do j=1,ndim
                xp(j,i) = xp(j,i) + (len_phys(j)/REAL(nproc,mk))*&
                &         randnb((ndim+1)*i-(ndim-j))
            enddo
        enddo

        ! do local mapping
        call ppm_map_part_partial(topoid,xp,np,info)
        call ppm_map_part_push(rcp,np,info)
        call ppm_map_part_push(wp,pdim,np,info)
        call ppm_map_part_send(np,newnp,info)
        call ppm_map_part_pop(wp,pdim,np,newnp,info)
        call ppm_map_part_pop(rcp,np,newnp,info)
        call ppm_map_part_pop(xp,ndim,np,newnp,info)
        np=newnp

        call ppm_topo_check(topoid,xp,np,ok,info)
        assert_true(ok)
        
    end test

    test ghost_get
        ! test ghost get and how it treats pbc
        ! Ticket #32
        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_util_dbg
        integer                         :: npart = 1
        integer                         :: newnpart
        integer                         :: mpart
        real(mk),dimension(:,:),pointer :: p
        real(mk), parameter             :: gl = 0.1_mk
    
        allocate(p(ndim,npart))
        p(1,1) = 0.05_mk
        p(2,1) = 0.05_mk
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        call ppm_map_part_global(topoid,p,npart,info)
        call ppm_map_part_send(npart,newnpart,info)
        call ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart

        call ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,1,1,info,p,npart)

        call ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)
        
        call ppm_topo_check(topoid,p,npart,ok,info)
        print *, npart,mpart
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,2,1,info,p,npart,mpart)


    end test

    test getsymbc
        ! tests symmetric boundary conditions and ghost get
        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_util_dbg
        integer                         :: npart = 1
        integer                         :: newnpart
        integer                         :: mpart
        real(mk),dimension(:,:),pointer :: p
        real(mk), parameter             :: gl = 0.1_mk
    
        allocate(p(ndim,npart))
        p(1,1) = 0.05_mk
        p(2,1) = 0.05_mk
        bcdef(1:6) = ppm_param_bcdef_symmetry

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,p,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               gl,cost,info)

        call ppm_map_part_global(topoid,p,npart,info)
        call ppm_map_part_send(npart,newnpart,info)
        call ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart

        call ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        call ppm_dbg_print_d(topoid,gl,1,1,info,p,npart)

        call ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)
        
        call ppm_topo_check(topoid,p,npart,ok,info)
        print *, npart,mpart
        assert_true(ok)
        call ppm_dbg_print_d(topoid,gl,2,1,info,p,npart,mpart)


    end test


end test_suite
