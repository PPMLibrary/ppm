test_suite ppm_module_part_modify



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
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:  ),pointer :: rcp => NULL()
real(mk),dimension(:,:),pointer :: wp => NULL()
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
real(mk),dimension(:  ),pointer :: p_h => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
real(mk),dimension(:  ),pointer :: ghostlayer => NULL()
integer, dimension(:  ),pointer :: ghostsize => NULL()
integer                         :: i,j,k,sum1,sum2
integer                         :: p_i
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
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
        tol = epsilon(1.0_mk)
        tolexp = int(log10(epsilon(1.0_mk)))
        
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

    test ghost_get_add
        ! test ghost get followed by adding particles

        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        real(mk), parameter             :: gl = 0.1_mk
        real(mk),dimension(:,:),pointer :: xpn => NULL()
        real(mk),dimension(:,:),pointer :: wpn => NULL()
        real(mk),dimension(:),pointer   :: rcpn => NULL()
        integer                         :: np_added=0,mp_new=0,np_new=0
        integer                         :: npart = 60000,i

        !----------------
        ! create particles
        !----------------

        xp = 0.0_mk
        rcp = 0.0_mk
        min_rcp = 0.01_mk
        max_rcp = 0.2_mk 

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
        do i=1,npart
            do j=1,ndim
                xp(j,i) = min_phys(j)+&
                len_phys(j)*randnb((ndim+1)*i-(ndim-j))
            enddo
            rcp(i) = min_rcp + (max_rcp-min_rcp)*randnb((ndim+1)*i-ndim)
            do j=1,pdim
                wp(j,i) = rcp(i)*REAL(j,MK)
            enddo
        enddo
        xp(1,1) = 0.02_mk
        xp(2,1) = 0.2_mk

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        !decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        call ppm_mktopo(topoid,xp,npart,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_rcp,cost,info)

        call ppm_map_part_global(topoid,xp,npart,info)
        call ppm_map_part_push(rcp,npart,info)
        call ppm_map_part_push(wp,pdim,npart,info)
        call ppm_map_part_send(npart,newnp,info)
        call ppm_map_part_pop(wp,pdim,npart,newnp,info)
        call ppm_map_part_pop(rcp,npart,newnp,info)
        call ppm_map_part_pop(xp,ndim,npart,newnp,info)
        npart=newnp

        call ppm_topo_check(topoid,xp,npart,ok,info)

        assert_true(ok)

        call ppm_map_part_ghost_get(topoid,xp,ndim,npart,0,gl,info)
        call ppm_map_part_push(rcp,npart,info)
        call ppm_map_part_push(wp,pdim,npart,info)
        call ppm_map_part_send(npart,mp,info)
        call ppm_map_part_pop(wp,pdim,npart,mp,info)
        call ppm_map_part_pop(rcp,npart,mp,info)
        call ppm_map_part_pop(xp,ndim,npart,mp,info)
        
        call ppm_topo_check(topoid,xp,npart,ok,info)
        assert_true(ok)

        np_added = 1
        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added))
        xpn(1,1)= 0.05_mk
        xpn(2,1)= 0.05_mk
        wpn(1,1)= 7._mk
        wpn(2,1)= 17._mk
        rcpn(1)= 0.33_mk

        !add one particle in the corner (it should generate 3 ghosts)
        call ppm_part_modify_add(topoid,xp,npart,mp,xpn,np_added,np_new,0,gl,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(wp,pdim,npart,mp,wpn,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(rcp,npart,mp,rcpn,info)
        assert_true(info.eq.0)
        call ppm_part_modify_send(info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(rcp,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(wp,pdim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(xp,ndim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)

        mp = mp_new

        !modify properties
        do i=1,np_new
          rcp(i) = 1._mk*i
          wp(1,i) = 1._mk*i
          wp(2,i) = 2._mk*i
        enddo
        rcp(1)=1.7_mk
        rcp(np_new)=27._mk

        !update ghosts without calling ghost_get
        call ppm_map_part_push(rcp,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_push(wp,pdim,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_send(np_new,mp,info)
        assert_true(info.eq.0)
        assert_true(np_new.eq.npart+1)
        call ppm_map_part_pop(wp,pdim,np_new,mp,info)
        assert_true(info.eq.0)
        call ppm_map_part_pop(rcp,np_new,mp,info)
        assert_true(info.eq.0)
        assert_true(rcp(np_new+1).eq.1.7_mk)
        assert_true(rcp(mp-2).eq.27._mk)
        assert_true(rcp(mp-1).eq.27._mk)
        assert_true(rcp(mp).eq.27._mk)

        deallocate(xpn,wpn,rcpn)

        npart = np_new

        !lets add some more and do it again
        np_added = 2011
        allocate(xpn(ndim,np_added),wpn(pdim,np_added),rcpn(np_added))
        call random_number(xpn)
        call random_number(wpn)
        call random_number(rcpn)

        call ppm_part_modify_add(topoid,xp,npart,mp,xpn,np_added,np_new,0,gl,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(wp,pdim,npart,mp,wpn,info)
        assert_true(info.eq.0)
        call ppm_part_modify_push(rcp,npart,mp,rcpn,info)
        assert_true(info.eq.0)
        call ppm_part_modify_send(info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(rcp,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(wp,pdim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)
        call ppm_part_modify_pop(xp,ndim,npart,mp,np_new,mp_new,info)
        assert_true(info.eq.0)

        mp = mp_new

        !modify properties
        do i=1,np_new
          rcp(i) = 3.14_mk
          wp(1,i) = 7.14_mk
          wp(2,i) = 8.14_mk
        enddo

        !update ghosts without calling ghost_get
        call ppm_map_part_push(rcp,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_push(wp,pdim,np_new,info)
        assert_true(info.eq.0)
        call ppm_map_part_send(np_new,mp,info)
        assert_true(info.eq.0)
        assert_true(np_new.eq.npart+np_added)
        call ppm_map_part_pop(wp,pdim,np_new,mp,info)
        assert_true(info.eq.0)
        call ppm_map_part_pop(rcp,np_new,mp,info)
        assert_true(info.eq.0)
        assert_true(rcp(np_new+1).eq.3.14_mk)
        assert_true(rcp(mp).eq.3.14_mk)

        deallocate(xpn,wpn,rcpn)

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
        real(mk),dimension(:,:),pointer :: p => NULL()
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
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,2,1,info,p,npart,mpart)


    end test

    test pbc_and_map_load
        ! test ghost get map load and periodic BC

        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_util_dbg
        integer                         :: npart = 8
        integer                         :: newnpart
        integer                         :: mpart
        real(mk),dimension(:,:),pointer :: p => NULL()
        real(mk),dimension(:)  ,pointer :: w => NULL()
        real(mk),dimension(2)           :: check
        real(mk), parameter             :: gl = 0.1_mk
    
        allocate(p(ndim,npart),w(npart))
        p(1,1) = 0.05_mk  ! left
        p(2,1) = 0.5_mk
        p(1,2) = 0.95_mk  ! right
        p(2,2) = 0.5_mk
        p(1,3) = 0.5_mk   ! bottom
        p(2,3) = 0.05_mk
        p(1,4) = 0.5_mk   ! top
        p(2,4) = 0.95_mk
        
        p(1,5) = 0.05_mk  ! left-bottom
        p(2,5) = 0.05_mk
        p(1,6) = 0.05_mk  ! left-top
        p(2,6) = 0.95_mk
        p(1,7) = 0.95_mk  ! right-bottom
        p(2,7) = 0.05_mk
        p(1,8) = 0.95_mk  ! right-top
        p(2,8) = 0.95_mk

        w(:) = 1.0_mk

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
        call ppm_map_part_push(w,npart,info)
        call ppm_map_part_send(npart,newnpart,info)
        call ppm_map_part_pop(w,npart,newnpart,info)
        call ppm_map_part_pop(p,ndim,npart,newnpart,info)
        npart=newnpart
        call ppm_topo_check(topoid,p,npart,ok,info)

        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,1,1,info,p,npart)

        call ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        call ppm_map_part_push(w,npart,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(w,npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)
        
        call ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,2,1,info,p,npart,mpart)
        call ppm_map_part_store(info)

        call ppm_map_part_ghost_put(topoid,info)
        call ppm_map_part_push(w,npart,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_ghost_pop(w,1,npart,mpart,info)

        call ppm_map_part_load(info)
        call ppm_map_part_push(p,ndim,npart,info,.TRUE.)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)

        call ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,3,1,info,p,npart,mpart)

        assert_equal(mpart-npart,4+3*4) ! check number of ghosts

        ! now go through all particles and try to find their ghosts  

        check(1) = 1.05_mk  ! left
        check(2) = 0.5_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        
        check(1) = -0.05_mk  ! right
        check(2) =  0.5_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = 0.5_mk   ! bottom
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  0.5_mk   ! top
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
       

        check(1) = 1.05_mk  ! left-bottom
        check(2) = 0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = 0.05_mk
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = 1.05_mk
        check(2) = 1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) =  1.05_mk  ! left-top
        check(2) =  0.95_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  1.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = -0.05_mk  ! right-bottom
        check(2) =  0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.95_mk
        check(2) =  1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = -0.05_mk
        check(2) =  1.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))

        check(1) = -0.05_mk  ! right-top
        check(2) =  0.95_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) =  0.95_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))
        check(1) = -0.05_mk
        check(2) = -0.05_mk
        assert_true(found_ghost(p(:,npart+1:mpart),mpart-npart,check))



        deallocate(w)


    end test

    function found_ghost(ghosts,n,cp)
    real(mk),dimension(2,n),intent(in) :: ghosts
    integer                ,intent(in) :: n
    real(mk),dimension(2)  ,intent(in) :: cp
    logical                 :: found_ghost

    integer                 :: i
    integer                 :: found

    found_ghost = .false.

    found = 0
    do i=1,n
        if((abs(ghosts(1,i) - cp(1)).lt.tol).and. &
        &  (abs(ghosts(2,i) - cp(2)).lt.tol)) then
            found = found + 1
        endif
    enddo
    if (found.eq.1) then
        found_ghost = .true.
    endif

    end function found_ghost

    test getsymbc
        ! tests symmetric boundary conditions and ghost get
        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_util_dbg
        integer                         :: npart = 1
        integer                         :: newnpart
        integer                         :: mpart
        real(mk),dimension(:,:),pointer :: p => NULL()
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
        !call ppm_dbg_print_d(topoid,gl,1,1,info,p,npart)

        call ppm_map_part_ghost_get(topoid,p,ndim,npart,0,gl,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)
        
        call ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,2,1,info,p,npart,mpart)


    end test
    
    test get_mixed_bc
        ! tests symmetric and periodic BC conditions and ghost get
        use ppm_module_typedef
        use ppm_module_mktopo
        use ppm_module_topo_check
        use ppm_module_util_dbg
        integer                         :: npart = 256
        integer                         :: snpart = 16
        integer                         :: newnpart
        integer                         :: mpart
        real(mk),dimension(:,:),pointer :: p => NULL()
        real(mk), parameter             :: gl = 0.1_mk
        real(mk)                        :: h
    
        h = len_phys(1)/(snpart)
        allocate(p(ndim,npart))
        k = 0
        do i=1,snpart
            do j=1,snpart
                k = k + 1
                p(1,k) = h/2 + (i-1)*h
                p(2,k) = h/2 + (j-1)*h
            enddo
        enddo
        npart = k
        bcdef(1:2) = ppm_param_bcdef_periodic
        bcdef(3:4) = ppm_param_bcdef_symmetry
        bcdef(5:6) = ppm_param_bcdef_freespace
        !bcdef(1:6) = ppm_param_bcdef_periodic

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

        call ppm_map_part_ghost_get(topoid,p,ndim,npart,1,gl,info)
        call ppm_map_part_send(npart,mpart,info)
        call ppm_map_part_pop(p,ndim,npart,mpart,info)
        assert_equal(mpart-npart,104)
        call ppm_topo_check(topoid,p,npart,ok,info)
        assert_true(ok)
        !call ppm_dbg_print_d(topoid,gl,1,1,info,p,npart,mpart)


    end test


end test_suite
