test_suite ppm_module_mktopo

USE ppm_module_typedef
USE ppm_module_topo
USE ppm_module_map
USE ppm_module_util_dbg
USE ppm_module_data
        
 !------------------------------------------------------------
 ! This test suite tests the topology creation in 2D
 ! All possible decompositions in 2D are tested with the 
 ! following set of particles:
 ! 1) Random in [0,1] with random cutoff [0.05]
 ! 2) Random origin in [-8,8], periodic, the closer to this origin
 !    the smaller the cutoff in this dimension (inhomogenous case)
 ! 3) Same as 2) but proc 0 has almost all particles
 ! 
 !    decomp = ppm_param_decomp_tree
 !    decomp = ppm_param_decomp_pruned_cell
 !    decomp = ppm_param_decomp_bisection
 !    decomp = ppm_param_decomp_xpencil
 !    decomp = ppm_param_decomp_ypencil
 !    decomp = ppm_param_decomp_cuboid
 ! 
 ! The following 4 things are checked for all topologies:
 ! 1) Proc's particles are in its subs
 ! 2) Each particle has the particles to interact with in neighborhood
 ! 3) The subs's sizes are at least as big as the maximum cutoff
 ! 4) Proc has the ghost particles it needs
 ! 
 !------------------------------------------------------------


#ifdef __MPI
    INCLUDE "mpif.h"
#endif

!----------------------
! global variables
!----------------------
integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: skin = 0._mk
integer,parameter               :: ndim = 2
integer                         :: decomp,assig,tolexp
real(mk)                        :: cutoff,min_ghost_req=0.025_mk, max_ghost_req=1.40_mk,so_far,total
integer                         :: info,comm,rank,nproc, topoid
integer                         :: np,mp,newnp,seedsize
real(mk),dimension(:,:),pointer :: xp
real(mk),dimension(:,:),pointer :: ghost_req
real(mk),dimension(:  ),pointer :: min_phys,max_phys, len_phys
integer                         :: i,j,k,sum1,sum2
integer, dimension(6)           :: bcdef
real(mk),dimension(:  ),pointer :: cost
integer                         :: isymm = 0
logical                         :: ok,has_one_way
integer                         :: np_sqrt, np_temp, tot_sum, ix, iy
integer, dimension(:), pointer  :: so_sum
integer, dimension(3)           :: now
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
real(mk), dimension(3)          :: offset

    init

        USE ppm_module_typedef
        USE ppm_module_init
               
        NULLIFY(xp,ghost_req)

#ifdef __MPI
        comm = mpi_comm_world
        CALL mpi_comm_rank(comm,rank,info)
        CALL mpi_comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        tolexp = INT(LOG10(EPSILON(1._mk)))+10
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

        so_far = 0.0_mk
        total = 31.0_mk

    end init


    finalize
        USE ppm_module_finalize
        CALL ppm_finalize(info)
    end finalize


    setup
        IF (rank.eq.0) print*,'Starting ',int(so_far) ,'th test --- ', 100*so_far/total , '% done!'
        so_far = so_far + 1.0_mk
        CALL flush()

        ! Do anything before each test case
        call itime(now)
        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
           seed(i)=now(1)+10+i*i*(rank+1)*now(3)
        enddo
        call random_seed(put=seed)

    end setup
        

    teardown
        DEALLOCATE(xp,ghost_req,stat=info)
        DEALLOCATE(min_phys,max_phys,len_phys)
        DEALLOCATE(seed)
        NULLIFY(xp,ghost_req,min_phys,max_phys,len_phys)
    end teardown

    test decomp_tree_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test decomp_pruned_cell_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_pruned_cell
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test bisection_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test xpencil_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test ypencil_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_ypencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test cuboid_random
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

        np = 1000
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)

        cutoff = 0.05_mk
        has_one_way = .FALSE.
        CALL random_number(xp)

        DO i = 1,np
         ghost_req(1,i) = cutoff
         ghost_req(2,i) = cutoff
        ENDDO
        
        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               cutoff,cost,info)
        Assert_Equal(info,0)

        IF (debug .EQ. 1) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
        ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

    end test

    test decomp_tree_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test decomp_pruned_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_pruned_cell
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test bisection_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test xpencil_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test ypencil_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_ypencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

   test cuboid_example_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test decomp_tree_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test decomp_pruned_cell_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_pruned_cell
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_bisection_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_xpencil_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_ypencil_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_ypencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_cuboid_example_not
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .FALSE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_tree_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         DEALLOCATE(randnb,so_sum)

    end test

    test decomp_pruned_cell_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_pruned_cell
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_bisection_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_xpencil_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_ypencil_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_ypencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_cuboid_example_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        np_sqrt = 16
        np_temp = np_sqrt*np_sqrt
        np = 4*np_temp
        ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
        ALLOCATE(so_sum(ndim),stat=info)
        
        ALLOCATE(randnb(ndim*np),stat=info)
        CALL random_number(randnb)

        has_one_way = .TRUE.
        
        tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
      
        IF (rank .eq. 0) THEN
           offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
           offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
           DO i = 1,ppm_nproc-1
              CALL mpi_send(offset(1),1,ppm_mpi_kind,i,0,comm,info)
              CALL mpi_send(offset(2),1,ppm_mpi_kind,i,0,comm,info)
           ENDDO
        ELSE
              CALL mpi_recv(offset(1),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
              CALL mpi_recv(offset(2),1,ppm_mpi_kind,0,0,comm,MPI_STATUS_IGNORE,info)
        ENDIF

         so_sum(1) = 0.0_mk
         do ix=1,np_sqrt
            so_sum(1) = so_sum(1) + ix
            so_sum(2) = 0.0_mk

            do iy=1,np_sqrt
               so_sum(2) = so_sum(2) + iy

               ! set positions of particles, s.t. lower ix,iy are closer together
               ! including a random distortion            
               xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
      &                                  len_phys(1)/REAL(tot_sum,MK) +  &
      &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
      &                                  - (ix*len_phys(1)/(tot_sum))/2)

               do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                  xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
               enddo

               xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
      &                                  len_phys(2)/REAL(tot_sum,MK) +  &
      &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
      &                                  - (iy*len_phys(2)/(tot_sum))/2)


               do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
               enddo
               do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                  xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
               enddo
               ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
               ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
               ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

            enddo

         enddo

         ! take mirror
         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do ix=1,np_sqrt
            do iy=1,np_sqrt
               xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
               xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
               ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
               ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
            enddo
         enddo

         do i=1,np   
               xp(1,i) = xp(1,i) + offset(1)
               xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
               
               xp(2,i) = xp(2,i) + offset(2)
               xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
         enddo
        
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        DEALLOCATE(randnb,so_sum)

    end test
    
    test decomp_tree_example_single_max
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

         has_one_way = .TRUE.
         
        IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)

         
            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            

            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo
         
            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)   
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF
    end test

    test decomp_tree_example_single_has
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------

         has_one_way = .TRUE.
         
         IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)

         
            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            

            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo
         
            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)   
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_tree
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

         IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test

    test decomp_pruned_cell_example_single
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        has_one_way = .TRUE.
        IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)

            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            

            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk
       

            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_pruned_cell
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test
    
    test decomp_bisection_example_single
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        has_one_way = .TRUE.
        IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)
            
            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            

            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_bisection
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test
    
    test decomp_xpencil_example_single
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        has_one_way = .FALSE.

        IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)
            
            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            
            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_xpencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test
    
    test decomp_ypencil_example_single
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        has_one_way = .TRUE.
        IF (rank .eq. 0) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)

            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            

            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)
        ENDIF
        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_ypencil
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               max_ghost_req,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test
    
    test decomp_cuboid_example_single
        


        ALLOCATE(min_phys(ndim),max_phys(ndim),len_phys(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        bcdef(1:6) = ppm_param_bcdef_periodic

        !----------------------
        ! create particles
        !----------------------
        has_one_way = .TRUE.
        IF (rank .EQ. 0 ) THEN
            np_sqrt = 16
            np_temp = np_sqrt*np_sqrt
            np = 4*np_temp
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            ALLOCATE(so_sum(ndim),stat=info)
            
            ALLOCATE(randnb(ndim*np),stat=info)
            CALL random_number(randnb)

            tot_sum = np_sqrt*(np_sqrt+1)/2 + 1
            
            
            offset(1) = ((-1)**(INT(randnb(1)*10))) * randnb(2) * 7.9999999_mk
            offset(2) = ((-1)**(INT(randnb(3)*10))) * randnb(4) * 7.9999999_mk


            so_sum(1) = 0.0_mk
            do ix=1,np_sqrt
               so_sum(1) = so_sum(1) + ix
               so_sum(2) = 0.0_mk

               do iy=1,np_sqrt
                  so_sum(2) = so_sum(2) + iy

                  ! set positions of particles, s.t. lower ix,iy are closer together
                  ! including a random distortion            
                  xp(1,(ix-1)*np_sqrt + iy) = (min_phys(1)+REAL(so_sum(1),MK)* &
         &                                  len_phys(1)/REAL(tot_sum,MK) +  &
         &                                  (ix*len_phys(1)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy) & 
         &                                  - (ix*len_phys(1)/(tot_sum))/2)

                  do while (xp(1,(ix-1)*np_sqrt + iy) > max_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(1,(ix-1)*np_sqrt + iy) < min_phys(1))
                     xp(1,(ix-1)*np_sqrt + iy) = xp(1,(ix-1)*np_sqrt + iy) + (len_phys(1)/(tot_sum))
                  enddo

                  xp(2,(ix-1)*np_sqrt + iy) = (min_phys(2)+REAL(so_sum(2),MK)* &
         &                                  len_phys(2)/REAL(tot_sum,MK) +  &
         &                                  (iy*len_phys(2)/(tot_sum)) * randnb((ix-1)*np_sqrt + iy + np_temp) & 
         &                                  - (iy*len_phys(2)/(tot_sum))/2)


                  do while (xp(2,(ix-1)*np_sqrt + iy) > max_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) - (len_phys(1)/(tot_sum))
                  enddo
                  do while (xp(2,(ix-1)*np_sqrt + iy) < min_phys(2))
                     xp(2,(ix-1)*np_sqrt + iy) = xp(2,(ix-1)*np_sqrt + iy) + (len_phys(2)/(tot_sum))
                  enddo
                  ! set ghost_req of particles for x and y dim, s.t. lower ix,iy have smaller
         ghost_req(1,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(1),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)
         ghost_req(2,(ix-1)*np_sqrt + iy) = min_ghost_req + REAL(so_sum(2),MK)/REAL(tot_sum,MK)*(max_ghost_req-min_ghost_req)

               enddo

            enddo

            ! take mirror
            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp) = xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do ix=1,np_sqrt
               do iy=1,np_sqrt
                  xp(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(1,(ix-1)*np_sqrt + iy)
                  xp(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = - xp(2,(ix-1)*np_sqrt + iy)
                  ghost_req(1,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(1,(ix-1)*np_sqrt + iy)
                  ghost_req(2,(ix-1)*np_sqrt + iy + np_temp + np_temp + np_temp) = ghost_req(2,(ix-1)*np_sqrt + iy)
               enddo
            enddo

            do i=1,np   
                  xp(1,i) = xp(1,i) + offset(1)
                  xp(1,i) = xp(1,i) - 2*len_phys(1)*int(xp(1,i)/len_phys(1))
                  
                  xp(2,i) = xp(2,i) + offset(2)
                  xp(2,i) = xp(2,i) - 2*len_phys(2)*int(xp(2,i)/len_phys(2))
            enddo
        ELSE
            np = 1
            ALLOCATE(xp(ndim,np),ghost_req(ndim,np),stat=info)
            xp(1,1) = 0.0_mk
            xp(2,1) = 0.0_mk
            ghost_req(1,1) = 0.0_mk
            ghost_req(2,1) = 0.0_mk
            ALLOCATE(randnb(ndim*np),stat=info)
        ENDIF

        min_phys(1:ndim) = -8.0_mk
        max_phys(1:ndim) = 8.0_mk
        len_phys(1:ndim) = max_phys-min_phys

        !----------------------
        ! make topology
        !----------------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal

        topoid = 0

        CALL ppm_mktopo(topoid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
        &               ghost_req,has_one_way,cost,info)
        Assert_Equal(info,0)

        IF (debug .GT. 0) THEN
           CALL ppm_dbg_print(topoid,cutoff,1,1,info,xp,np)
        ENDIF

        !----------------------
        ! global mapping
        !----------------------
        CALL ppm_map_part_global(topoid,xp,np,info)
        CALL ppm_map_part_push(ghost_req,ndim,np,info)
        CALL ppm_map_part_send(np,newnp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,np,newnp,info)
        CALL ppm_map_part_pop(xp,ndim,np,newnp,info)

        !----------------------
        ! ghost mapping
        !----------------------
        CALL ppm_map_part_ghost_get(topoid,xp,ndim,newnp,isymm,info)
        CALL ppm_map_part_push(ghost_req,ndim,newnp,info)
        CALL ppm_map_part_send(newnp,mp,info)
        CALL ppm_map_part_pop(ghost_req,ndim,newnp,mp,info)
        CALL ppm_map_part_pop(xp,ndim,newnp,mp,info)

        !----------------------
        ! make the tests
        !----------------------

        ! 1. Are all particle this proc have in its subs
        CALL ppm_topo_check(topoid,xp,newnp,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check failed'
        ENDIF
        ! 2. check if all particles to interact with are in neighboring box
        !         for each particle determine interacting particles and
        !         check if they are in own box or neighbors
        CALL ppm_topo_check_neigh(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_neigh failed'
        ENDIF
         ! 3. minboxsizes inside, requires all neighbors correct
        !         check for each box if it fulfills ghost req of each particle inside
        !         this function also checks for has_one way case
        CALL ppm_topo_check_minbox(topoid,xp,ghost_req,newnp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_minbox failed'
        ENDIF
        ! 4. check if we have all ghost particles
        !         check for all particles on this proc if it has the ghost particles
        !         it needs. collect all particles from other procs and check
        CALL ppm_map_check_ghosts(topoid,xp,ghost_req,newnp,mp,has_one_way,ok,info)
        Assert_True(ok)
        IF (.NOT. ok) THEN
           write(*,*) '[',rank,'] topo_check_ghosts failed'
        ENDIF

        IF (rank .eq. 0) THEN
            DEALLOCATE(randnb,so_sum)
         ENDIF

    end test

    test xp_null_cuboid
      ! ppm trac ticket #93 (https://ppm.inf.ethz.ch/trac/ticket/93)
      use ppm_module_util_dbg

      integer :: meshid = -1

      deallocate(xp)
      xp => null()

      nm = 32 * nproc

      call ppm_mktopo(topoid, meshid,     &
                      xp, 0,              &
                      decomp, assig,      &
                      min_phys, max_phys, &
                      bcdef, ghostsize,  &
                      cost, Nm, info)

      Assert_Equal(info, 0)

!       call ppm_dbg_print(topoid,0.0_mk,1,1,info)

    end test

end test_suite
